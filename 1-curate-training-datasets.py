from collections import defaultdict
import click
import copy
import functools
import json
import logging
from multiprocessing import Pool
from openff.qcsubmit.results import (
    OptimizationResultCollection, TorsionDriveResultCollection
)
from openff.qcsubmit.results.filters import (
    ConformerRMSDFilter, ConnectivityFilter, ElementFilter, HydrogenBondFilter,
    RecordStatusFilter, ResultRecordFilter, SMARTSFilter, SMILESFilter
)
from openff.toolkit.topology import Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils import UndefinedStereochemistryError
from pathlib import Path
from qcportal import FractalClient
from qcportal.models import TorsionDriveRecord
from qcportal.models.records import RecordStatusEnum
import random
from tempfile import NamedTemporaryFile
from tqdm import tqdm


class UndefinedStereoFilter(ResultRecordFilter):

    def _filter_function(self, result, record, molecule) -> bool:

        has_stereochemistry = True

        molecule = copy.deepcopy(molecule)
        molecule._conformers = [molecule.conformers[0]]

        try:
            with NamedTemporaryFile(suffix = '.sdf') as tmp_file:

                molecule.to_file(tmp_file.name, 'SDF')
                molecule.from_file(tmp_file.name)

        except UndefinedStereochemistryError:
            has_stereochemistry = False

        return has_stereochemistry


@click.command()
@click.option(
    '-d',
    '--dataset_dir',
    default = 'training-datasets',
    show_default = True,
    type = click.STRING,
    help = 'Directory path to which training datasets will be written.',
)
def main(dataset_dir):

    logging.getLogger('openff').setLevel(logging.ERROR)

    Path(dataset_dir).mkdir(parents = True, exist_ok = True)

    default_filters = [
        RecordStatusFilter(status = RecordStatusEnum.complete),
        ConnectivityFilter(tolerance = 1.2),
        UndefinedStereoFilter(),
    ]

    # Pull down the TorsionDrive and Optimization training datasets and filter
    # out records that are not complete or that contain intramolecular hydrogen
    # bonds
    client = FractalClient()

    # Small molecule TorsionDrives
    torsiondrive_dataset = TorsionDriveResultCollection.from_server(
        client = client,
        datasets = [
            'OpenFF Gen 2 Torsion Set 1 Roche 2',
            'OpenFF Gen 2 Torsion Set 2 Coverage 2',
            'OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2',
            'OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2',
            'OpenFF Gen 2 Torsion Set 5 Bayer 2',
            'OpenFF Gen 2 Torsion Set 6 supplemental 2',
        ],
        spec_name = 'default',
    )

    # Drop record ids with inconsistent optimization histories or that cause
    # failures in ForceBalance
    torsiondrive_dataset.entries[client.address] = [
        entry for entry in torsiondrive_dataset.entries[client.address]
        if entry.record_id not in [
            '6098580',
            '2703504',
            '2703505',
            '18045478',
            # SMIRNOFF Coverage torsions set inconsistent IDs
            '2703253',
            '2703343',
            '2703386',
            '2703439',
            '2703449',
            '2703545',
            '2703546',
            '2703616',
            # from Gen3 set probably
            '35045000',
            # ForceBalance failures added by Chapin Cavender for protein fits
            '2703525',
            '2703603',
            '18045477',
            '18536106',
        ]
    ]

    torsiondrive_dataset = torsiondrive_dataset.filter(
        *default_filters,
        HydrogenBondFilter(method = 'baker-hubbard'),
        ElementFilter(
            # The elements supported by SMIRNOFF. Excluding iodine here since we
            # don't have iodine torsions and any record with iodine is tainted
            # for the above datasets because of the auxiliary basis set issue.
            allowed_elements=['H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br']
        )
    )

    with open(
        Path(dataset_dir, 'small-molecule-torsiondrive-dataset.json'), 'w'
    ) as json_file:

        json_file.write(torsiondrive_dataset.json())

    # Protein TorsionDrives. Don't filter incomplete records or records with
    # hydrogen bonds for these.
    #protein_torsiondrive_dataset = TorsionDriveResultCollection.from_server(
    #    client=FractalClient(),
    #    datasets = [
    #        'OpenFF Protein Dipeptide 2-D TorsionDrive v2.1',
    #        'OpenFF Protein Capped 1-mer Sidechains v1.2',
    #    ],
    #    spec_name='default',
    #)
    # Hack to avoid filtering incomplete results from
    # TorsionDriveResultCollection.from_datasets()
    from openff.qcsubmit.results import TorsionDriveResult
    protein_datasets = [
        client.get_collection('TorsionDriveDataset', dataset_name)
        for dataset_name in [
            'OpenFF Protein Dipeptide 2-D TorsionDrive v2.1',
            'OpenFF Protein Capped 1-mer Sidechains v1.2',
        ]
    ]
    result_records = defaultdict(dict)
    for dataset in protein_datasets:
        query = dataset.query("default")
        result_records[dataset.client.address].update(
            {
                query[entry.name].id: TorsionDriveResult(
                    record_id=query[entry.name].id,
                    cmiles=entry.attributes[
                        "canonical_isomeric_explicit_hydrogen_mapped_smiles"
                    ],
                    inchi_key=entry.attributes.get("fixed_hydrogen_inchi_key")
                    or Molecule.from_mapped_smiles(
                        entry.attributes[
                            "canonical_isomeric_explicit_hydrogen_mapped_smiles"
                        ],
                        allow_undefined_stereo=True,
                    ).to_inchikey(fixed_hydrogens=True),
                )
                for entry in dataset.data.records.values()
                if entry.name in query
            }
        )
    protein_torsiondrive_dataset = TorsionDriveResultCollection(
        entries={
            address: [*entries.values()]
            for address, entries in result_records.items()
        }
    )

    protein_torsiondrive_dataset = protein_torsiondrive_dataset.filter(
        ConnectivityFilter(tolerance = 1.2),
        UndefinedStereoFilter(),
        ElementFilter(
            # The elements supported by SMIRNOFF. Excluding iodine here since we
            # don't have iodine torsions and any record with iodine is tainted
            # for the above datasets because of the auxiliary basis set issue.
            allowed_elements=['H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br']
        )
    )

    # Each dihedral scan for each amino acid should be weighted equally.
    # Separate TorsionDrive records by weight to be applied based on number
    # of protomers/tautomers and number of conformations of constrained
    # non-driven dihedrals.
    protein_torsiondrive_datasets = dict()
    protein_torsiondrive_datasets['weight1'] = TorsionDriveResultCollection(
        entries = {
            address: [
                entry for entry in entries
                if entry.record_id in [
                    '99478323', # Backbone ALA-rotamer-1
                    '99478356', # Backbone ARG-rotamer-1
                    '99478358', # Backbone ASN-rotamer-1
                    '99478380', # Backbone GLN-rotamer-1
                    '99478419', # Backbone GLY-rotamer-1
                    '99478504', # Backbone ILE-rotamer-1
                    '99478505', # Backbone LEU-rotamer-1
                    '99478599', # Backbone MET-rotamer-1
                    '99478630', # Backbone PHE-rotamer-1
                    '99478688', # Backbone PRO-rotamer-1
                    '99478689', # Backbone SER-rotamer-1
                    '99478690', # Backbone THR-rotamer-1
                    '99478719', # Backbone TRP-rotamer-1
                    '99478752', # Backbone TYR-rotamer-1
                    '99478784', # Backbone VAL-rotamer-1
                ]
            ]
            for address, entries in protein_torsiondrive_dataset.entries.items()
        }
    )

    protein_torsiondrive_datasets['weight2'] = TorsionDriveResultCollection(
        entries = {
            address: [
                entry for entry in entries
                if entry.record_id in [
                    '99478357', # Backbone ASH-rotamer-1
                    '99478359', # Backbone ASP-rotamer-1
                    '99478377', # Backbone CYS-rotamer-1
                    '99478378', # Backbone CYX-rotamer-1
                    '99478379', # Backbone GLH-rotamer-1
                    '99478418', # Backbone GLU-rotamer-1
                    '99478597', # Backbone LYN-rotamer-1
                    '99478598', # Backbone LYS-rotamer-1
                    '103670809', # Sidechain ARG-alpha
                    '103670810', # Sidechain ARG-beta
                    '103670813', # Sidechain ASN-alpha
                    '103670814', # Sidechain ASN-beta
                    '103670823', # Sidechain GLN-alpha
                    '103670824', # Sidechain GLN-beta
                    '103670833', # Sidechain ILE-alpha
                    '103670834', # Sidechain ILE-beta
                    '103670835', # Sidechain LEU-alpha
                    '103670836', # Sidechain LEU-beta
                    '103670841', # Sidechain MET-alpha
                    '103670842', # Sidechain MET-beta
                    '103670843', # Sidechain PHE-alpha
                    '103670844', # Sidechain PHE-beta
                    '103670845', # Sidechain SER-alpha
                    '103670846', # Sidechain SER-beta
                    '103670847', # Sidechain THR-alpha
                    '103670848', # Sidechain THR-beta
                    '103670849', # Sidechain TRP-alpha
                    '103670850', # Sidechain TRP-beta
                    '103670851', # Sidechain TYR-alpha
                    '103670852', # Sidechain TYR-beta
                    '103670853', # Sidechain VAL-alpha
                    '103670854', # Sidechain VAL-beta
                ]
            ]
            for address, entries in protein_torsiondrive_dataset.entries.items()
        }
    )

    protein_torsiondrive_datasets['weight3'] = TorsionDriveResultCollection(
        entries = {
            address: [
                entry for entry in entries
                if entry.record_id in [
                    '99478420', # Backbone HID-rotamer-1
                    '99478421', # Backbone HIE-rotamer-1
                    '99478450', # Backbone HIP-rotamer-1
                ]
            ]
            for address, entries in protein_torsiondrive_dataset.entries.items()
        }
    )

    protein_torsiondrive_datasets['weight4'] = TorsionDriveResultCollection(
        entries = {
            address: [
                entry for entry in entries
                if entry.record_id in [
                    '103670811', # Sidechain ASH-alpha
                    '103670812', # Sidechain ASH-beta
                    '103670815', # Sidechain ASP-alpha
                    '103670816', # Sidechain ASP-beta
                    '103670817', # Sidechain CYS-alpha
                    '103670818', # Sidechain CYS-beta
                    '103670819', # Sidechain CYX-alpha
                    '103670820', # Sidechain CYX-beta
                    '103670821', # Sidechain GLH-alpha
                    '103670822', # Sidechain GLH-beta
                    '103670825', # Sidechain GLU-alpha
                    '103670826', # Sidechain GLU-beta
                    '103670837', # Sidechain LYN-alpha
                    '103670838', # Sidechain LYN-beta
                    '103670839', # Sidechain LYS-alpha
                    '103670840', # Sidechain LYS-beta
                ]
            ]
            for address, entries in protein_torsiondrive_dataset.entries.items()
        }
    )

    protein_torsiondrive_datasets['weight6'] = TorsionDriveResultCollection(
        entries = {
            address: [
                entry for entry in entries
                if entry.record_id in [
                    '103670827', # Sidechain HID-alpha
                    '103670828', # Sidechain HID-beta
                    '103670829', # Sidechain HIE-alpha
                    '103670830', # Sidechain HIE-beta
                    '103670831', # Sidechain HIP-alpha
                    '103670832', # Sidechain HIP-beta
                ]
            ]
            for address, entries in protein_torsiondrive_dataset.entries.items()
        }
    )

    for weight, dataset in protein_torsiondrive_datasets.items():

        with open(
            Path(dataset_dir, f'protein-{weight}-torsiondrive-dataset.json'),
            'w'
        ) as json_file:

            json_file.write(dataset.json())

    # Get optimization datasets without iodine, then get datasets with iodine,
    # then combine and filter
    optimization_dataset = OptimizationResultCollection.from_server(
        client = FractalClient(),
        datasets = [
            # Small molecule optimization datasets
            'OpenFF Gen 2 Opt Set 1 Roche',
            'OpenFF Gen 2 Opt Set 2 Coverage',
            'OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy',
            'OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy',
            'OpenFF Gen 2 Opt Set 5 Bayer',
            # Protein optimization datasets
            'OpenFF Protein Capped 1-mers 3-mers Optimization Dataset v1.0',
        ],
        spec_name = 'default',
    )

    optimization_dataset = optimization_dataset.filter(
        ElementFilter(
            # The elements supported by SMIRNOFF. Excluding iodine here since we
            # don't have iodine torsions and any record with iodine is tainted
            # for the above datasets because of the auxiliary basis set issue.
            # New sets added below in optimization_dataset_iodine have iodine
            # containing molecules that are safe.
            allowed_elements=['H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br']
        ),
    )

    iodine_optimization_dataset = OptimizationResultCollection.from_server(
        client = FractalClient(),
        datasets = [
            'OpenFF Gen2 Optimization Dataset Protomers v1.0',
            'OpenFF Iodine Chemistry Optimization Dataset v1.0',
        ],
        spec_name = 'default',
    )

    iodine_optimization_dataset = iodine_optimization_dataset.filter(
        ElementFilter(
            # The elements supported by SMIRNOFF, including iodine
            allowed_elements=[
                'H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I',
            ]
        ),
    )

    optimization_dataset.entries[client.address].extend(
        iodine_optimization_dataset.entries[client.address]
    )

    optimization_dataset.entries[client.address] = [
        entry for entry in optimization_dataset.entries[client.address]
        if entry.record_id not in [
            '2002949',
            '2002950',
            # ForceBalance failures added by Chapin Cavender for protein fits
            '18434120',
        ]
    ]

    optimization_dataset = optimization_dataset.filter(
        *default_filters,
        ConformerRMSDFilter(max_conformers = 10),
    )

    with open(Path(dataset_dir, 'optimization-dataset.json'), 'w') as json_file:
        json_file.write(optimization_dataset.json())


if __name__ == "__main__":
    main()

