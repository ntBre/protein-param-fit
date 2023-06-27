import click
from collections import defaultdict
import functools
import json
import logging
from multiprocessing import Pool
from openff.qcsubmit.results import (
    OptimizationResultCollection, TorsionDriveResultCollection
)
from openff.qcsubmit.results.filters import SMARTSFilter
from openff.toolkit.typing.engines.smirnoff import ForceField
from pathlib import Path
from qcportal.models import TorsionDriveRecord
from tqdm import tqdm


# Assign parameters for a single molecule
def label_molecule(record_and_molecule, force_field):

    record, molecule = record_and_molecule
    full_labels = force_field.label_molecules(molecule.to_topology())[0]

    molecule_labels = {'smiles': molecule.to_smiles(mapped = True)}

    for parameter_type in full_labels:

        molecule_labels[parameter_type] = {
            str(indices): parameter.id
            for indices, parameter in full_labels[parameter_type].items()
        }

    if isinstance(record, TorsionDriveRecord):
        molecule_labels['dihedral_indices'] = record.keywords.dihedrals

    return molecule_labels


# Assign parameters to molecules in a training dataset
def label_dataset(training_set, force_field, output_path, n_processes):

    dataset_labels = []

    with Pool(n_processes) as pool:

        for molecule_label in tqdm(
            pool.imap(
                functools.partial(label_molecule, force_field = force_field),
                training_set.to_records()
            ),
            total = training_set.n_results
        ):

            dataset_labels.append(molecule_label)

    with open(output_path, 'w') as out_file:
        json.dump(dataset_labels, out_file)

    return dataset_labels


# Select parameters for fitting based on coverage
def select_parameters(
    dataset_labels, force_field, output_path, parameter_types
):

    # Print out coverage information.
    coverage = defaultdict(int)

    for molecule_labels in dataset_labels:

        parameter_ids = set()

        for parameter_type in parameter_types:

            parameter_labels = molecule_labels[parameter_type]

            if 'dihedral_indices' in molecule_labels:

                dihedral_indices = [
                    {*indices[1:3]}
                    for indices in molecule_labels['dihedral_indices']
                ]

                for str_indices, parameter_id in parameter_labels.items():

                
                    indices = {int(i) for i in str_indices.split(',')[1:3]}
                    if indices in dihedral_indices:
                        parameter_ids.add(parameter_id)

            else:

                for parameter_id in parameter_labels.values():
                    parameter_ids.add(parameter_id)

        for parameter_id in parameter_ids:
            coverage[parameter_id] += 1

    # Save out the SMIRKS which should be trained against this set.
    with open(output_path, 'w') as out_file:

        selected_parameters = defaultdict(list)

        for parameter_type in parameter_types:

            for parameter_id, count in coverage.items():

                found_parameters = force_field.get_parameter_handler(
                    parameter_type
                ).get_parameter({'id': parameter_id})

                if (
                    (count < 5 or len(found_parameters) == 0)
                    and 'Protein' not in parameter_id
                ):

                    print(
                        f'Skipping parameter {parameter_id} with {count} '
                        'targets'
                    )
                    continue

                print(f'Fitting parameter {parameter_id} with {count} targets')

                selected_parameters[parameter_type].append(
                    found_parameters[0].smirks
                )

        json.dump(selected_parameters, out_file)


@click.command()
@click.option(
    '-d',
    '--dataset_dir',
    default = Path('..', 'training-datasets'),
    show_default = True,
    type = click.STRING,
    help = 'Directory path to training datasets.',
)
@click.option(
    '-i',
    '--input_ff',
    default = 'initial-force-field.offxml',
    show_default = True,
    type = click.STRING,
    help = 'File path to the initial force field (small molecule + protein '
        'library charges).',
)
@click.option(
    '-n',
    '--n_processes',
    default = 8,
    show_default = True,
    type = click.INT,
    help = 'Number of parallel processes with Pool.',
)
@click.option(
    '-s',
    '--smirks_dir',
    default = 'training-smirks',
    show_default = True,
    type = click.STRING,
    help = 'Directory path to which training SMIRKS will be written.',
)
def main(dataset_dir, input_ff, n_processes, smirks_dir):

    logging.basicConfig(level=logging.INFO)

    initial_force_field = ForceField(input_ff)

    # Filtering out unusual chemistries
    smarts_to_exclude = [
        '[#8+1:1]=[#7:2]',
        # Changed from [#6:2] to [*:2] by Chapin Cavender for protein fits
        #'[#15:1]=[#6:2]',
        '[#15:1]=[*:2]',
        '[#16+1:1]~[*:2]',
        '[*:1]=[#15:2]-[#7:3]~[*:4]',
        '[#17:1]~[#1:2]',
    ]

    optimization_dataset = OptimizationResultCollection.parse_file(
        Path(dataset_dir, 'optimization-elf10-filter-dataset.json')
    )

    optimization_dataset = optimization_dataset.filter(
        SMARTSFilter(smarts_to_exclude = smarts_to_exclude)
    )

    with open(
        Path(smirks_dir, 'optimization-smarts-filter-dataset.json'), 'w'
    ) as json_file:

        json_file.write(optimization_dataset.json())

    optimization_labels = label_dataset(
        optimization_dataset, initial_force_field,
        Path(smirks_dir, 'optimization-labels.json'), n_processes = n_processes
    )

    select_parameters(
        optimization_labels, initial_force_field,
        Path(smirks_dir, 'optimization-bond-smirks.json'),
        parameter_types = ['Bonds'],
    )

    select_parameters(
        optimization_labels, initial_force_field,
        Path(smirks_dir, 'optimization-angle-smirks.json'),
        parameter_types = ['Angles'],
    )

    # Concatenate small molecule and protein TorsionDrive datasets together
    torsiondrive_dataset = TorsionDriveResultCollection.parse_file(
        Path(
            dataset_dir, 'small-molecule-torsiondrive-elf10-filter-dataset.json'
        )
    )

    torsiondrive_dataset = torsiondrive_dataset.filter(
        SMARTSFilter(smarts_to_exclude = smarts_to_exclude)
    )

    with open(
        Path(
            smirks_dir, 'small-molecule-torsiondrive-smarts-filter-dataset.json'
        ),
        'w'
    ) as json_file:

        json_file.write(torsiondrive_dataset.json())

    for weight in [1, 2, 3, 4, 6]:

        protein_torsiondrive_dataset = TorsionDriveResultCollection.parse_file(
            Path(
                dataset_dir,
                f'protein-weight{weight}-torsiondrive-elf10-filter-dataset.json'
            )
        )

        protein_torsiondrive_dataset = protein_torsiondrive_dataset.filter(
            SMARTSFilter(smarts_to_exclude = smarts_to_exclude)
        )

        with open(
            Path(
                smirks_dir,
                f'protein-weight{weight}-torsiondrive-smarts-filter-dataset.json'
            ),
            'w'
        ) as json_file:

            json_file.write(protein_torsiondrive_dataset.json())

        for address in torsiondrive_dataset.entries.keys():

            torsiondrive_dataset.entries[address].extend(
                protein_torsiondrive_dataset.entries[address]
            )

    torsiondrive_labels = label_dataset(
        torsiondrive_dataset, initial_force_field,
        Path(smirks_dir, 'torsiondrive-labels.json'), n_processes = n_processes
    )

    # Select parameters for concatenated TorsionDrive datasets
    select_parameters(
        torsiondrive_labels, initial_force_field,
        Path(smirks_dir, 'torsiondrive-torsion-smirks.json'),
        parameter_types = ['ProperTorsions'],
    )


if __name__ == "__main__":
    main()

