import click
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
import openmm
from openmm import unit
from pathlib import Path
from residue_smirks import RESIDUE_SMIRKS


@click.command()
@click.option(
    '-f',
    '--force_field',
    default = 'protein-library-charges.offxml',
    show_default = True,
    type = click.STRING,
    help = 'File path to the force field containing protein library charges.',
)
@click.option(
    '-l',
    '--library_charge_file',
    default = Path('mean-charges', 'mean-charges.dat'),
    show_default = True,
    type = click.STRING,
    help = 'File path to the protein library charges.',
)
@click.option(
    '-s',
    '--smiles_file',
    default = Path(
        '..', '..', 'software', 'qca-dataset-submission', 'submissions',
        '2022-05-30-OpenFF-Protein-Capped-1-mers-3-mers-Optimization',
        'capped_1-mers_3-mers.smi'
    ),
    show_default = True,
    type = click.STRING,
    help = 'File path to the list of peptide SMILES strings.',
)
def main(force_field, library_charge_file, smiles_file):

    ff = ForceField('openff-2.0.0.offxml', force_field)

    # Get list of SMILES for peptides to test
    peptide_smiles = []
    with open(smiles_file, 'r') as smi:
        for line in smi:
            peptide_smiles.append(line.strip())

    # List of residues
    residues = [
        'ALA', 'ARG', 'ASH', 'ASN', 'ASP', 'CYS', 'CYX', 'GLH', 'GLN', 'GLU',
        'GLY', 'HID', 'HIE', 'HIP', 'ILE', 'LEU', 'LYN', 'LYS', 'MET', 'PHE',
        'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    ]

    # Read library charges
    library_charges = dict()
    with open(library_charge_file, 'r') as charge_file:
        for line in charge_file:

            residue, atom, charge = line.split()
            key = residue.upper()

            if key in residues or key in ['ACE', 'NME']:

                if key not in library_charges:
                    library_charges[key] = []

                library_charges[key].append(float(charge))

    # Loop over capped 1-mers
    for residue_idx, residue in enumerate(residues):

        offmol = Molecule.from_smiles(peptide_smiles[residue_idx])
        omm_system = ff.create_openmm_system(offmol.to_topology())
        nb_force = [
            force for force in omm_system.getForces()
            if isinstance(force, openmm.NonbondedForce)
        ][0]

        for residue_name in ['ACE', residue, 'NME']:

            residue_atoms = offmol.chemical_environment_matches(
                RESIDUE_SMIRKS[residue_name], unique = True
            )[0]

            for charge_idx, atom_idx in enumerate(residue_atoms):

                library_charge = library_charges[residue_name][charge_idx]
                openmm_charge = nb_force.getParticleParameters(
                    atom_idx)[0].value_in_unit(unit.elementary_charge)

                print(
                    f'ACE-{residue}-NME {residue_name} {charge_idx:2d} '
                    f'{offmol.atoms[atom_idx].element.symbol} '
                    f'{library_charge:11.8f} {openmm_charge:11.8f}'
                )


if __name__ == "__main__":
    main()

