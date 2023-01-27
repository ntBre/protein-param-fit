from aa_library_charges import INDEX_TO_AA
from mean_charges import AA_WEIGHTS
import numpy
from openff.toolkit.topology import Molecule
from openmm import unit
import os
from polymetrizer.tests.smiles import ACE, NME


def main():

    ace_smiles = '[C:1](=[O:2])[C:3]([H:4])([H:5])[H:6]'
    nme_smiles = '[N:1]([H:2])[C:3]([H:4])([H:5])[H:6]'

    ace_name_map = ['C1', 'O2', 'C3', 'H4', 'H5', 'H6']
    nme_name_map = ['N1', 'H2', 'C3', 'H4', 'H5', 'H6']

    charge_dir = 'aa_library_charges'

    sum_ace_charges = numpy.zeros(6)
    sum_nme_charges = numpy.zeros(6)
    sum_weights = 0

    with open(
        os.path.join('mean_charges', f'ace_charges.dat'), 'w'
    ) as ace_out_file:

        with open(
            os.path.join('mean_charges', f'nme_charges.dat'), 'w'
        ) as nme_out_file:

            for charge_dir_file in os.listdir(charge_dir):

                if (
                    'Ace' not in charge_dir_file or 'Nme' not in charge_dir_file
                    or not charge_dir_file.endswith('.smi')
                ):

                    continue

                with open(
                    os.path.join(charge_dir, charge_dir_file), 'r'
                ) as smiles_file:

                    mol_name, mol_smiles = smiles_file.readline().split()

                # Split molecule name at hyphens to get list of residues
                residues = mol_name.split('-')[1:-1]
                residue_str = ' '.join(residues)

                mol = Molecule.from_mapped_smiles(
                    mol_smiles, allow_undefined_stereo = True
                )
                ace_atoms = list(
                    mol.chemical_environment_matches(ace_smiles, unique=True)[0]
                )
                nme_atoms = list(
                    mol.chemical_environment_matches(nme_smiles, unique=True)[0]
                )

                mol_charges = numpy.loadtxt(
                    os.path.join(charge_dir, f'{mol_name}.dat')
                )

                if mol_name == 'Arg_Ace-Val-Arg-Arg-Cys-Val-Nme':
                    weight = 1.0
                else:
                    weight = numpy.prod([AA_WEIGHTS[res] for res in residues])
                sum_ace_charges += weight * mol_charges[ace_atoms]
                sum_nme_charges += weight * mol_charges[nme_atoms]
                sum_weights += weight

                ace_charges = [f'{q:11.8f}' for q in mol_charges[ace_atoms]]
                ace_charge_str = ' '.join(ace_charges)
                nme_charges = [f'{q:11.8f}' for q in mol_charges[nme_atoms]]
                nme_charge_str = ' '.join(nme_charges)

                ace_out_file.write(residue_str + ' ' + ace_charge_str + '\n')
                nme_out_file.write(residue_str + ' ' + nme_charge_str + '\n')

    # Get mean charge over all oligomers
    ace_mean_charges = sum_ace_charges / sum_weights
    nme_mean_charges = sum_nme_charges / sum_weights

    # Distribute non-integer charge
    ace_total_charge = numpy.sum(ace_mean_charges)
    if ace_total_charge != 0.0:
        ace_mean_charges -= ace_total_charge / 6
    nme_total_charge = numpy.sum(nme_mean_charges)
    if nme_total_charge != 0.0:
        nme_mean_charges -= nme_total_charge / 6

    with open(
        os.path.join('mean_charges', 'mean_charges.dat'), 'a'
    ) as out_file:

        for i in range(6):

            out_file.write(
                f'Ace {ace_name_map[i]:3s} {ace_mean_charges[i]:11.8f}\n'
            )

        for i in range(6):

            out_file.write(
                f'Nme {nme_name_map[i]:3s} {nme_mean_charges[i]:11.8f}\n'
            )


if __name__ == "__main__":
    main()

