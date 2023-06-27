from aa_library_charges import INDEX_TO_AA
from mean_charges import AA_WEIGHTS
import numpy
from openff.toolkit.topology import Molecule
from openmm import unit
from pathlib import Path
from polymetrizer.tests.smiles import ACE, NME


def main():

    ace_smiles = '[C:1](=[O:2])[C:3]([H:4])([H:5])[H:6]'
    nh2_smiles = '[NX3][CX4][CX3](=O)[N:1]([H:2])[H:3]'
    nme_smiles = '[N:1]([H:2])[C:3]([H:4])([H:5])[H:6]'

    ace_name_map = ['C1', 'O2', 'C3', 'H4', 'H5', 'H6']
    nh2_name_map = ['N1', 'H2', 'H3']
    nme_name_map = ['N1', 'H2', 'C3', 'H4', 'H5', 'H6']

    charge_dir = Path('aa_library_charges')

    sum_ace_charges = numpy.zeros(6)
    sum_nh2_charges = numpy.zeros(3)
    sum_nme_charges = numpy.zeros(6)

    sum_nh2_weights = 0
    sum_nme_weights = 0

    with open(
        Path('mean_charges', f'nh2-ace-charges.dat'), 'w'
    ) as ace_out_file:

        with open(
            Path('mean_charges', f'nh2-nh2-charges.dat'), 'w'
        ) as nh2_out_file:

            with open(
                Path('mean_charges', f'nh2-nme-charges.dat'), 'w'
            ) as nme_out_file:

                for charge_dir_file in charge_dir.iterdir():

                    if (
                        not charge_dir_file.name.endswith('.smi') or not (
                            charge_dir_file.name.startswith('Nme')
                            or charge_dir_file.name.startswith('Nh2')
                        )
                    ):

                        continue

                    with open(charge_dir_file, 'r') as smiles_file:
                        mol_name, mol_smiles = smiles_file.readline().split()

                    # Split molecule name at hyphens to get list of residues
                    residues = mol_name.split('-')[1:-1]
                    residue_str = ' '.join(residues)

                    mol = Molecule.from_mapped_smiles(
                        mol_smiles, allow_undefined_stereo = True
                    )
                    ace_atoms = mol.chemical_environment_matches(
                        ace_smiles, unique = True
                    )
                    nh2_atoms = mol.chemical_environment_matches(
                        nh2_smiles, unique = True
                    )
                    nme_atoms = mol.chemical_environment_matches(
                        nme_smiles, unique = True
                    )

                    mol_charges = numpy.loadtxt(
                        Path(charge_dir, f'{mol_name}.dat')
                    )

                    weight = numpy.prod([AA_WEIGHTS[res] for res in residues])

                    ace_charges = mol_charges[list(ace_atoms[0])]
                    sum_ace_charges += weight * ace_charges
                    ace_charge_str = ' '.join(
                        [f'{q:11.8f}' for q in ace_charges]
                    )
                    ace_out_file.write(f'{residue_str} {ace_charge_str}\n')

                    if charge_dir_file.name.startswith('Nh2'):

                        nh2_charges = mol_charges[list(nh2_atoms[0])]
                        sum_nh2_charges += weight * nh2_charges
                        sum_nh2_weights += weight
                        nh2_charge_str = ' '.join(
                            [f'{q:11.8f}' for q in nh2_charges]
                        )
                        nh2_out_file.write(f'{residue_str} {nh2_charge_str}\n')

                    else:

                        nme_charges = mol_charges[list(nme_atoms[0])]
                        sum_nme_charges += weight * nme_charges
                        sum_nme_weights += weight
                        nme_charge_str = ' '.join(
                            [f'{q:11.8f}' for q in nme_charges]
                        )
                        nme_out_file.write(f'{residue_str} {nme_charge_str}\n')

    # Get mean charge over all oligomers
    ace_mean_charges = sum_ace_charges / (sum_nh2_weights + sum_nme_weights)
    nh2_mean_charges = sum_nh2_charges / sum_nh2_weights
    nme_mean_charges = sum_nme_charges / sum_nme_weights

    # Distribute non-integer charge
    ace_total_charge = numpy.sum(ace_mean_charges)
    if ace_total_charge != 0.0:
        ace_mean_charges -= ace_total_charge / len(ace_mean_charges)

    nh2_total_charge = numpy.sum(nh2_mean_charges)
    if nh2_total_charge != 0.0:
        nh2_mean_charges -= nh2_total_charge / len(nh2_mean_charges)

    nme_total_charge = numpy.sum(nme_mean_charges)
    if nme_total_charge != 0.0:
        nme_mean_charges -= nme_total_charge / len(nme_mean_charges)

    with open(Path('mean_charges', 'nh2-mean-charges.dat'), 'a') as out_file:

        for i in range(len(ace_mean_charges)):

            out_file.write(
                f'Ace {ace_name_map[i]:3s} {ace_mean_charges[i]:11.8f}\n'
            )

        for i in range(len(nh2_mean_charges)):

            out_file.write(
                f'Nh2 {nh2_name_map[i]:3s} {nh2_mean_charges[i]:11.8f}\n'
            )

        for i in range(len(nme_mean_charges)):

            out_file.write(
                f'Nme {nme_name_map[i]:3s} {nme_mean_charges[i]:11.8f}\n'
            )


if __name__ == "__main__":
    main()

