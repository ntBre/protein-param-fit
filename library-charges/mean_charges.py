import click
import numpy
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from openmm import unit
import os
from polymetrizer import Monomer
from polymetrizer.tests.smiles import AMINO_ACIDS_ALL_STATES, ACE, NME
from residue_smirks import RESIDUE_SMIRKS
from terminal_library_charges import (
    neutral_c_term, neutral_n_term, terminal_smiles
)


AA_WEIGHTS = {
    'Ala': 1.0,
    'Arg': 1.0,
    'Ash': 0.5,
    'Asn': 1.0,
    'Asp': 0.5,
    'Cys': 0.5,
    'Cyx': 0.5,
    'Glh': 0.5,
    'Gln': 1.0,
    'Glu': 0.5,
    'Gly': 1.0,
    'Hid': 1.0 / 3,
    'Hie': 1.0 / 3,
    'Hip': 1.0 / 3,
    'Ile': 1.0,
    'Leu': 1.0,
    'Lyn': 0.5,
    'Lys': 0.5,
    'Met': 1.0,
    'Phe': 1.0,
    'Pro': 1.0,
    'Ser': 1.0,
    'Thr': 1.0,
    'Trp': 1.0,
    'Tyr': 1.0,
    'Val': 1.0,
}

RESIDUE_WEIGHTS = dict()
for residue in AA_WEIGHTS:

    RESIDUE_WEIGHTS[residue] = AA_WEIGHTS[residue]

    for terminal_prefix in ['N0', 'N', 'C0', 'C']:
        RESIDUE_WEIGHTS[terminal_prefix + residue] = AA_WEIGHTS[residue]

RESIDUE_WEIGHTS['Ace'] = 1.0
RESIDUE_WEIGHTS['Nme'] = 1.0


# Compute mean charges for a residue
def compute_mean_charges(
    residue_name, residue_offmol, residue_index, charge_dir, mean_dir
):

    # Get atom map to residue SMIRKS with Amber atom ordering
    residue_smirks = RESIDUE_SMIRKS[residue_name.upper()]
    smirks_matches = residue_offmol.chemical_environment_matches(
        residue_smirks, unique = True
    )

    if len(smirks_matches) == 0:

        raise ValueError(
            f'Molecule {residue_offmol.to_smiles(mapped = True)} does not '
            f'match {residue_name} residue SMIRKS {residue_smirks}'
        )

    # Get atom indices of residue from polymetrizer SMILES. Since the oligomer
    # was built from the residue of interest first, that residue will have the
    # lowest indices.
    atom_idx = 0
    atom_map = dict()

    for i, atom in enumerate(residue_offmol.atoms):
        if atom.atomic_number:

            atom_idx += 1
            atom_map[i] = atom_idx

    num_atoms = atom_idx
    residue_offmol.properties['atom_map'] = atom_map

    sum_charges = numpy.zeros(num_atoms)
    sum_weights = 0

    cap = residue_name in ['Ace', 'Nme']

    with open(
        os.path.join(mean_dir, f'{residue_name.lower()}-charges.dat'), 'w'
    ) as out_file:

        for charge_dir_file in os.listdir(charge_dir):

            # Loop over files containing molecule SMILES in charge_dir. For
            # caps, loop over files containing both Ace and Nme. For other
            # residues, loop over files starting with the residue name.
            if (
                not charge_dir_file.endswith('.smi') or not (
                    charge_dir_file.startswith(f'{residue_name}_') or (
                        cap and 'Ace' in charge_dir_file
                        and 'Nme' in charge_dir_file
                    )
                )
            ):

                continue

            with open(
                os.path.join(charge_dir, charge_dir_file), 'r'
            ) as smiles_file:

                mol_name, mol_smiles = smiles_file.readline().split()

            # Split molecule name at hyphens to get list of residues
            residues = mol_name.split('_')[-1].split('-')
            residue_str = ' '.join(residues)

            if residues[residue_index] != residue_name:

                raise ValueError(
                    f'Residue {residue_index} is {residues[residue_index]} '
                    f'instead of {residue_name}'
                )

            mol_charges = numpy.loadtxt(
                os.path.join(charge_dir, f'{mol_name}.dat')
            )

            # Weight for Ace-Val-Arg-Arg-Cys-Val-Nme is 1.0 because ELF10 failed
            # for Ace-Val-Arg-Arg-Cyx-Val-Nme
            if mol_name == 'Arg_Ace-Val-Arg-Arg-Cys-Val-Nme':
                weight = 1.0
            else:
                weight = numpy.prod([RESIDUE_WEIGHTS[res] for res in residues])

            if cap:

                offmol = Molecule.from_mapped_smiles(
                    mol_smiles, allow_undefined_stereo = True
                )
                cap_matches = offmol.chemical_environment_matches(
                    residue_smirks, unique = True
                )
                residue_charges = mol_charges[list(cap_matches[0])]

            else:
                residue_charges = mol_charges[:num_atoms]

            sum_charges += weight * residue_charges
            sum_weights += weight

            charge_str = ' '.join([f'{q:11.8f}' for q in residue_charges])
            out_file.write(residue_str + ' ' + charge_str + '\n')

    # Get mean charge over all oligomers
    mean_charges = sum_charges / sum_weights

    # Distribute non-integer charge
    expected_charge = residue_offmol.total_charge.value_in_unit(
        unit.elementary_charge
    )
    total_charge = numpy.sum(mean_charges)
    charge_diff = expected_charge - total_charge

    print(
        f'Residue {residue_name} Expected charge {expected_charge:.2f} e '
        f'Total charge {total_charge:.2f} e'
    )

    if charge_diff != 0:
        mean_charges += charge_diff / num_atoms

    # Write mean charges to file
    with open(os.path.join(mean_dir, 'mean-charges.dat'), 'a') as out_file:
        for smirks_idx, atom_idx in enumerate(smirks_matches[0]):

            element_symbol = residue_offmol.atoms[atom_idx].element.symbol

            if cap:
                charge_idx = smirks_idx
            else:
                charge_idx = residue_offmol.properties['atom_map'][atom_idx] - 1

            out_file.write(
                f'{residue_name} {element_symbol}{smirks_idx + 1:<2d} '
                f'{mean_charges[charge_idx]:11.8f}\n'
            )


@click.command()
@click.option(
    '-c',
    '--charge_dir',
    default = 'aa_library_charges',
    show_default = True,
    type = click.STRING,
    help = 'Directory path containing ELF10 charges for molecules.',
)
@click.option(
    '-m',
    '--mean_dir',
    default = 'mean-charges',
    show_default = True,
    type = click.STRING,
    help = 'Directory path to which mean charges will be written.',
)
def main(charge_dir, mean_dir):

    # Compute mean charges in order of SMIRNOFF library charge parameters, i.e.
    # general to specific. Start with caps: Ace and Nme.
    ace_residue = Molecule.from_smiles(ACE, allow_undefined_stereo = True)
    compute_mean_charges('Ace', ace_residue, 0, charge_dir, mean_dir)

    nme_residue = Molecule.from_smiles(NME, allow_undefined_stereo = True)
    compute_mean_charges('Nme', nme_residue, 6, charge_dir, mean_dir)

    # Mostly alphabetical, but make sure 1) CYX is before CYS and 2) HIP is
    # after HID and HIE
    aa_order = [
        'Ala', 'Arg', 'Ash', 'Asn', 'Asp', 'Cyx', 'Cys', 'Glh', 'Gln', 'Glu',
        'Gly', 'Hid', 'Hie', 'Hip', 'Ile', 'Leu', 'Lyn', 'Lys', 'Met', 'Phe',
        'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'
    ]

    # Main chain residues, i.e. Y in Ace-Val-X-Y-Z-Val-Nme
    for aa_name in aa_order:

        central_residue = Molecule.from_smiles(
            AMINO_ACIDS_ALL_STATES[aa_name], allow_undefined_stereo = True
        )

        compute_mean_charges(aa_name, central_residue, 3, charge_dir, mean_dir)

    # N terminal residues, i.e. X in X-Y-Val-Nme. Neutral must come before
    # charged.
    for aa_name in aa_order:

        terminal_name = f'N0{aa_name}'
        terminal_monomer = Monomer.from_smiles(
            AMINO_ACIDS_ALL_STATES[aa_name], name = terminal_name
        )
        terminal_residue = terminal_monomer.substitute(
            neutral_n_term, r_self = 1, r_other = 8
        ).to_openff()

        compute_mean_charges(
            terminal_name, terminal_residue, 0, charge_dir, mean_dir
        )

    # Charged N terminal residues
    for aa_name in aa_order:

        terminal_name = f'N{aa_name}'
        terminal_residue = Molecule.from_smiles(
            terminal_smiles[terminal_name], allow_undefined_stereo = True
        )

        compute_mean_charges(
            terminal_name, terminal_residue, 0, charge_dir, mean_dir
        )

    # C terminal residues, i.e. X in Ace-Val-Y-X. Do neutral first to mimic N
    #  terminal residues.
    for aa_name in aa_order:

        terminal_name = f'C0{aa_name}'
        terminal_monomer = Monomer.from_smiles(
            AMINO_ACIDS_ALL_STATES[aa_name], name = terminal_name
        )
        terminal_residue = terminal_monomer.substitute(
            neutral_c_term, r_self = 2, r_other = 8
        ).to_openff()

        compute_mean_charges(
            terminal_name, terminal_residue, 3, charge_dir, mean_dir
        )

    # Charged C terminal residues
    for aa_name in aa_order:

        terminal_name = f'C{aa_name}'
        terminal_residue = Molecule.from_smiles(
            terminal_smiles[terminal_name], allow_undefined_stereo = True
        )

        compute_mean_charges(
            terminal_name, terminal_residue, 3, charge_dir, mean_dir
        )


if __name__ == "__main__":
    main()

