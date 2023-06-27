from aa_library_charges import INDEX_TO_AA
from amber_atom_name_map import AMBER_ATOM_NAME_MAP, POLYMETRIZER_ATOM_NAME_MAP
import click
from matplotlib import pyplot
import numpy
from openff.toolkit.topology import Molecule
from openmm import unit
import os
from polymetrizer.tests.smiles import AMINO_ACIDS_ALL_STATES
import seaborn


pyplot.style.use('dark_background')
seaborn.set_palette("colorblind")


# Compute mean charges over X and Z for Ace-Val-X-Y-Z-Val-Nme
def get_mean_charges(
    aa_index, charge_dir = 'aa_library_charges', mean_dir = 'mean_charges'
):

    aa_name = INDEX_TO_AA[aa_index]
    central_residue = Molecule.from_smiles(
        AMINO_ACIDS_ALL_STATES[aa_name], allow_undefined_stereo = True
    )

    # Get atom indices of central residue in Ace-Val-X-Y-Z-Val-Nme.
    # Since the oligomer was built from the center out, the central residue will
    # have the lowest indices.
    atom_idx = 0
    atom_map = dict()
    name_map = dict()

    for i, atom in enumerate(central_residue.atoms):
        if atom.atomic_number:

            atom_idx += 1
            atom_map[i] = atom_idx
            name_map[atom_idx - 1] = f'{atom.element.symbol}{atom_idx:d}'

    num_atoms = atom_idx
    central_residue.properties['atom_map'] = atom_map

    sum_charges = numpy.zeros(num_atoms)
    sum_weights = 0

    with open(
        os.path.join(mean_dir, f'{aa_name.lower()}_charges.dat'), 'w'
    ) as out_file:

        for charge_dir_file in os.listdir(charge_dir):

            if (
                not charge_dir_file.startswith(aa_name + '_')
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

            if residues[2] != aa_name:

                raise ValueError(
                    f'Central residue {residues[2]} is not {aa_name}'
                )

            mol_charges = numpy.loadtxt(
                os.path.join(charge_dir, f'{mol_name}.dat')
            )

            # Weight for Ace-Val-Arg-Arg-Cys-Val-Nme is 1.0 because ELF10 failed
            # for Ace-Val-Arg-Arg-Cyx-Val-Nme
            if mol_name == 'Arg_Ace-Val-Arg-Arg-Cys-Val-Nme':
                weight = 1.0
            else:
                weight = numpy.prod([AA_WEIGHTS[res] for res in residues])

            sum_charges += weight * mol_charges[:num_atoms]
            sum_weights += weight

            central_charges = [f'{q:11.8f}' for q in mol_charges[:num_atoms]]
            central_charge_str = ' '.join(central_charges)

            out_file.write(residue_str + ' ' + central_charge_str + '\n')

    # Get mean charge over all oligomers
    mean_charges = sum_charges / sum_weights

    # Distribute non-integer charge
    expected_charge = central_residue.total_charge.value_in_unit(
        unit.elementary_charge
    )
    total_charge = numpy.sum(mean_charges)
    if total_charge != expected_charge:
        mean_charges += (expected_charge - total_charge) / num_atoms

    with open(
        os.path.join('mean_charges', 'mean_charges.dat'), 'a'
    ) as out_file:

        for i in range(num_atoms):

            out_file.write(
                f'{aa_name} {name_map[i]:3s} {mean_charges[i]:11.8f}\n'
            )


@click.command()
@click.option(
    "-a",
    "--amber",
    "amber_charges_path",
    default="amber_ff99_charges",
    show_default=True,
    type=click.STRING,
    help="The path to the file containing the Amber charges.",
)
@click.option(
    "-l",
    "--library",
    "library_charges_path",
    default=os.path.join("mean_charges", "mean_charges.dat"),
    show_default=True,
    type=click.STRING,
    help="The path to the file containing the library charges.",
)
def main(library_charges_path, amber_charges_path):

    # Read reference Amber charges into a dictionary
    amber_label = 'Amber ff99'
    charges = {amber_label: dict()}
    with open(amber_charges_path, 'r') as amber_charges_file:

        for line in amber_charges_file:

            fields = line.split()
            if fields[0] == "#":
                continue

            res_name = fields[0]
            if res_name == 'CYM':
                continue

            if res_name not in charges[amber_label]:
                charges[amber_label][res_name] = dict()

            atom_name = fields[1].strip('\"')

            if (
                res_name in AMBER_ATOM_NAME_MAP
                and atom_name in AMBER_ATOM_NAME_MAP[res_name]
            ):

                atom_name = AMBER_ATOM_NAME_MAP[res_name][atom_name] 

            charges[amber_label][res_name][atom_name] = float(fields[3])

    # Read library charges into a dictionary
    library_label = 'ELF10'
    charges[library_label] = dict()
    with open(library_charges_path, 'r') as library_charges_file:

        for line in library_charges_file:

            fields = line.split()

            res_name = fields[0].upper()
            if res_name not in charges[library_label]:
                charges[library_label][res_name] = dict()

            atom_name = POLYMETRIZER_ATOM_NAME_MAP[res_name][fields[1]]

            if atom_name in charges[amber_label][res_name]:
                charges[library_label][res_name][atom_name] = float(fields[2])
            else:
                print(f'Atom {atom_name} in residue {res_name} not in Amber')

    # Plot library charges and Amber charges
    charge_types = [*charges]
    bar_width = 0.8 / len(charge_types)
    bar_location_offset = (len(charge_types) - 1) / 2

    for res_name in charges[amber_label]:

        if res_name not in charges[library_label]:
            continue

        atoms = [*charges[amber_label][res_name]]
        tick_locations = numpy.arange(len(atoms))
        figure = pyplot.figure(figsize = tuple(4.25 * x for x in (1, 0.75)))

        for i, charge_type in enumerate(charge_types):

            bar_locations = (
                tick_locations + (i - bar_location_offset) * bar_width
            )
            bar_heights = [
                charges[charge_type][res_name][atom] for atom in atoms
            ]

            pyplot.bar(
                bar_locations, bar_heights, width=bar_width, label=charge_type
            )

        pyplot.xticks(tick_locations, labels=atoms, fontsize=8)
        pyplot.ylim(-1, 1)
        pyplot.yticks(numpy.arange(-1, 1.1, 0.2))
        pyplot.ylabel('Charge (e)')
        pyplot.legend(
            bbox_to_anchor = (0, 1, 1, 0), loc = 'lower left', ncol = 2
        )

        pyplot.savefig(os.path.join('plots', f'{res_name.lower()}_charges.pdf'))
        pyplot.close(figure)


if __name__ == "__main__":
    main()

