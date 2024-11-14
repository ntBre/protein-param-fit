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


RESIDUE_SMILES = {
    'ALA': 'N[C@@H](C)C(=O)',
    'ARG': 'N[C@@H](CCC[NH+]=C(N)N)C(=O)',
    'ASH': 'N[C@@H](CC(=O)O)C(=O)',
    'ASN': 'N[C@@H](CC(=O)N)C(=O)',
    'ASP': 'N[C@@H](CC(=O)[O-])C(=O)',
    'CYS': 'N[C@@H](CS)C(=O)',
    'CYX': 'N[C@@H](CSSC[C@@H](C(=O)NC)NC(=O)C)C(=O)',
    'GLH': 'N[C@@H](CCC(=O)O)C(=O)',
    'GLN': 'N[C@@H](CCC(=O)N)C(=O)',
    'GLU': 'N[C@@H](CCC(=O)[O-])C(=O)',
    'GLY': 'NCC(=O)',
    'HID': 'N[C@@H](CC1=CN=CN1)C(=O)',
    'HIE': 'N[C@@H](CC1=CNC=N1)C(=O)',
    'HIP': 'N[C@@H](CC1=CNC=[NH+]1)C(=O)',
    'ILE': 'N[C@@H]([C@@H](C)CC)C(=O)',
    'LEU': 'N[C@@H](CC(C)C)C(=O)',
    'LYN': 'N[C@@H](CCCCN)C(=O)',
    'LYS': 'N[C@@H](CCCC[NH3+])C(=O)',
    'MET': 'N[C@@H](CCSC)C(=O)',
    'PHE': 'N[C@@H](Cc1ccccc1)C(=O)',
    'PRO': 'N1CCC[C@H]1C(=O)',
    'SER': 'N[C@@H](CO)C(=O)',
    'THR': 'N[C@@H]([C@H](O)C)C(=O)',
    'TRP': 'N[C@@H](CC1=CNc2c1cccc2)C(=O)',
    'TYR': 'N[C@@H](Cc1ccc(cc1)O)C(=O)',
    'VAL': 'N[C@@H](C(C)C)C(=O)',
    'ACE': 'CC(=O)',
    'NME': 'NC',
}


pyplot.style.use('dark_background')
seaborn.set_palette("colorblind")


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

