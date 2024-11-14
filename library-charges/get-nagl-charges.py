from amber_atom_name_map import AMBER_ATOM_NAME_MAP
import click
import numpy
from openff.toolkit import ForceField, Molecule
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from openff.units import unit
import pandas


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


@click.command()
@click.option(
    "-a",
    "--amber-charges-path",
    default="amber-ff99-charges.dat",
    show_default=True,
    type=click.STRING,
    help="The path to the file containing the Amber charges.",
)
@click.option(
    "-l",
    "--library-charges-path",
    default="protein-library-charges.offxml",
    show_default=True,
    type=click.STRING,
    help="The path to the file containing the library charges.",
)
@click.option(
    "-n",
    "--nagl-model",
    default="openff-gnn-am1bcc-0.1.0-rc.1.pt",
    show_default=True,
    type=click.STRING,
    help="The file name or file path of the NAGL model to use for inference.",
)
def main(amber_charges_path, library_charges_path, nagl_model):

    # Read reference Amber charges into a dictionary
    amber_charges = dict()
    with open(amber_charges_path, 'r') as amber_charges_file:

        for line in amber_charges_file:

            fields = line.split()
            if fields[0] == "#":
                continue

            res_name = fields[0]
            if res_name == 'CYM':
                continue

            if res_name not in amber_charges:
                amber_charges[res_name] = list()

            atom_name = fields[1].strip('\"')

            if (
                res_name in AMBER_ATOM_NAME_MAP
                and atom_name in AMBER_ATOM_NAME_MAP[res_name]
            ):

                atom_name = AMBER_ATOM_NAME_MAP[res_name][atom_name] 

            amber_charges[res_name].append((atom_name, float(fields[3])))

    # Load library charge force field
    ff = ForceField(library_charges_path)
    library_charge_handler = ff['LibraryCharges']

    charges = list()
    for residue in RESIDUE_SMILES.keys():
        if residue in {"ACE", "NME"}:
            continue

        library_charge_parameter = library_charge_handler.get_parameter(
            {"id": f"Protein-{residue}"}
        )[0]
        residue_smirks = (
            "[#6X4]-[#6X3](=[#8])-[#7X3]-[#6X4]-[#6X3](=[#8])-[#7X3]-"
            + library_charge_parameter.smirks
            + "-[#6X4]-[#6X3](=[#8])-[#7X3]-[#6X4]-[#6X3](=[#8])-[#7X3]-[#6X4]"
        )

        ala_smiles = RESIDUE_SMILES["ALA"]
        ace_smiles = RESIDUE_SMILES["ACE"]
        nme_smiles = RESIDUE_SMILES["NME"]
        res_smiles = RESIDUE_SMILES[residue]
        mol_smiles = (
            f"{ace_smiles}{ala_smiles}{ala_smiles}{res_smiles}{ala_smiles}"
            + f"{ala_smiles}{nme_smiles}"
        )

        elf10_mol = Molecule.from_smiles(mol_smiles)
        elf10_atom_map = elf10_mol.chemical_environment_matches(
            residue_smirks,
            unique=True,
        )[0]
        elf10_mol.assign_partial_charges("am1bccelf10")

        nagl_mol = Molecule.from_smiles(mol_smiles)
        nagl_atom_map = nagl_mol.chemical_environment_matches(
            residue_smirks,
            unique=True,
        )[0]
        nagl_mol.assign_partial_charges(
            nagl_model,
            toolkit_registry=NAGLToolkitWrapper()
        )

        for i, (atom_name, amber_charge) in enumerate(amber_charges[residue]):
            library_charge = library_charge_parameter.charge[i]
            elf10_charge = elf10_mol.partial_charges[elf10_atom_map[i]]
            nagl_charge = nagl_mol.partial_charges[nagl_atom_map[i]]
            charges.append(
                {
                    "Residue": residue,
                    "Atom": atom_name,
                    "Amber ff99": amber_charge,
                    "Library ELF10": library_charge / unit.elementary_charge,
                    "OpenEye ELF10": elf10_charge / unit.elementary_charge,
                    "NAGL 0.1.0rc1": nagl_charge / unit.elementary_charge,
                }
            )

    pandas.DataFrame(charges).to_csv("nagl-charges.dat")

if __name__ == "__main__":
    main()

