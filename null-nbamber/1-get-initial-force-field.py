from pathlib import Path
import re

import click
from openff.toolkit import ForceField
from openff.units import unit


@click.command()
@click.option(
    "-a",
    "--amber-nonbonded-path",
    "amber_nonbonded_path",
    default="NBAmber.dat",
    show_default=True,
    type=click.STRING,
    help="Path to force field containing protein library charges.",
)
@click.option(
    "-i",
    "--input-force-field",
    "input_force_field",
    default="openff_unconstrained-2.1.0.offxml",
    show_default=True,
    type=click.STRING,
    help="Name of or path to initial small molecule force field.",
)
@click.option(
    "-l",
    "--library-charge-force-field",
    "library_charge_force_field",
    default=Path("library-charges", "protein-library-charges.offxml"),
    show_default=True,
    type=click.STRING,
    help="Path to force field containing protein library charges.",
)
@click.option(
    "-o",
    "--output-force-field",
    "output_force_field",
    default="initial-force-field.offxml",
    show_default=True,
    type=click.STRING,
    help="File path to which the output force field will be written.",
)
def main(
    amber_nonbonded_path: str,
    input_force_field: str,
    library_charge_force_field: str,
    output_force_field: str
):
    # Read initial small molecule force field
    force_field = ForceField(input_force_field)

    # Read force field with protein library charge SMIRKS
    library_charge_smirks_ff = ForceField(library_charge_force_field)

    # Read Amber nonbonded parameters
    amber_nonbonded_params = dict()
    with open(Path(amber_nonbonded_path), "r") as nonbonded_file:
        for line in nonbonded_file:
            fields = line.split()
            residue_name = fields[0]
            if residue_name not in amber_nonbonded_params:
                amber_nonbonded_params[residue_name] = {
                    "Charge": list(),
                    "LJ": list(),
                }
            amber_nonbonded_params[residue_name]["Charge"].append(float(fields[2]))
            lj_tuple = (fields[1], float(fields[3]), float(fields[4]))
            amber_nonbonded_params[residue_name]["LJ"].append(lj_tuple)

    # Add Amber nonbonded parameters to protein library charge SMIRKS
    library_charge_handler = force_field["LibraryCharges"]
    vdW_handler = force_field["vdW"]
    for parameter in library_charge_smirks_ff["LibraryCharges"]:
        # Skip neutral terminal residues N0 and C0
        if "0" in parameter.id:
            continue

        # Skip residues not in Amber's library charges
        residue_name = parameter.id.split("-")[1]
        if residue_name in {"NASH", "NGLH", "NLYN", "CASH", "CGLH", "CLYN"}:
            continue

        # Add Amber library charges
        library_charge_dict = {"id": f"Protein-{residue_name}", "smirks": parameter.smirks}
        for index, charge in enumerate(amber_nonbonded_params[residue_name]["Charge"]):
            library_charge_dict[f"charge{index + 1}"] = charge * unit.elementary_charge
        library_charge_handler.add_parameter(library_charge_dict)

        # Add Amber Lennard-Jones parameters
        for index, (atom_name, R_min, epsilon) in enumerate(amber_nonbonded_params[residue_name]["LJ"]):
            new_smirks = re.sub(
                r"\[(.+?):([0-9]+?)\]",
                lambda m: f"[{m.group(1)}:1]" if m.group(2) == str(index + 1) else f"[{m.group(1)}]",
                parameter.smirks,
            )
            vdw_dict = {
                "id": f"Protein-LJ-{residue_name}-{atom_name}",
                "smirks": new_smirks,
                "rmin_half": R_min * unit.angstrom,
                "epsilon": epsilon * unit.kilocalorie_per_mole,
            }
            vdW_handler.add_parameter(vdw_dict)

    force_field.to_file(output_force_field)


if __name__ == "__main__":
    main()
