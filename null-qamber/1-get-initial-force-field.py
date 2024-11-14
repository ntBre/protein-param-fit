from pathlib import Path

import click
from openff.toolkit import ForceField
from openff.units import unit


@click.command()
@click.option(
    "-a",
    "--amber-charges-path",
    "amber_charges_path",
    default="QAmber.dat",
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
    amber_charges_path: str,
    input_force_field: str,
    library_charge_force_field: str,
    output_force_field: str
):
    # Read initial small molecule force field
    force_field = ForceField(input_force_field)

    # Read force field with protein library charge SMIRKS
    library_charge_smirks_ff = ForceField(library_charge_force_field)

    # Read Amber library charges
    amber_charges = dict()
    with open(Path(amber_charges_path), "r") as charge_file:
        for line in charge_file:
            fields = line.split()
            residue_name = fields[0]
            if residue_name not in amber_charges:
                amber_charges[residue_name] = list()
            amber_charges[residue_name].append(float(fields[2]))

    # Add Amber library charges to protein library charge SMIRKS
    library_charge_handler = force_field["LibraryCharges"]
    for parameter in library_charge_smirks_ff["LibraryCharges"]:
        # Skip neutral terminal residues N0 and C0
        if "0" in parameter.id:
            continue
        library_charge_dict = {"id": parameter.id, "smirks": parameter.smirks}
        residue_name = parameter.id.split("-")[1]
        #Skip residues not in Amber's library charges
        if residue_name in {"NASH", "NGLH", "NLYN", "CASH", "CGLH", "CLYN"}:
            continue
        for index, charge in enumerate(amber_charges[residue_name]):
            library_charge_dict[f"charge{index+1}"] = charge * unit.elementary_charge
        library_charge_handler.add_parameter(library_charge_dict)

    force_field.to_file(output_force_field)


if __name__ == "__main__":
    main()
