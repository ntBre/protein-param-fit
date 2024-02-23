from pathlib import Path

import click
from openff.toolkit import ForceField


@click.command()
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
    input_force_field: str,
    library_charge_force_field: str,
    output_force_field: str
):
    # Combine initial small molecule force field and protein library charges
    force_field = ForceField(input_force_field, library_charge_force_field)
    force_field.to_file(output_force_field)


if __name__ == "__main__":
    main()
