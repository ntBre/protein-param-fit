from pathlib import Path

import click
from openff.toolkit import ForceField


@click.command()
@click.option(
    "-i",
    "--input_ff",
    default=Path("forcebalance", "result", "optimize", "force-field.offxml"),
    show_default=True,
    type=click.STRING,
    help="File path to initial force field from output of ForceBalance.",
)
@click.option(
    "-o",
    "--output_ff",
    default="final-force-field.offxml",
    show_default=True,
    type=click.STRING,
    help="File path to which the output force field will be written.",
)
def main(input_ff, output_ff):

    # Load initial force field from ForceBalance with cosmetic attributes
    force_field = ForceField(input_ff, allow_cosmetic_attributes=True)

    # Write force field without cosmetic attributes from ForceBalance
    unconstrained_ff = output_ff.replace(".offxml", "_unconstrained.offxml")
    force_field.to_file(unconstrained_ff, discard_cosmetic_attributes=True)

    # Add constraints to covalent bonds involving hydrogen
    constraint_handler = force_field["Constraints"]
    covalent_h_bond_constraint = {
        "id": "c1",
        "smirks": "[#1:1]-[*:2]",
    }
    constraint_handler.add_parameter(covalent_h_bond_constraint, before=0)

    # Write force field including covalent hydrogen bond constraint without
    # cosmetic attributes from ForceBalance
    force_field.to_file(output_ff, discard_cosmetic_attributes=True)


if __name__ == "__main__":
    main()
