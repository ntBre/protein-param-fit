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
    "-n",
    "--nagl-model",
    "nagl_model",
    default="openff-gnn-am1bcc-0.1.0-rc.2.pt",
    show_default=True,
    type=click.STRING,
    help="Name of or path to the NAGL model for assigning partial charges.",
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
    nagl_model: str,
    output_force_field: str
):
    # Replace ToolkitAM1BCC handler with ChargeIncremementModel using NAGL
    force_field = ForceField(input_force_field)
    force_field.deregister_parameter_handler("ToolkitAM1BCC")
    force_field.get_parameter_handler(
        "ChargeIncrementModel",
        {
            "version": 0.3,
            "partial_charge_method": nagl_model,
        },
    )
    force_field.to_file(output_force_field)


if __name__ == "__main__":
    main()
