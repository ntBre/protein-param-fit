import click
from openff.toolkit.typing.engines.smirnoff import ForceField
from pathlib import Path


@click.command()
@click.option(
    '-i',
    '--input_ff',
    default = Path('forcebalance', 'result', 'optimize', 'force-field.offxml'),
    show_default = True,
    type = click.STRING,
    help = 'File path to initial force field from output of ForceBalance.',
)
@click.option(
    '-o',
    '--output_ff',
    default = 'final-force-field.offxml',
    show_default = True,
    type = click.STRING,
    help = 'File path to which the output force field will be written.',
)
def main(input_ff, output_ff):

    # Load initial force field from ForceBalance with cosmetic attributes
    force_field = ForceField(input_ff, allow_cosmetic_attributes = True)

    # Write force field without cosmetic attributes
    force_field.to_file(output_ff, discard_cosmetic_attributes = True)


if __name__ == "__main__":
    main()

