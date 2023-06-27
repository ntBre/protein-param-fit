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
    '-n',
    '--nh2_library_charge_ff',
    default = Path(
        '..', 'library-charges', 'protein-nh2-library-charges.offxml'
    ),
    show_default = True,
    type = click.STRING,
    help = 'File path to which the output force field will be written.',
)
@click.option(
    '-o',
    '--output_ff',
    default = 'final-force-field-nh2.offxml',
    show_default = True,
    type = click.STRING,
    help = 'File path to which the output force field will be written.',
)
def main(input_ff, nh2_library_charge_ff, output_ff):

    # Load initial force field from ForceBalance with cosmetic attributes
    force_field = ForceField(input_ff, allow_cosmetic_attributes = True)

    # Load force field with NH2 library charges
    nh2_force_field = ForceField(nh2_library_charge_ff)

    # Copy NH2 library charges to new force field after ACE library charges
    nh2_handler = nh2_force_field.get_parameter_handler('LibraryCharges')
    reference_parameter = nh2_handler.get_parameter({'id': 'Protein-NH2'})[0]

    new_parameter = {
        'id': reference_parameter.id,
        'smirks': reference_parameter.smirks,
    }

    for i, charge in enumerate(reference_parameter.charge):
        new_parameter[f'charge{i + 1}'] = charge

    library_handler = force_field.get_parameter_handler('LibraryCharges')
    ace_parameter = library_handler.get_parameter({'id': 'Protein-ACE'})[0]
    library_handler.add_parameter(new_parameter, after = ace_parameter.smirks)

    # Replace -[#7X3]-[#6X4] with SMARTS that matches amide cap for all protein
    # residues except ACE
    for parameter in library_handler.parameters:

        if (
            parameter.id.startswith('Protein-')
            and parameter.smirks.endswith('-[#7X3]-[#6X4]')
            and parameter.id != 'Protein-ACE'
        ):

            # Truncate the last 14 characters, i.e. '-[#7X3]-[#6X4]'
            truncate_position = len(parameter.smirks) - 14

            parameter.smirks = (
                parameter.smirks[:truncate_position] + '-[#7X3;$(*-[#6X4]),H2]'
            )

    # Write force field without cosmetic attributes from ForceBalance
    unconstrained_ff = output_ff.replace('.offxml', '_unconstrained.offxml')
    force_field.to_file(unconstrained_ff, discard_cosmetic_attributes = True)

    # Add constraints to covalent bonds involving hydrogen
    constraint_handler = force_field.get_parameter_handler('Constraints')
    covalent_h_bond_constraint = {
        'id': 'c1',
        'smirks': '[#1:1]-[*:2]',
    }
    constraint_handler.add_parameter(covalent_h_bond_constraint, before = 0)

    # Write force field including covalent hydrogen bond constraint without
    # cosmetic attributes from ForceBalance
    force_field.to_file(output_ff, discard_cosmetic_attributes = True)


if __name__ == "__main__":
    main()

