import click
from openff.toolkit.typing.engines.smirnoff import ForceField
from pathlib import Path


@click.command()
@click.option(
    '-i',
    '--input_ff',
    default = Path('..', 'sage-behara.offxml'),
    show_default = True,
    type = click.STRING,
    help = 'File path to initial small molecule force field.',
)
@click.option(
    '-o',
    '--output_ff',
    default = 'initial-force-field.offxml',
    show_default = True,
    type = click.STRING,
    help = 'File path to which the output force field will be written.',
)
def main(input_ff, output_ff):

    # Combine initial small molecule force field and protein library charges
    force_field = ForceField(input_ff)

    # Change LibraryCharge parameters for ions to use "id" instead of "name"
    library_charge_handler = force_field.get_parameter_handler('LibraryCharges')
    for parameter in library_charge_handler.parameters:

        if parameter.name != None:

            parameter.id = parameter.name
            parameter.name = None

    # Define new torsion types for delocalized charges that are assigned
    # incorrect torsions
    # (see https://github.com/openforcefield/openff-forcefields/issues/56)
    delocalized_charge_torsion_types = [
        # Amidinium
        {
            'label': 't18b',
            'reference_id': 't18',
            'smirks': '[*:1]-[#6X4:2]-[#6X3:3](~!@[#7X3])~!@[#7X3:4]',
        },
        # Carboxylate
        {
            'label': 't18a',
            'reference_id': 't18',
            'smirks': '[*:1]-[#6X4:2]-[#6X3:3](~[#8X1])~[#8X1:4]',
        },
        {
            'label': 't19a',
            'reference_id': 't19',
            'smirks': '[#1:1]-[#6X4:2]-[#6X3:3](~[#8X1])~[#8X1:4]',
        },
        {
            'label': 't31a',
            'reference_id': 't31',
            'smirks': '[#6X3:1]-[#6X4;r3:2]-[#6X3:3](~[#8X1])~[#8X1:4]',
        },
        {
            'label': 't42a',
            'reference_id': 't42',
            'smirks': '[#6X4;r3:1]-;@[#6X4;r3:2]-[#6X3:3](~[#8X1])~[#8X1:4]',
        },
        {
            'label': 't48a',
            'reference_id': 't48',
            'smirks': '[#6X3:1]=[#6X3:2]-[#6X3:3](~[#8X1])~[#8X1:4]',
        },
        # Nitro from sp2 carbon (nitro from sp3 carbon is covered by t65)
        {
            'label': 't82a',
            'reference_id': 't82',
            'smirks': '[*:1]-[#6X3:2]-[#7X3:3](~[#8X1])~[#8X1:4]',
        },
        {
            'label': 't83a',
            'reference_id': 't83',
            'smirks': '[*:1]=,:[#6X3:2]-[#7X3:3](~[#8X1])~[#8X1:4]',
        },
        # Guanidinium
        {
            'label': 't87a',
            'reference_id': 't87',
            'smirks': '[*:1]-[#7X3:2]~!@[#6X3:3](~!@[#7X3])~!@[#7X3:4]',
        },
    ]

    # Guanidinium bonds
    delocalized_charge_bond_types = [
        {
            'label': 'b13a',
            'reference_id': 'b13',
            'smirks': '[#6X3:1](~!@[#7X3])(~!@[#7X3])~!@[#7X3:2]',
        },
    ]

    # Add delocalized charge torsions to force field
    torsion_handler = force_field.get_parameter_handler('ProperTorsions')

    for torsion_type in delocalized_charge_torsion_types:

        # Copy parameters from an existing torsion type
        reference_parameter = torsion_handler.get_parameter(
            {'id': torsion_type['reference_id']}
        )[0]

        new_parameter = {
            'id': torsion_type['label'],
            'smirks': torsion_type['smirks'],
        }

        for i, n in enumerate(reference_parameter.periodicity):

            new_parameter[f'periodicity{i+1}'] = n
            new_parameter[f'k{i+1}'] = reference_parameter.k[i]
            new_parameter[f'phase{i+1}'] = reference_parameter.phase[i]
            new_parameter[f'idivf{i+1}'] = reference_parameter.idivf[i]

        torsion_handler.add_parameter(
            new_parameter, after = reference_parameter.smirks
        )

    # Add delocalized charge bonds to force field
    bond_handler = force_field.get_parameter_handler('Bonds')

    for bond_type in delocalized_charge_bond_types:

        # Copy parameters from an existing torsion type
        reference_parameter = bond_handler.get_parameter(
            {'id': bond_type['reference_id']}
        )[0]

        new_parameter = {
            'id': bond_type['label'],
            'smirks': bond_type['smirks'],
            'length': reference_parameter.length,
            'k': reference_parameter.k
        }

        bond_handler.add_parameter(
            new_parameter, after = reference_parameter.smirks
        )

    force_field.to_file(output_ff)


if __name__ == "__main__":
    main()

