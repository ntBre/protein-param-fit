import click
import json
from openff.toolkit.topology import Molecule
from pathlib import Path


@click.command()
@click.option(
    '-l',
    '--label_dir',
    type = click.STRING,
    default = 'sage-behara-smirks',
    show_default = True,
    help = 'Directory path to the labeld molecules.',
)
def main(label_dir):

    with open(Path(label_dir, 'optimization-labels.json'), 'r') as json_file:
        labels = json.load(json_file)

    with open(Path(label_dir, 'torsiondrive-labels.json'), 'r') as json_file:
        labels.extend(json.load(json_file))

    # SMIRKS for functional groups with delocalized charges
    delocalized_charge_smirks = {
        'Amidinium': {
            'vdW': '[#6X3](~!@[#7X3])~!@[#7X3:1]',
            'Bonds': '[#6X3:1](~!@[#7X3])~!@[#7X3:2]',
            'Angles': '[*:1]~[#6X3:2](~!@[#7X3])~!@[#7X3:3]',
            'ProperTorsions': '[*:1]~[*:2]~[#6X3:3](~!@[#7X3])~!@[#7X3:4]',
        },
        'Carboxylate': {
            'vdW': '[#6X3](~[#8X1])~[#8X1:1]',
            'Bonds': '[#6X3:1](~[#8X1])~[#8X1:2]',
            'Angles': '[*:1]~[#6X3:2](~[#8X1])~[#8X1:3]',
            'ProperTorsions': '[*:1]~[*:2]~[#6X3:3](~[#8X1])~[#8X1:4]',
        },
        'Nitro': {
            'vdW': '[#7X3](~[#8X1])~[#8X1:1]',
            'Bonds': '[#7X3:1](~[#8X1])~[#8X1:2]',
            'Angles': '[*:1]~[#7X3:2](~[#8X1])~[#8X1:3]',
            'ProperTorsions': '[*:1]~[*:2]~[#7X3:3](~[#8X1])~[#8X1:4]',
        },
        'Phosphate': {
            'vdW': '[#15X3](~[#8X1])~[#8X1:1]',
            'Bonds': '[#15X3:1](~[#8X1])~[#8X1:2]',
            'Angles': '[*:1]~[#15X3:2](~[#8X1])~[#8X1:3]',
            'ProperTorsions': '[*:1]~[*:2]~[#15X3:3](~[#8X1])~[#8X1:4]',
        },
        'Sulfate': {
            'vdW': '[#15X3](~[#8X1])~[#8X1:1]',
            'Bonds': '[#15X3:1](~[#8X1])~[#8X1:2]',
            'Angles': '[*:1]~[#15X3:2](~[#8X1])~[#8X1:3]',
            'ProperTorsions': '[*:1]~[*:2]~[#15X3:3](~[#8X1])~[#8X1:4]',
        },
    }

    # Find chemically equivalent atoms that are assigned different parameters
    processed_smiles = []
    for molecule_labels in labels:

        if molecule_labels['smiles'] in processed_smiles:
            continue

        processed_smiles.append(molecule_labels['smiles'])
        offmol = Molecule.from_mapped_smiles(
            molecule_labels['smiles'], allow_undefined_stereo = True
        )
        different_assignments = []

        for func_group, parameter_smirks in delocalized_charge_smirks.items():

            for parameter_type, smirks in parameter_smirks.items():

                parameter_ids = molecule_labels[parameter_type]

                matches = offmol.chemical_environment_matches(smirks)

                # Only interested in two or more matches to compare
                if len(matches) < 2:
                    continue

                # Get dict of parameter ids
                match_ids = dict()

                for match in matches:

                    str_match = str(match)

                    if str_match in parameter_ids:
                        match_ids[match] = parameter_ids[str_match]
                    else:
                        match_ids[match] = parameter_ids[str(match[::-1])]

                # Print matches and ids if parameter ids are not all equal
                if len(set(match_ids.values())) > 1:

                    for match, parameter_id in match_ids.items():

                        indices = ' '.join(f'{int(i) + 1:2d}' for i in match)
                        different_assignments.append(
                            f'{func_group:11s} {parameter_type:14s} '
                            f'{indices:>11s} {parameter_id}'
                        )

        if len(different_assignments) > 0:

            print('Molecule', offmol.to_smiles(mapped = True))

            for line in different_assignments:
                print(line)


if __name__ == "__main__":
    main()

