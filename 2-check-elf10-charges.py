import click
from openff.qcsubmit.results import (
    OptimizationResultCollection, TorsionDriveResultCollection
)
from openff.toolkit.topology import Molecule
from pathlib import Path
from tqdm import tqdm

def remove_elf_failures(training_set, output_to_save):

    entries = list(
        training_set.entries['https://api.qcarchive.molssi.org:443/']
    )

    index_to_delete = []
    cmiles = []
    failed_cmiles = []

    for i, entry in tqdm(enumerate(entries)):

        if entry.cmiles in failed_cmiles:

            print(entry.record_id)
            index_to_delete.append(i)

        elif entry.cmiles in cmiles:
            pass

        else:

            cmiles.append(entry.cmiles)
            mol = Molecule.from_mapped_smiles(
                entry.cmiles, allow_undefined_stereo = True
            )

            try:
                mol.assign_partial_charges(partial_charge_method='am1bccelf10')

            except:

                print(entry.record_id)
                failed_cmiles.append(entry.cmiles)
                index_to_delete.append(i)

    for i in sorted(index_to_delete, reverse=True):
        del training_set.entries['https://api.qcarchive.molssi.org:443/'][i]

    with open(output_to_save, 'w') as json_file:
        json_file.write(training_set.json())


@click.command()
@click.option(
    '-d',
    '--dataset_dir',
    default = 'training-datasets',
    show_default = True,
    type = click.STRING,
    help = 'Directory path to which training datasets will be written.',
)
def main(dataset_dir):

    torsion_datasets = [
        'small-molecule', 'protein-weight1', 'protein-weight2',
        'protein-weight3', 'protein-weight4', 'protein-weight6',
    ]

    for dataset in torsion_datasets:

        torsion_training_set = TorsionDriveResultCollection.parse_file(
            Path(dataset_dir, f'{dataset}-torsiondrive-dataset.json')
        )

        remove_elf_failures(
            torsion_training_set,
            Path(
                dataset_dir,
                f'{dataset}-torsiondrive-elf10-filter-dataset.json'
            )
        )

    optimization_training_set = OptimizationResultCollection.parse_file(
        Path(dataset_dir, 'optimization-dataset.json')
    )

    remove_elf_failures(
        optimization_training_set,
        Path(dataset_dir, 'optimization-elf10-filter-dataset.json')
    )


if __name__ == "__main__":
    main()

