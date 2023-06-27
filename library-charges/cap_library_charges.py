from aa_library_charges import (
    INDEX_TO_AA, ace, nme, compute_oligomer_charges, extend_n_terminus,
    extend_c_terminus
)
import click
import os
from polymetrizer import Monomer
from polymetrizer.cap import Cap
from polymetrizer.tests.smiles import AMINO_ACIDS_ALL_STATES


@click.command()
@click.argument(
    'a1_index',
    type=click.INT,
)
@click.argument(
    'a2_index',
    type=click.INT,
)
@click.argument(
    'a3_index',
    type=click.INT,
)
@click.argument(
    'cap_name',
    type=click.STRING,
)
def main(a1_index, a2_index, a3_index, cap_name):

    cap_name = cap_name.capitalize()

    if cap_name == 'Nh2':
        cap = Cap.from_smiles('[*:9]N', name = 'Nh2')
    elif cap_name == 'Nme':
        cap = nme
    else:
        raise ValueError(f'Argument `cap_name` must be one of "Nh2" or "Nme".')

    # Build Ace-Val-Cyx-Val-Nme disulfide linkage
    build_disulfide_linkage = False
    for aa_index in [a1_index, a2_index, a3_index]:
        if aa_index == 0 or aa_index == 7:
            build_disulfide_linkage = True
            break

    if build_disulfide_linkage:

        cyx = Monomer.from_smiles(AMINO_ACIDS_ALL_STATES['Cyx'])
        val = Monomer.from_smiles(AMINO_ACIDS_ALL_STATES['Val'])

        disulfide_linkage = cyx.substitute(val, r_self=1, r_other=2)
        disulfide_linkage.substitute(ace, r_self=1, r_other=6, inplace=True)
        disulfide_linkage.substitute(val, r_self=2, r_other=1, inplace=True)
        disulfide_linkage.substitute(nme, r_self=2, r_other=7, inplace=True)

    else:
        disulfide_linkage = None

    # Get list of indices for each position. 0 indicates loop over all indices.
    flanking_indices = []

    for flanking_index in [a2_index, a1_index]:
        if flanking_index == 0:
            flanking_indices.append(range(1, len(INDEX_TO_AA.keys()) + 1))
        else:
            flanking_indices.append([flanking_index])

    # Build molecules using polymetrizer, starting with the terminal residue.
    terminal_name = INDEX_TO_AA[a3_index]
    terminal_monomer = Monomer.from_smiles(
        AMINO_ACIDS_ALL_STATES[terminal_name], name=terminal_name
    )
    terminal_residue = terminal_monomer.substitute(cap, r_self=2, r_other=9)

    if terminal_name == 'Cyx':

        monomer = terminal_residue.substitute(
            disulfide_linkage, r_self=3, r_other=3
        )

    else:

        monomer = terminal_residue

    # Build toward Ace-capped terminus
    for flanking_index_1 in flanking_indices[0]:

        flanking_name_1 = INDEX_TO_AA[flanking_index_1]

        dimer = extend_n_terminus(
            monomer, flanking_name_1, disulfide_linkage
        )

        for flanking_index_2 in flanking_indices[1]:

            flanking_name_2 = INDEX_TO_AA[flanking_index_2]

            # Check for existing file to allow batch re-submits
            mol_name = (
                f'{cap_name}_Ace-{flanking_name_2}-{flanking_name_1}-'
                f'{terminal_name}-{cap_name}'
            )

            file_name = os.path.join(
                'aa_library_charges', f'{mol_name}.dat'
            )

            if os.path.isfile(file_name):
                continue

            trimer = extend_n_terminus(
                dimer, flanking_name_2, disulfide_linkage
            )
            trimer.substitute(ace, r_self=1, r_other=2, inplace=True)

            # Compute charges
            compute_oligomer_charges(trimer, mol_name, file_name)


if __name__ == "__main__":
    main()

