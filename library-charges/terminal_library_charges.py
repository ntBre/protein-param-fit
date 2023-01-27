from aa_library_charges import (
    INDEX_TO_AA, ace, nme, compute_oligomer_charges, extend_n_terminus,
    extend_c_terminus
)
import click
import os
from polymetrizer import Monomer
from polymetrizer.cap import Cap
from polymetrizer.tests.smiles import AMINO_ACIDS_ALL_STATES


terminal_smiles = {
    'CAla': '[H][N]([C@]([H])([C](=[O])[O-])[C]([H])([H])[H])[*:1]',
    'CArg': '[H][N]([H])[C](=[N+]([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1])[N]([H])[H]',
    'CAsn': '[H][N]([H])[C](=[O])[C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1]',
    'CAsp': '[H][N]([C@]([H])([C](=[O])[O-])[C]([H])([H])[C](=[O])[O-])[*:1]',
    'CAsh': '[H][N]([C@]([H])([C](=[O])[O-])[C]([H])([H])[C](=[O])[O][H])[*:1]',
    'CCys': '[H][S][C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1]',
    'CCyx': '[S]([*:3])[C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1]',
    'CGlh': '[H][N]([C@]([H])([C](=[O])[O-])[C]([H])([H])[C]([H])([H])[C](=[O])[O][H])[*:1]',
    'CGln': '[H][N]([H])[C](=[O])[C]([H])([H])[C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1]',
    'CGlu': '[H][N]([C@]([H])([C](=[O])[O-])[C]([H])([H])[C]([H])([H])[C](=[O])[O-])[*:1]',
    'CGly': '[H][N]([C]([H])([H])[C](=[O])[O-])[*:1]',
    'CHid': '[H][c]1[n]([H])[c]([C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1])[c]([H])[n]1',
    'CHie': '[H][c]1[n][c]([C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1])[c]([H])[n]1[H]',
    'CHip': '[H][c]1[n+]([H])[c]([C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1])[c]([H])[n]1[H]',
    'CIle': '[H][N]([C@]([H])([C](=[O])[O-])[C@@]([H])([C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[H])[*:1]',
    'CLeu': '[H][N]([C@]([H])([C](=[O])[O-])[C]([H])([H])[C]([H])([C]([H])([H])[H])[C]([H])([H])[H])[*:1]',
    'CLyn': '[H][N]([C@]([H])([C](=[O])[O-])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[N]([H])[H])[*:1]',
    'CLys': '[H][N]([C@]([H])([C](=[O])[O-])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[N+]([H])([H])[H])[*:1]',
    'CMet': '[H][N]([C@]([H])([C](=[O])[O-])[C]([H])([H])[C]([H])([H])[S][C]([H])([H])[H])[*:1]',
    'CPhe': '[H][c]1[c]([H])[c]([H])[c]([C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1])[c]([H])[c]1[H]',
    'CPro': '[H][C]1([H])[N]([*:1])[C@]([H])([C](=[O])[O-])[C]([H])([H])[C]1([H])[H]',
    'CSer': '[H][O][C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1]',
    'CThr': '[H][O][C@]([H])([C]([H])([H])[H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1]',
    'CTrp': '[H][c]1[c]([H])[c]([H])[c]2[c]([c]1[H])[c]([C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1])[c]([H])[n]2[H]',
    'CTyr': '[H][O][c]1[c]([H])[c]([H])[c]([C]([H])([H])[C@@]([H])([C](=[O])[O-])[N]([H])[*:1])[c]([H])[c]1[H]',
    'CVal': '[H][N]([C@]([H])([C](=[O])[O-])[C]([H])([C]([H])([H])[H])[C]([H])([H])[H])[*:1]',
    'NAla': '[H][N+]([C@]([H])([C](=[O])[*:2])[C]([H])([H])[H])([H])[H]',
    'NArg': '[H][N]([H])[C](=[N+]([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H])[N]([H])[H]',
    'NAsn': '[H][N]([H])[C](=[O])[C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H]',
    'NAsp': '[H][N+]([C@]([H])([C](=[O])[*:2])[C]([H])([H])[C](=[O])[O-])([H])[H]',
    'NAsh': '[H][N+]([C@]([H])([C](=[O])[*:2])[C]([H])([H])[C](=[O])[O][H])([H])[H]',
    'NCys': '[H][S][C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H]',
    'NCyx': '[S]([*:3])[C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H]',
    'NGlh': '[H][N+]([C@]([H])([C](=[O])[*:2])[C]([H])([H])[C]([H])([H])[C](=[O])[O][H])([H])[H]',
    'NGln': '[H][N]([H])[C](=[O])[C]([H])([H])[C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H]',
    'NGlu': '[H][N+]([C@]([H])([C](=[O])[*:2])[C]([H])([H])[C]([H])([H])[C](=[O])[O-])([H])[H]',
    'NGly': '[H][N+]([C]([H])([H])[C](=[O])[*:2])([H])[H]',
    'NHid': '[H][c]1[n]([H])[c]([C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H])[c]([H])[n]1',
    'NHie': '[H][c]1[n][c]([C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H])[c]([H])[n]1[H]',
    'NHip': '[H][c]1[n+]([H])[c]([C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H])[c]([H])[n]1[H]',
    'NIle': '[H][N+]([C@]([H])([C](=[O])[*:2])[C@@]([H])([C]([H])([H])[H])[C]([H])([H])[C]([H])([H])[H])([H])[H]',
    'NLeu': '[H][N+]([C@]([H])([C](=[O])[*:2])[C]([H])([H])[C]([H])([C]([H])([H])[H])[C]([H])([H])[H])([H])[H]',
    'NLyn': '[H][N+]([C@]([H])([C](=[O])[*:2])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[N]([H])[H])([H])[H]',
    'NLys': '[H][N+]([C@]([H])([C](=[O])[*:2])[C]([H])([H])[C]([H])([H])[C]([H])([H])[C]([H])([H])[N+]([H])([H])[H])([H])[H]',
    'NMet': '[H][N+]([C@]([H])([C](=[O])[*:2])[C]([H])([H])[C]([H])([H])[S][C]([H])([H])[H])([H])[H]',
    'NPhe': '[H][c]1[c]([H])[c]([H])[c]([C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H])[c]([H])[c]1[H]',
    'NPro': '[H][C]1([H])[N+]([H])([H])[C@]([H])([C](=[O])[*:2])[C]([H])([H])[C]1([H])[H]',
    'NSer': '[H][O][C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H]',
    'NThr': '[H][O][C@]([H])([C]([H])([H])[H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H]',
    'NTrp': '[H][c]1[c]([H])[c]([H])[c]2[c]([c]1[H])[c]([C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H])[c]([H])[n]2[H]',
    'NTyr': '[H][O][c]1[c]([H])[c]([H])[c]([C]([H])([H])[C@@]([H])([C](=[O])[*:2])[N+]([H])([H])[H])[c]([H])[c]1[H]',
    'NVal': '[H][N+]([C@]([H])([C](=[O])[*:2])[C]([H])([C]([H])([H])[H])[C]([H])([H])[H])([H])[H]',
}


neutral_n_term = Cap.from_smiles('[H][*:8]', name = 'H')
neutral_c_term = Cap.from_smiles("[*:9][O][H]", name = 'OH')


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
    'terminus',
    type=click.STRING,
)
@click.argument(
    'charge_state',
    type=click.STRING,
)
def main(a1_index, a2_index, a3_index, terminus, charge_state):

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

    if terminus == 'N':
        flanking_index_order = [a2_index, a3_index]
    elif terminus == 'C':
        flanking_index_order = [a2_index, a1_index]
    else:
        raise ValueError('"terminus" argument must be "N" or "C"')

    for flanking_index in flanking_index_order:
        if flanking_index == 0:
            flanking_indices.append(range(1, len(INDEX_TO_AA.keys()) + 1))
        else:
            flanking_indices.append([flanking_index])

    # Build molecules using polymetrizer, starting with the terminal residue.
    terminal_name = INDEX_TO_AA[a1_index if terminus == 'N' else a3_index]

    if charge_state == 'charged':

        # Build charged terminus, i.e. NH3+ or COO-
        terminal_residue = Monomer.from_smiles(
            terminal_smiles[f'{terminus}{terminal_name}'], name=terminal_name
        )

    elif charge_state == 'neutral':

        # Build neutral terminus, i.e. NH2 or COOH
        terminal_monomer = Monomer.from_smiles(
            AMINO_ACIDS_ALL_STATES[terminal_name], name=terminal_name
        )

        # Add neutral cap to main chain residue, i.e. H for N terminus and OH
        # for C terminus
        if terminus == 'N':

            terminal_residue = terminal_monomer.substitute(
                neutral_n_term, r_self=1, r_other=8
            )

        else:

            terminal_residue = terminal_monomer.substitute(
                neutral_c_term, r_self=2, r_other=9
            )

    else:

        raise ValueError(
            '"charge_state" argument must be "charged" or "neutral"'
        )

    if terminal_name == 'Cyx':

        monomer = terminal_residue.substitute(
            disulfide_linkage, r_self=3, r_other=3
        )

    else:

        monomer = terminal_residue

    # Build toward Ace- or Nme-capped terminus
    for flanking_index_1 in flanking_indices[0]:

        flanking_name_1 = INDEX_TO_AA[flanking_index_1]

        if terminus == 'N':

            dimer = extend_c_terminus(
                monomer, flanking_name_1, disulfide_linkage
            )

        else:

            dimer = extend_n_terminus(
                monomer, flanking_name_1, disulfide_linkage
            )

        for flanking_index_2 in flanking_indices[1]:

            flanking_name_2 = INDEX_TO_AA[flanking_index_2]

            # Check for existing file to allow batch re-submits
            if charge_state == 'charged':
                if terminus == 'N':

                    mol_name = (
                        f'N{terminal_name}_N{terminal_name}-{flanking_name_1}'
                        f'-{flanking_name_2}-Nme'
                    )

                else:

                    mol_name = (
                        f'C{terminal_name}_Ace-{flanking_name_2}'
                        f'-{flanking_name_1}-C{terminal_name}'
                    )

            else:
                if terminus == 'N':

                    mol_name = (
                        f'N0{terminal_name}_N0{terminal_name}-{flanking_name_1}'
                        f'-{flanking_name_2}-Nme'
                    )

                else:

                    mol_name = (
                        f'C0{terminal_name}_Ace-{flanking_name_2}'
                        f'-{flanking_name_1}-C0{terminal_name}'
                    )


            file_name = os.path.join(
                'aa_library_charges', f'{mol_name}.dat'
            )

            if os.path.isfile(file_name):
                continue

            if terminus == 'N':

                trimer = extend_c_terminus(
                    dimer, flanking_name_2, disulfide_linkage
                )
                trimer.substitute(nme, r_self=2, r_other=1, inplace=True)

            else:

                trimer = extend_n_terminus(
                    dimer, flanking_name_2, disulfide_linkage
                )
                trimer.substitute(ace, r_self=1, r_other=2, inplace=True)

            # Compute charges
            compute_oligomer_charges(trimer, mol_name, file_name)


if __name__ == "__main__":
    main()

