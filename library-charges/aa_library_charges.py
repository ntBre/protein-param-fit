import click
import numpy
from openff.toolkit.typing.engines.smirnoff import ForceField
import os
from polymetrizer import Monomer
from polymetrizer.cap import Cap
from polymetrizer.tests.smiles import AMINO_ACIDS_ALL_STATES, ACE, NME
from rdkit.Chem import AllChem


INDEX_TO_AA = {
    1: 'Ala',
    2: 'Arg',
    3: 'Ash',
    4: 'Asn',
    5: 'Asp',
    6: 'Cys',
    7: 'Cyx',
    8: 'Glh',
    9: 'Gln',
    10: 'Glu',
    11: 'Gly',
    12: 'Hid',
    13: 'Hie',
    14: 'Hip',
    15: 'Ile',
    16: 'Leu',
    17: 'Lyn',
    18: 'Lys',
    19: 'Met',
    20: 'Phe',
    21: 'Pro',
    22: 'Ser',
    23: 'Thr',
    24: 'Trp',
    25: 'Tyr',
    26: 'Val',
}


ace = Cap.from_smiles(ACE, name="Ace")
nme = Cap.from_smiles(NME, name="Nme")

# Compute ELF10 charges for an oligomer
def compute_oligomer_charges(oligomer, mol_name, file_name):

    # OpenFF territory
    # Due to the ordered nature of the oligomer graph,
    # we know that the central atoms are first
    offmol = oligomer.to_openff()

    # Optimize to reduce error in AM1 optimization by OpenEye as this
    # optimization is atom-order dependent might need to comment this out if
    # there are stereo errors
    try:

        offmol.generate_conformers(n_conformers=500)  # rec_confs

        rdmol = offmol.to_rdkit()
        AllChem.MMFFOptimizeMoleculeConfs(rdmol, numThreads=0, maxIters=1000)
        opt = type(offmol).from_rdkit(rdmol, allow_undefined_stereo=True)

        offmol._conformers = []
        for conformer in opt._conformers:
            offmol._add_conformer(conformer)

    except:
        pass

    # Assign charges and save
    offmol.assign_partial_charges('am1bccelf10')

    smiles_file = os.path.join('aa_library_charges', f'{mol_name}.smi')
    offmol.properties.pop('atom_map', None)
    smiles = offmol.to_smiles(mapped=True)

    with open(smiles_file, 'w') as f:
        f.write(f'{mol_name} {smiles}')

    charges = offmol.partial_charges._value
    numpy.savetxt(file_name, charges)
    print(f'Saved to {file_name}')


# Extend a polypeptide chain at the N terminus, possibly add a disulfide linkage
def extend_n_terminus(oligomer, monomer_name, disulfide_linkage=None):

    monomer = Monomer.from_smiles(
        AMINO_ACIDS_ALL_STATES[monomer_name], name=monomer_name
    )

    new_oligomer = oligomer.substitute(monomer, r_self=1, r_other=2)

    if monomer_name == 'Cyx':
        new_oligomer.substitute(
            disulfide_linkage, r_self=3, r_other=3, inplace=True
        )

    return new_oligomer


# Extend a polypeptide chain at the C terminus, possibly add a disulfide linkage
def extend_c_terminus(oligomer, monomer_name, disulfide_linkage=None):

    monomer = Monomer.from_smiles(
        AMINO_ACIDS_ALL_STATES[monomer_name], name=monomer_name
    )

    new_oligomer = oligomer.substitute(monomer, r_self=2, r_other=1)

    if monomer_name == 'Cyx':
        new_oligomer.substitute(
            disulfide_linkage, r_self=3, r_other=3, inplace=True
        )

    return new_oligomer


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
    'a4_index',
    type=click.INT,
)
@click.argument(
    'a5_index',
    type=click.INT,
)
def main(a1_index, a2_index, a3_index, a4_index, a5_index):

    # Build Ace-Val-Cyx-Val-Nme disulfide linkage
    build_disulfide_linkage = False
    for aa_index in [a1_index, a2_index, a3_index, a4_index, a5_index]:
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
    for flanking_index in [a1_index, a2_index, a4_index, a5_index]:
        if flanking_index == 0:
            flanking_indices.append(range(1, len(INDEX_TO_AA.keys()) + 1))
        else:
            flanking_indices.append([flanking_index])

    # Build molecule with polymetrizer, starting with central residue.
    a3_name = INDEX_TO_AA[a3_index]
    a3 = Monomer.from_smiles(AMINO_ACIDS_ALL_STATES[a3_name], name=a3_name)

    if a3_name == 'Cyx':
        monomer = a3.substitute(disulfide_linkage, r_self=3, r_other=3)
    else:
        monomer = a3

    # Build toward N terminus.
    for flanking_index_2 in flanking_indices[1]:

        a2_name = INDEX_TO_AA[flanking_index_2]
        dimer = extend_n_terminus(monomer, a2_name, disulfide_linkage)

        for flanking_index_1 in flanking_indices[0]:

            a1_name = INDEX_TO_AA[flanking_index_1]
            trimer = extend_n_terminus(dimer, a1_name, disulfide_linkage)

            # N terminal cap
            trimer.substitute(ace, r_self=1, r_other=2, inplace=True)

            # Build toward C terminus
            for flanking_index_4 in flanking_indices[2]:

                a4_name = INDEX_TO_AA[flanking_index_4]
                tetramer = extend_c_terminus(trimer, a4_name, disulfide_linkage)

                # C terminal residue
                for flanking_index_5 in flanking_indices[3]:

                    a5_name = INDEX_TO_AA[flanking_index_5]

                    # Check for existing file to allow batch re-submits
                    mol_name = (
                        f'{a3_name}_Ace-{a1_name}-{a2_name}-{a3_name}-{a4_name}'
                        f'-{a5_name}-Nme'
                    )
                    file_name = os.path.join(
                        'aa_library_charges', f'{mol_name}.dat'
                    )

                    if os.path.isfile(file_name):
                        continue

                    pentamer = extend_c_terminus(
                        tetramer, a5_name, disulfide_linkage
                    )

                    # C terminal cap
                    pentamer.substitute(nme, r_self=2, r_other=1, inplace=True)

                    # Compute charges
                    compute_oligomer_charges(pentamer, mol_name, file_name)


if __name__ == "__main__":
    main()

