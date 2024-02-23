import click
from openff.toolkit.typing.engines.smirnoff import ForceField
from openmm import unit
from pathlib import Path
from residue_smirks import RESIDUE_SMIRKS


@click.command()
@click.option(
    '-l',
    '--library_charges',
    default = Path('mean-charges', 'mean-charges.dat'),
    show_default = True,
    type = click.STRING,
    help = 'File path to the protein library charges.',
)
@click.option(
    '-n',
    '--nh2-library_charges',
    default = Path('mean-charges', 'nh2-mean-charges.dat'),
    show_default = True,
    type = click.STRING,
    help = 'File path to the NH2 cap library charges.',
)
@click.option(
    '-o',
    '--output_ff',
    default = 'protein-library-charges.offxml',
    show_default = True,
    type = click.STRING,
    help = 'File path to which the output force field will be written.',
)
def main(library_charges, nh2_library_charges, output_ff):

    nh2_smirks = '[#7X3]-[#6X4]-[#6X3](=[#8])-[#7:1](-[#1:2])-[#1:3]'

    # Dictionary of SMIRNOFF LibraryCharge parameter dictionaries, indexed by
    # residue name
    library_charge_dicts = {}

    # List of residues in desired order of general to specific
    residues = []

    # Read library charges
    with open(library_charges, 'r') as library_charge_file:
        for line in library_charge_file:

            residue, atom, charge = line.split()
            residue = residue.upper()

            if residue not in residues:

                # Add NH2 cap immediately before NME cap
                if residue == 'NME':
                    with open(nh2_library_charges, 'r') as nh2_library_charge_file:
                        for nh2_line in nh2_library_charge_file:

                            nh2_residue, nh2_atom, nh2_charge = nh2_line.split()

                            if nh2_residue != 'Nh2':
                                continue

                            nh2_residue = 'NH2'

                            if nh2_residue not in residues:

                                residues.append(nh2_residue)
                                charge_idx = 0
                                library_charge_dicts[nh2_residue] = {
                                    'smirks': nh2_smirks,
                                    'id': f'Protein-{nh2_residue}',
                                }

                            charge_idx += 1
                            library_charge_dicts[nh2_residue][f'charge{charge_idx}'] = (
                                nh2_charge * unit.elementary_charge
                            )

                smirks = RESIDUE_SMIRKS[residue]

                # Replace element symbols with atomic number
                smirks = smirks.replace('H', '#1')
                smirks = smirks.replace('C', '#6')
                smirks = smirks.replace('N', '#7')
                smirks = smirks.replace('O', '#8')
                smirks = smirks.replace('S', '#16')

                # Add non-tagged peptide bond atoms to non-terminal residue
                # boundaries to avoid assigning charges to non-peptides
                if (
                    not (len(residue) > 3 and residue[0] == 'N')
                    and residue != 'ACE'
                ):

                    smirks = '[#6X4]-[#6X3](=[#8])-' + smirks

                if (
                    not (len(residue) > 3 and residue[0] == 'C')
                    and residue != 'NME'
                ):

                    smirks = smirks + '-[#7X3;$(*-[#6X4]),H2]'

                residues.append(residue)
                charge_idx = 0
                library_charge_dicts[residue] = {
                    'smirks': smirks,
                    'id': f'Protein-{residue}',
                }

            charge_idx += 1
            library_charge_dicts[residue][f'charge{charge_idx}'] = (
                charge * unit.elementary_charge
            )

    # Initialize the force field
    force_field = ForceField()
    force_field.aromaticity_model = 'OEAroModel_MDL'

    # Initialize empty parameter handlers
    for smirnoff_tag in [
        'Bonds', 'Angles', 'ProperTorsions', 'ImproperTorsions', 'vdW',
        'Electrostatics'
    ]:

        handler = force_field.get_parameter_handler(smirnoff_tag)

    # Add library charge parameters
    handler = force_field.get_parameter_handler('LibraryCharges')
    for residue in residues:
        handler.add_parameter(library_charge_dicts[residue])

    force_field.to_file(output_ff)


if __name__ == "__main__":
    main()

