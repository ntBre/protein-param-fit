from pathlib import Path

import click
from openff.toolkit import ForceField
from openff.units import unit


@click.command()
@click.option(
    "-i",
    "--input-force-field",
    "input_force_field",
    default=Path("..", "null-nagl", "initial-force-field.offxml"),
    show_default=True,
    type=click.STRING,
    help="File path to initial force field with protein library charges.",
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
@click.option(
    "-r",
    "--reference-force-field",
    "reference_force_field",
    default=Path("..", "null-nagl", "initial-force-field.offxml"),
    show_default=True,
    type=click.STRING,
    help="Reference force field for periodicities and initial parameter values.",
)
def main(input_force_field, output_force_field, reference_force_field):
    # Initial force field with protein library charges
    force_field = ForceField(input_force_field)

    # Reference force field from which periodicities and initial parameter
    # values are taken.
    reference_force_field = ForceField(reference_force_field)

    # Define new protein-specific torsion types by label, protein-specific
    # SMIRKS, and ID of reference parameter in order from general to specific
    protein_torsion_types = [
        {
            "label": "Protein-phi-general",
            "reference_id": "t66",
            "smirks": (
                "[#6X4]-[#6X3:1](=[#8])-[#7X3:2]-[#6X4:3]-[#6X3:4](=[#8])"
                "-[#7X3]-[#6X4]"
            ),
        },
        {
            "label": "Protein-phi-sidechain",
            "reference_id": "t67",
            "smirks": (
                "[#6X4]-[#6X3:1](=[#8])-[#7X3:2]-[#6X4:3](-[#6X3](=[#8])-[#7X3]"
                "-[#6X4])-[#6X4:4]"
            ),
        },
        {
            "label": "Protein-phi-beta-branched",
            "reference_id": "t67",
            "smirks": (
                "[#6X4]-[#6X3:1](=[#8])-[#7X3:2]-[#6X4:3](-[#6X3](=[#8])-[#7X3]"
                "-[#6X4])-[#6X4:4](-[!X1])-[!X1]"
            ),
        },
        {
            "label": "Protein-phi-ring",
            "reference_id": "t67",
            "smirks": (
                "[#6X4]-[#6X3:1](=[#8])-[#7X3r5:2]-[#6X4r5:3](-[#6X3](=[#8])"
                "-[#7X3]-[#6X4])-[#6X4r5:4]"
            ),
        },
        {
            "label": "Protein-psi-general",
            "reference_id": "t22",
            "smirks": (
                "[#6X4]-[#6X3](=[#8])-[#7X3:1]-[#6X4:2]-[#6X3:3](=[#8])"
                "-[#7X3:4]-[#6X4]"
            ),
        },
        {
            "label": "Protein-psi-sidechain",
            "reference_id": "t23",
            "smirks": (
                "[#6X4]-[#6X3](=[#8])-[#7X3]-[#6X4:2](-[#6X3:3](=[#8])-[#7X3:4]"
                "-[#6X4])-[#6X4:1]"
            ),
        },
        {
            "label": "Protein-psi-beta-branched",
            "reference_id": "t23",
            "smirks": (
                "[#6X4]-[#6X3](=[#8])-[#7X3]-[#6X4:2](-[#6X3:3](=[#8])-[#7X3:4]"
                "-[#6X4])-[#6X4:1](-[!X1])-[!X1]"
            ),
        },
        # For chi1 we want long aliphatic sidechains (ARG, LYN, LYS, PRO) to get
        # the general t1 parameter, so carve out that exception in the general
        # chi1 SMIRKS. For canonical amino acids, this general parameter should
        # hit GLH, GLN, GLU, and MET N-CA-CB-CG.
        {
            "label": "Protein-chi1-aliphatic",
            "reference_id": "t1",
            "smirks": (
                "[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3]-[#6X4:4]" "-[!#1;!#6,!X4]"
            ),
        },
        # SER N-CA-CB-OG
        {
            "label": "Protein-chi1-gamma-hydroxyl",
            "reference_id": "t1",
            "smirks": "[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3]-[#8X2H1:4]",
        },
        # CYS N-CA-CB-SG
        {
            "label": "Protein-chi1-gamma-thiol",
            "reference_id": "t1",
            "smirks": "[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3]-[#16X2H1:4]",
        },
        # CYX N-CA-CB-SG
        {
            "label": "Protein-chi1-gamma-disulfide",
            "reference_id": "t1",
            "smirks": ("[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3]-[#16X2:4]-[#16X2]"),
        },
        # LEU N-CA-CB-CG
        {
            "label": "Protein-chi1-gamma-branched",
            "reference_id": "t1",
            "smirks": (
                "[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3]-[#6X4:4](-[!X1])" "-[!X1]"
            ),
        },
        # ASH and ASN N-CA-CB-CG
        {
            "label": "Protein-chi1-gamma-carbonyl",
            "reference_id": "t1",
            "smirks": "[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3]-[#6X3:4]=[#8]",
        },
        # ASP N-CA-CB-CG
        {
            "label": "Protein-chi1-gamma-carboxylate",
            "reference_id": "t1",
            "smirks": (
                "[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3]-[#6X3:4](~[#8X1])" "~[#8X1]"
            ),
        },
        # PHE and TYR N-CA-CB-CG
        {
            "label": "Protein-chi1-gamma-aromatic",
            "reference_id": "t1",
            "smirks": "[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3]-[cr6:4]",
        },
        # HID, HIE, and HIP N-CA-CB-CG
        {
            "label": "Protein-chi1-gamma-imidazole",
            "reference_id": "t1",
            "smirks": (
                "[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3]-[#6X3r5:4]~[#7r5]"
                "~[#6r5]~[#7r5]~[#6r5]"
            ),
        },
        # TRP N-CA-CB-CG
        {
            "label": "Protein-chi1-gamma-indole",
            "reference_id": "t1",
            "smirks": (
                "[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3]-[#6X3r5:4]~[#6x3]"
                "~[#6x3]~[#7r5]~[#6r5]"
            ),
        },
        # THR N-CA-CB-OG1
        {
            "label": "Protein-chi1-beta-branched-hydroxyl",
            "reference_id": "t1",
            "smirks": ("[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3](-[!X1])-[#8X2H1:4]"),
        },
        # ILE N-CA-CB-CG2, THR N-CA-CB-CG2, VAL N-CA-CB-CG1, VAL N-CA-CB-CG2
        {
            "label": "Protein-chi1-beta-branched-methyl",
            "reference_id": "t1",
            "smirks": ("[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3](-[!X1])-[#6X4H3:4]"),
        },
        # ILE N-CA-CB-CG1
        {
            "label": "Protein-chi1-beta-branched-methylene",
            "reference_id": "t1",
            "smirks": (
                "[#7X3:1]-[#6X4:2](-[#6X3]=[#8])-[#6X4:3](-[!X1])-[#6X4:4]" "-[!X1]"
            ),
        },
        # For chi2, long aliphatic sidechains (ARG, LYN, LYS, PRO) will get the
        # general t2 parameter.
        # MET CA-CB-CG-SD
        {
            "label": "Protein-chi2-delta-thioether",
            "reference_id": "t1",
            "smirks": ("[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[#6X4:3]-[#16X2:4]"),
        },
        # GLH and GLN CA-CB-CG-CD
        {
            "label": "Protein-chi2-delta-carbonyl",
            "reference_id": "t1",
            "smirks": ("[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[#6X4:3]-[#6X3:4]=[#8]"),
        },
        # GLU CA-CB-CG-CD
        {
            "label": "Protein-chi2-delta-carboxylate",
            "reference_id": "t1",
            "smirks": (
                "[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[#6X4:3]-[#6X3:4]"
                "(~[#8X1])~[#8X1]"
            ),
        },
        # LEU CA-CB-CG-CD1 and CA-CB-CG-CD2
        {
            "label": "Protein-chi2-gamma-branched-methyl",
            "reference_id": "t2",
            "smirks": (
                "[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[#6X4:3](-[!X1])" "-[#6X4H3:4]"
            ),
        },
        # For ASH and ASN, CA-CB-CG-OD1 (carbonyl oxygen) should get the general
        # t18 parameter, and specific parameters should be fit for the ASH
        # hydroxyl oxygen and the ASN amide nitrogen.
        # ASH CA-CB-CG-OD2
        {
            "label": "Protein-chi2-gamma-carboxyl",
            "reference_id": "t17",
            "smirks": (
                "[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[#6X3:3](=[#8])" "-[#8X2H1:4]"
            ),
        },
        # ASN CA-CB-CG-ND2
        {
            "label": "Protein-chi2-gamma-amide",
            "reference_id": "t23",
            "smirks": (
                "[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[#6X3:3](=[#8])" "-[#7X3:4]"
            ),
        },
        # ASP CA-CB-CG-CD
        {
            "label": "Protein-chi2-gamma-carboxylate",
            "reference_id": "t18a",
            "smirks": (
                "[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[#6X3:3](~[#8X1])" "~[#8X1:4]"
            ),
        },
        # PHE and TYR CA-CB-CG-CD
        {
            "label": "Protein-chi2-gamma-aromatic",
            "reference_id": "t17",
            "smirks": "[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[cr6:3]:[cr6:4]",
        },
        # For His, HIE and HIP CA-CB-CG-CD2 should get the general t18
        # parameter, and specific parameters should be fit for HID CD2, HIE ND1,
        # and HID and HIP ND1
        # HID CA-CB-CG-CD2
        {
            "label": "Protein-chi2-gamma-imidazole-delta-tautomer",
            "reference_id": "t18",
            "smirks": (
                "[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[#6X3r5:3]~[#6r5:4]"
                "~[#7H0r5]~[#6r5]~[#7H1r5]"
            ),
        },
        # HIE CA-CB-CG-ND1
        {
            "label": "Protein-chi2-gamma-imidazole-epsilon-tautomer",
            "reference_id": "t17",
            "smirks": (
                "[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[#6X3r5:3]~[#7H0r5:4]"
                "~[#6r5]~[#7H1r5]~[#6r5]"
            ),
        },
        # HID and HIP CA-CB-CG-ND1
        {
            "label": "Protein-chi2-gamma-imidazolium",
            "reference_id": "t23",
            "smirks": (
                "[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[#6X3r5:3]~[#7H1r5:4]"
                "~[#6r5]~[#7r5]~[#6r5]"
            ),
        },
        # For TRP, CA-CB-CG-CD1 (bonded to nitrogen NE1) should get the general
        # t18 parameter, and a specific parameter should be fit for the ring
        # junction carbon.
        # TRP CA-CB-CG-CD2
        {
            "label": "Protein-chi2-gamma-indole",
            "reference_id": "t17",
            "smirks": (
                "[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2]-[#6X3r5:3]~[#6x3:4]"
                "~[#6x3]~[#7r5]~[#6r5]"
            ),
        },
        # ILE CA-CB-CG1-CD1
        {
            "label": "Protein-chi2-beta-branched-ethyl",
            "reference_id": "t2",
            "smirks": (
                "[#7X3]-[#6X4:1](-[#6X3]=[#8])-[#6X4:2](-[!X1])-[#6X4:3]" "-[#6X4H3:4]"
            ),
        },
    ]

    # Add protein-specific torsion types to force field at end of proper
    # torsions so that these take priority over small molecule parameters
    torsion_handler = force_field["ProperTorsions"]
    reference_torsion_handler = reference_force_field["ProperTorsions"]

    for torsion_type in protein_torsion_types:
        # Copy parameters from an existing torsion type
        reference_parameters = reference_torsion_handler.get_parameter(
            {"id": torsion_type["reference_id"]}
        )

        if len(reference_parameters) == 0:
            raise ValueError(
                "Cannot find reference parameter " f'{torsion_type["reference_id"]}'
            )

        else:
            reference_parameter = reference_parameters[0]

        new_parameter = {
            "id": torsion_type["label"],
            "smirks": torsion_type["smirks"],
        }

        # First write k = 0.0 and phase = 0.0 up to periodicity 4
        for i in range(1, 5):
            new_parameter[f"periodicity{i}"] = i
            new_parameter[f"k{i}"] = 0.0 / unit.mole * unit.kilocalorie
            new_parameter[f"phase{i}"] = 0.0 * unit.degree
            new_parameter[f"idivf{i}"] = 1.0

        # Then overwrite periodicities present in the reference parameter
        for i, n in enumerate(reference_parameter.periodicity):
            new_parameter[f"periodicity{n}"] = n
            new_parameter[f"k{n}"] = reference_parameter.k[i]
            new_parameter[f"phase{n}"] = reference_parameter.phase[i]
            new_parameter[f"idivf{n}"] = reference_parameter.idivf[i]

        torsion_handler.add_parameter(new_parameter)

    force_field.to_file(output_force_field)


if __name__ == "__main__":
    main()
