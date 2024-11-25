import os
import json
from pathlib import Path

import click
from openff.bespokefit.optimizers.forcebalance import ForceBalanceInputFactory
from openff.bespokefit.schema.fitting import OptimizationSchema, OptimizationStageSchema
from openff.bespokefit.schema.optimizers import ForceBalanceSchema
from openff.bespokefit.schema.smirnoff import (
    AngleHyperparameters,
    AngleSMIRKS,
    BondHyperparameters,
    ImproperTorsionHyperparameters,
    ProperTorsionHyperparameters,
    BondSMIRKS,
    ProperTorsionSMIRKS,
)
from openff.bespokefit.schema.targets import (
    AbInitioTargetSchema,
    OptGeoTargetSchema,
    TorsionProfileTargetSchema,
)
from openff.qcsubmit.results import (
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)
from openff.qcsubmit.results.filters import SMARTSFilter, SMILESFilter
from openff.toolkit import ForceField


def load_training_data(
    optimization_dataset: str,
    torsiondrive_dataset: str,
    smarts_to_exclude: str | None = None,
    smiles_to_exclude: str | None = None,
    verbose: bool = False
):
    if smarts_to_exclude is not None:
        exclude_smarts = Path(smarts_to_exclude).read_text().splitlines()
    else:
        exclude_smarts = []
    
    if smiles_to_exclude is not None:
        exclude_smiles = Path(smiles_to_exclude).read_text().splitlines()
    else:
        exclude_smiles = []

    torsion_training_set = TorsionDriveResultCollection.parse_file(torsiondrive_dataset)
    if verbose:
        print(f"Loaded torsion training set with {torsion_training_set.n_results} entries.")
    
    torsion_training_set = torsion_training_set.filter(
        SMARTSFilter(smarts_to_exclude=exclude_smarts),
        SMILESFilter(smiles_to_exclude=exclude_smiles),
    )
    
    if verbose:
        print(f"Filtered torsion training set to {torsion_training_set.n_results} entries.")
    
    optimization_training_set = OptimizationResultCollection.parse_file(optimization_dataset)
    if verbose:
        print(f"Loaded optimization training set with {optimization_training_set.n_results} entries.")
    optimization_training_set = optimization_training_set.filter(
        SMARTSFilter(smarts_to_exclude=exclude_smarts),
        SMILESFilter(smiles_to_exclude=exclude_smiles),
    )
    if verbose:
        print(f"Filtered optimization training set to {optimization_training_set.n_results} entries.")

    return torsion_training_set, optimization_training_set


@click.command()
@click.option(
    "--tag",
    type=str,
    default="forcebalance",
    help="The tag to use for the fitting run.",
)
@click.option(
    "--force-field-path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the force field to use. (offxml)",
)
@click.option(
    "--optimization-dataset-path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the optimization dataset to use. (JSON)",
)
@click.option(
    "--torsiondrive-dataset-path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the torsion dataset to use. (JSON)",
)
@click.option(
    "--valence-smirks-path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the valence parameters to optimize (JSON).",
)
@click.option(
    "--torsion-smirks-path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the torsions to optimize (JSON).",
)
@click.option(
    "--smarts-to-exclude",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    default=None,
    help=(
        "The path to a file containing a list of SMARTS patterns "
        "to exclude from the training set. "
        "The patterns should be separated by new lines."
    ),
)
@click.option(
    "--smiles-to-exclude",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    default=None,
    help=(
        "The path to a file containing a list of SMILES patterns "
        "to exclude from the training set. "
        "The patterns should be separated by new lines."
    ),
)
@click.option(
    "--protein-record-ids-path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    default=None,
    help=(
        "The path to a file containing a list of QCFractal record IDs "
        "corresponding to protein TorsionDrives. "
        "The record IDs should be separated by new lines."
    ),
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Whether to print verbose logging messages.",
)
@click.option(
    "--max-iterations",
    type=int,
    default=50,
    show_default=True,
    help="The maximum number of iterations to run the fitting for.",
)
@click.option(
    "--opt-geo-weight",
    type=float,
    default=0.01,
    show_default=True,
    help="The weight of optimized geometry targets in the ForceBalance "
        "objective function.",
)
@click.option(
    "--protein-torsiondrive-weight",
    type=float,
    default=1.0,
    show_default=True,
    help="The weight of protein TorsionDrive targets in the ForceBalance "
        "objective function.",
)
@click.option(
    "--port",
    type=int,
    default=55125,
    show_default=True,
    help="The port to run the server on.",
)
@click.option(
    "--torsiondrive-target-type",
    type=str,
    default="TorsionProfile",
    show_default=True,
    help="The ForceBalance target type for TorsionDrive QC data.",
)
@click.option(
    "--torsiondrive-weight",
    type=float,
    default=1.0,
    show_default=True,
    help="The weight of TorsionDrive targets in the ForceBalance objective "
        "function.",
)
@click.option(
    "--proper-torsion-prior",
    type=float,
    default=5.0,
    show_default=True,
    help="ProperTorsion prior for ForceBalance."
)
def main(
    tag: str,
    force_field_path: str,
    optimization_dataset_path: str,
    torsiondrive_dataset_path: str,
    valence_smirks_path: str,
    torsion_smirks_path: str,
    smarts_to_exclude: str | None,
    smiles_to_exclude: str | None,
    protein_record_ids_path: str | None,
    verbose: bool,
    max_iterations: int,
    opt_geo_weight: float,
    protein_torsiondrive_weight: float,
    port: int,
    torsiondrive_target_type: str,
    torsiondrive_weight: float,
    proper_torsion_prior: float,
):
    optimizer = ForceBalanceSchema(
        max_iterations=max_iterations,
        step_convergence_threshold=0.01,
        objective_convergence_threshold=0.1,
        gradient_convergence_threshold=0.1,
        n_criteria=2,
        initial_trust_radius=-1.0,
        finite_difference_h=0.01,
        extras={
            "wq_port": str(port),
            "asynchronous": "True",
            "search_tolerance": "0.1",
            "backup": "0",
            "retain_micro_outputs": "0",
        },
    )

    # Prepare QC datasets
    torsion_training_set, optimization_training_set = load_training_data(
        optimization_dataset=optimization_dataset_path,
        torsiondrive_dataset=torsiondrive_dataset_path,
        smarts_to_exclude=smarts_to_exclude,
        smiles_to_exclude=smiles_to_exclude,
        verbose=verbose
    )

    # Set up options for TorsionDrive QC data
    torsiondrive_extras = {"remote": "1"}
    if torsiondrive_target_type == "TorsionProfile":
        torsiondrive_target_schema_type = TorsionProfileTargetSchema
    elif torsiondrive_target_type == "AbInitio":
        torsiondrive_target_schema_type = AbInitioTargetSchema
        torsiondrive_extras["energy_asymmetry"] = 100.0
        torsiondrive_extras["energy_mode"] = "qm_minimum"
    else:
        raise ValueError(
            "Argument torsiondrive-target-type must be one of:"
            "\n    TorsionProfile\n    AbInitio"
        )

    # Set up TorsionDrive target schemas
    if protein_record_ids_path is None:
        torsion_profile_target_schemas = [
            torsiondrive_target_schema_type(
                reference_data=torsion_training_set,
                weight=torsiondrive_weight,
                attenuate_weights=True,
                energy_denominator=1.0,
                energy_cutoff=8.0,
                extras=torsiondrive_extras,
            )
        ]

    else:
        # Split TorsionDrive dataset into protein and small molecule datasets
        protein_record_ids = {
            int(record_id)
            for record_id in Path(protein_record_ids_path).read_text().splitlines()
        }
        client_address = list(torsion_training_set.entries.keys())[0]

        protein_entries = [
            entry
            for entry in torsion_training_set.entries[client_address]
            if entry.record_id in protein_record_ids
        ]
        small_molecule_entries = [
            entry
            for entry in torsion_training_set.entries[client_address]
            if entry.record_id not in protein_record_ids
        ]

        protein_torsion_training_set = TorsionDriveResultCollection(
            entries={client_address: protein_entries}
        )
        small_molecule_torsion_training_set = TorsionDriveResultCollection(
            entries={client_address: small_molecule_entries}
        )

        torsion_profile_target_schemas = [
            torsiondrive_target_schema_type(
                reference_data=protein_torsion_training_set,
                weight=protein_torsiondrive_weight,
                attenuate_weights=True,
                energy_denominator=1.0,
                energy_cutoff=8.0,
                extras=torsiondrive_extras,
            ),
            torsiondrive_target_schema_type(
                reference_data=small_molecule_torsion_training_set,
                weight=torsiondrive_weight,
                attenuate_weights=True,
                energy_denominator=1.0,
                energy_cutoff=8.0,
                extras=torsiondrive_extras,
            )
        ]

    targets = [
        *torsion_profile_target_schemas,
        OptGeoTargetSchema(
            reference_data=optimization_training_set,
            weight=opt_geo_weight,
            extras={"batch_size": 30, "remote": "1"},
            bond_denominator=0.05,
            angle_denominator=5.0,
            dihedral_denominator=10.0,
            improper_denominator=10.0,
        ),
    ]

    # a16, a17, a27, a35
    linear_angle_smirks = [
        "[*:1]~[#6X2:2]~[*:3]",  # a16
        "[*:1]~[#7X2:2]~[*:3]",  # a17
        "[*:1]~[#7X2:2]~[#7X1:3]",  # a27
        "[*:1]=[#16X2:2]=[*:3]",
    ]  # a35, this one anyways doesn't have a training target for ages

    with open(valence_smirks_path, "r") as json_file:
        valence_smirks = json.load(json_file)
    with open(torsion_smirks_path, "r") as json_file:
        torsion_smirks = json.load(json_file)

    target_parameters = []
    for smirks in valence_smirks["Angles"]:
        if smirks in linear_angle_smirks:
            parameter = AngleSMIRKS(smirks=smirks, attributes={"k"})
        else:
            parameter = AngleSMIRKS(smirks=smirks, attributes={"k", "angle"})
        target_parameters.append(parameter)
    
    for smirks in valence_smirks["Bonds"]:
        target_parameters.append(BondSMIRKS(smirks=smirks, attributes={"k", "length"}))
    
    force_field = ForceField(force_field_path)

    torsion_handler = force_field.get_parameter_handler("ProperTorsions")
    for smirks in torsion_smirks["ProperTorsions"]:
        original_k = torsion_handler.parameters[smirks].k
        attributes = {f"k{i + 1}" for i in range(len(original_k))}
        target_parameters.append(ProperTorsionSMIRKS(smirks=smirks, attributes=attributes))

    optimization_schema = OptimizationSchema(
        id=tag,
        initial_force_field=str(Path(force_field_path).resolve()),
        stages=[
            OptimizationStageSchema(
                optimizer=optimizer,
                targets=targets,
                parameters=target_parameters,
                parameter_hyperparameters=[
                    AngleHyperparameters(priors={"k": 100, "angle": 5}),
                    BondHyperparameters(priors={"k": 100, "length": 0.1}),
                    ProperTorsionHyperparameters(priors={"k": proper_torsion_prior}),
                    ImproperTorsionHyperparameters(priors={"k": 5}),
                ],
            )
        ]
    )

    with open(f"{optimization_schema.id}.json", "w") as json_file:
        json_file.write(optimization_schema.json(indent=2))

    # Generate the ForceBalance inputs
    ForceBalanceInputFactory.generate(
        optimization_schema.id,
        optimization_schema.stages[0],
        ForceField(optimization_schema.initial_force_field),
    )



if __name__ == "__main__":
    main()
