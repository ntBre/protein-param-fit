import click
import json
from openff.bespokefit.optimizers.forcebalance import ForceBalanceInputFactory
from openff.bespokefit.schema.fitting import (
    OptimizationSchema, OptimizationStageSchema
)
from openff.bespokefit.schema.optimizers import ForceBalanceSchema
from openff.bespokefit.schema.smirnoff import (
    AngleHyperparameters, AngleSMIRKS, BondHyperparameters, BondSMIRKS,
    ProperTorsionHyperparameters, ProperTorsionSMIRKS
)
from openff.bespokefit.schema.targets import (
    OptGeoTargetSchema, TorsionProfileTargetSchema
)
from openff.qcsubmit.results import (
    OptimizationResultCollection, TorsionDriveResultCollection
)
from openff.toolkit.typing.engines.smirnoff import ForceField
from pathlib import Path


@click.command()
@click.option(
    '-i',
    '--input_ff',
    default = 'initial-force-field.offxml',
    show_default = True,
    type = click.STRING,
    help = 'File path to the initial force field (small molecule + protein '
        'library charges).',
)
@click.option(
    '-s',
    '--smirks_dir',
    default = 'training-smirks',
    show_default = True,
    type = click.STRING,
    help = 'Directory path to training SMIRKS.',
)
def main(input_ff, smirks_dir):

    initial_force_field = ForceField(input_ff)

    optimization_training_set = OptimizationResultCollection.parse_file(
        Path(smirks_dir, 'optimization-smarts-filter-dataset.json')
    )

    torsion_training_set = dict()

    # Order labels so that 2-D protein datasets come before 1-D small molecule
    # datasets
    for label in [
        'protein-weight1', 'protein-weight2', 'protein-weight3',
        'protein-weight4', 'protein-weight6', 'small-molecule'
    ]:

        torsion_training_set[label] = TorsionDriveResultCollection.parse_file(
            Path(smirks_dir, f'{label}-torsiondrive-smarts-filter-dataset.json')
        )

    # Set weights for TorsionDrive datasets
    n_results = {
        label: dataset.n_results
        for label, dataset in torsion_training_set.items()
    }

    protein_degeneracy = {
        label: float(label[-1]) for label in n_results
        if label != 'small-molecule'
    }
    total_protein_weight = sum(
        [
            n_results[label] / protein_degeneracy[label]
            for label in protein_degeneracy
        ]
    )
    sm_to_protein = n_results['small-molecule'] / total_protein_weight

    torsion_weights = {
        label: sm_to_protein / protein_degeneracy[label]
        for label in protein_degeneracy
    }
    torsion_weights['small-molecule'] = 1.0

    for label, weight in torsion_weights.items():
        print(f'{label:15s} {n_results[label]:4d} {weight:6.2f}')

    # List of ForceBalance targets
    forcebalance_targets = [
        TorsionProfileTargetSchema(
            reference_data = torsion_training_set[label],
            weight = torsion_weights[label],
            energy_denominator = 5.0,
            energy_cutoff = 20.0,
            extras = {'remote': '1'},
        )
        for label in torsion_training_set
    ]
    forcebalance_targets.append(
        OptGeoTargetSchema(
            reference_data = optimization_training_set,
            weight = 0.1,
            #extras = {'batch_size': 30, 'remote': '1'},
            extras = {'batch_size': 1, 'remote': '1'},
        )
    )

    # Read parameters to train
    with open(Path(smirks_dir, 'optimization-bond-smirks.json'), 'r') as f:
        bond_smirks = json.load(f)
    with open(Path(smirks_dir, 'optimization-angle-smirks.json'), 'r') as f:
        angle_smirks = json.load(f)
    with open(Path(smirks_dir, 'torsiondrive-torsion-smirks.json'), 'r') as f:
        torsion_smirks = json.load(f)

    # a16, a17, a27, a35
    linear_angle_smirks = [
        '[*:1]~[#6X2:2]~[*:3]',  # a16
        '[*:1]~[#7X2:2]~[*:3]',  # a17
        '[*:1]~[#7X2:2]~[#7X1:3]',  # a27
        '[*:1]=[#16X2:2]=[*:3]',  # a35
    ]

    proper_torsion_handler = initial_force_field.get_parameter_handler(
        'ProperTorsions'
    )

    target_parameters = [
        *[
            BondSMIRKS(smirks=smirks, attributes={'k', 'length'})
            for smirks in bond_smirks['Bonds']
        ],
        *[
            AngleSMIRKS(smirks=smirks, attributes={'k', 'angle'})
            if smirks not in linear_angle_smirks
            else AngleSMIRKS(smirks=smirks, attributes={'k'})
            for smirks in angle_smirks['Angles']
        ],
        *[
            ProperTorsionSMIRKS(
                smirks=smirks,
                attributes={
                    f'k{i + 1}' for i in range(
                        len(proper_torsion_handler.parameters[smirks].k)
                    )
                },
            )
            for smirks in torsion_smirks['ProperTorsions']
        ],
    ]

    # Define the full schema for the optimization.
    optimization_schema = OptimizationSchema(
        id = 'forcebalance',
        initial_force_field = str(Path(input_ff).resolve()),
        # Define the optimizer / ForceBalance specific settings.
        stages = [
            OptimizationStageSchema(
                optimizer = ForceBalanceSchema(
                    max_iterations = 50,
                    step_convergence_threshold = 0.01,
                    objective_convergence_threshold = 0.1,
                    gradient_convergence_threshold = 0.1,
                    n_criteria = 2,
                    initial_trust_radius = -1.0,
                    extras = {
                        'asynchronous': 'True',
                        'search_tolerance': 1,
                        'wq_port': '55125',
                    },
                ),
                targets = forcebalance_targets,
                # Define the parameters to refit and the priors to place on them
                parameters = target_parameters,
                parameter_hyperparameters = [
                    AngleHyperparameters(priors = {'k': 100, 'length': 20}),
                    BondHyperparameters(priors = {'k': 100, 'length': 0.1}),
                    ProperTorsionHyperparameters(priors = {'k': 1}),
                ],
            )
        ],
    )

    with open(f'{optimization_schema.id}.json', 'w') as json_file:
        json_file.write(optimization_schema.json())

    # Generate the ForceBalance inputs
    ForceBalanceInputFactory.generate(
        optimization_schema.id,
        optimization_schema.stages[0],
        ForceField(optimization_schema.initial_force_field),
    )


if __name__ == "__main__":
    main()

