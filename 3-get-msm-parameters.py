import json
from collections import defaultdict
from pathlib import Path

import click
import numpy
import tqdm
from openff.qcsubmit.results import BasicResultCollection, OptimizationResultCollection
from openff.qcsubmit.results.filters import LowestEnergyFilter
from openff.toolkit import ForceField, Molecule
from openff.units import unit
from qcportal.optimization import OptimizationRecord
from qcportal.singlepoint import SinglepointRecord
from qubekit.bonded.mod_seminario import ModSeminario
from qubekit.molecules import Ligand


def calculate_parameters(
    qc_record: SinglepointRecord,
    molecule: Molecule,
    forcefield: ForceField,
) -> dict[str, dict[str, list[unit.Quantity]]]:
    """
    Calculate the modified seminario parameters for the given input molecule
    and store them by OFF SMIRKS.
    """

    mod_sem = ModSeminario()

    # Create the QUBEKit molecule. Atoms should be in the same order as the
    # OpenFF molecule
    qube_mol = Ligand.from_rdkit(molecule.to_rdkit(), name="offmol")
    qube_mol.hessian = qc_record.return_result

    # Calculate the Modified Seminario parameters and store in the molecule
    qube_mol = mod_sem.run(qube_mol)

    # Label assigned parameters for the OpenFF molecule
    labels = forcefield.label_molecules(molecule.to_topology())[0]

    # Loop over all bonds and angles and collect the results in OpenMM units
    all_parameters = {
        "bond_eq": defaultdict(list),
        "bond_k": defaultdict(list),
        "angle_eq": defaultdict(list),
        "angle_k": defaultdict(list),
    }

    for bond, parameter in labels["Bonds"].items():
        # bond is a tuple of the atom indinces the parameter is applied to
        qube_param = qube_mol.BondForce[bond]
        all_parameters["bond_eq"][parameter.smirks].append(qube_param.length)
        all_parameters["bond_k"][parameter.smirks].append(qube_param.k)

    for angle, parameter in labels["Angles"].items():
        qube_param = qube_mol.AngleForce[angle]
        all_parameters["angle_eq"][parameter.smirks].append(qube_param.angle)
        all_parameters["angle_k"][parameter.smirks].append(qube_param.k)

    return all_parameters


@click.command()
@click.option(
    "--initial-force-field",
    "initial_force_field",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the initial force field file (OFFXML).",
)
@click.option(
    "--output-force-field",
    "output_force_field",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the output force field file (OFFXML).",
)
@click.option(
    "--optimization-dataset",
    "optimization_dataset",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help="The path to the optimization dataset.",
)
@click.option(
    "--working-directory",
    "working_directory",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    required=False,
    help=(
        "The path to the working directory. "
        "Intermediate files are saved here if provided"
    ),
)
@click.option(
    "--verbose/--no-verbose",
    default=False,
    help="Enable verbose logging.",
)
def main(
    initial_force_field: str,
    output_force_field: str,
    optimization_dataset: str,
    working_directory: str | None = None,
    verbose: bool = False,
):
    dataset = OptimizationResultCollection.parse_file(optimization_dataset)

    # filter for lowest energy results
    filtered = dataset.filter(LowestEnergyFilter())

    # filter to only keep entries with hessians calculated
    hessian_set = filtered.to_basic_result_collection(driver="hessian")

    if working_directory is not None:
        hessian_file = Path(working_directory, "hessian-set.json")
        with open(hessian_file, "w") as json_file:
            json_file.write(hessian_set.json(indent=2))
        if verbose:
            print(f"Hessian set written to: {hessian_file}")

    if verbose:
        print(f"Found {hessian_set.n_results} hessian calculations")
        print(f"Found {hessian_set.n_molecules} hessian molecules")

    force_field = ForceField(initial_force_field, allow_cosmetic_attributes=True)

    records_and_molecules = list(hessian_set.to_records())

    if verbose:
        records_and_molecules = tqdm.tqdm(
            records_and_molecules,
            desc="Calculating parameters",
        )

    all_parameters = {
        "bond_eq": defaultdict(list),
        "bond_k": defaultdict(list),
        "angle_eq": defaultdict(list),
        "angle_k": defaultdict(list),
    }
    errored_records_and_molecules = []
    for record, molecule in records_and_molecules:
        try:
            parameters = calculate_parameters(record, molecule, force_field)
        except BaseException:
            errored_records_and_molecules.append((record, molecule))
            continue
        else:
            for key, values in parameters.items():
                for smirks, value in values.items():
                    all_parameters[key][smirks].extend(value)

    if working_directory is not None:
        seminario_file = Path(working_directory, "seminario-parameters.json")
        with open(seminario_file, "w") as json_file:
            json.dump(all_parameters, json_file, indent=2)

    if verbose:
        print(f"Found {len(errored_records_and_molecules)} errored calculations")

    if working_directory is not None:
        if len(errored_records_and_molecules):
            key = list(dataset.entries.keys())[0]
            opt_records_by_id = {
                record.record_id: record for record in hessian_set.entries[key]
            }
            records, _ = zip(*errored_records_and_molecules)
            errored_records = [opt_records_by_id[record.id] for record in records]
            errored_dataset = BasicResultCollection(entries={key: errored_records})
            error_file = Path(working_directory, "errored-dataset.json")
            with open(error_file, "w") as json_file:
                json_file.write(errored_dataset.json(indent=2))

            if verbose:
                print(f"Errored dataset written to: {error_file}")

    # Now we need to update the FF parameters
    kj_per_mol_per_nm_2 = unit.kilojoule_per_mole / unit.nanometer**2
    bond_handler = force_field.get_parameter_handler("Bonds")
    for smirks in all_parameters["bond_eq"]:
        bond = bond_handler.parameters[smirks]

        bond_length = numpy.mean(all_parameters["bond_eq"][smirks]) * unit.nanometer
        bond.length = bond_length.to(unit.angstrom)

        bond_k = numpy.mean(all_parameters["bond_k"][smirks]) * kj_per_mol_per_nm_2
        bond.k = bond_k.to(unit.kilocalorie_per_mole / unit.angstrom**2)

    kj_per_mol_per_rad_2 = unit.kilojoule_per_mole / unit.radian**2
    angle_handler = force_field.get_parameter_handler("Angles")
    for smirks in all_parameters["angle_eq"]:
        angle = angle_handler.parameters[smirks]

        angle_eq = numpy.mean(all_parameters["angle_eq"][smirks]) * unit.radian
        angle.angle = angle_eq.to(unit.degree)

        angle_k = numpy.mean(all_parameters["angle_k"][smirks]) * kj_per_mol_per_rad_2
        angle.k = angle_k.to(unit.kilocalorie_per_mole / unit.radian**2)

    force_field.to_file(output_force_field)


if __name__ == "__main__":
    main()
