## Manifest

- `3-get-initial-force-field.py`. Add protein library charges and new SMIRNOFF types to handle delocalized charge parameters.
- `4-select-fit-parameters.py`. Get coverage of force field SMIRNOFF types in training dataset.
- `5-setup-forcebalance-inputs.py`. Setup training targets and input files for ForceBalance.
- `final-force-field.offxml`. Output of ForceBalance parameter fit. This is the force field for the null model.
- `forcebalance`. Directory containing output from the first of two ForceBalance runs.
- `forcebalance-2`. Directory containing output from the second of two ForceBalance runs.
- `initial-force-field.offxml`. Input for ForceBalance parameter fit. Output of `3-get-initial-force-field.py`.
- `training-smirks`. Directory containing parameter labels for training datasets and parameter SMIRKS for fit parameters as JSON files.

