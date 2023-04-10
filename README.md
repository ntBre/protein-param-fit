# protein-param-fit

OpenFF force field parameter fits for proteins

## Manifest

- `1-curate-training-datasets.py`. Pulls down QC training datasets from QCArchive.
- `2-check-elf10-charges.py`. Filters molecules in the training dataset for which ELF10 charge assignment fails.
- `conda-envs`. Directory containing minimal and solved conda environments.
- `delocalized-charge-smirks`. Directory containing a force field with parameter SMIRKS for delocalized charges but no protein library charges.
- `library-charges`. Directory containing the derivation of library charges for protein residues.
- `null-model`. Directory containing the parameter fit for the null model, which uses Sage and protein QC training data but has no protein-specific SMIRNOFF types.
- `sage-behara.offxml`. Initial force field taken from Pavan Behara [here](https://github.com/MobleyLab/fitting-exp/blob/3aaf4c5d1ed757698c2bf4f910f4b074079c58e3/iter25/modified-force-field-trained-on-sage-targets.offxml)
- `test-delocalized-charge-assignments.py`. Identifies molecules in a QC dataset that are assigned different parameters to chemically equivalent atoms with delocalized charges.
- `training-datasets`. Directory containing QC training datasets as JSON files.

