#!/bin/bash

#python ../1-get-initial-force-field.py                                               \
#    --input-force-field          "openff_unconstrained-2.1.0.offxml"                 \
#    --library-charge-force-field "../library-charges/protein-library-charges.offxml" \
#    --output-force-field         "initial-force-field.offxml"

#python ../2-curate-training-datasets.py download-optimization                      \
#    --core-opt-dataset      "OpenFF Gen 2 Opt Set 1 Roche"                         \
#    --core-opt-dataset      "OpenFF Gen 2 Opt Set 2 Coverage"                      \
#    --core-opt-dataset      "OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy"            \
#    --core-opt-dataset      "OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy"        \
#    --core-opt-dataset      "OpenFF Gen 2 Opt Set 5 Bayer"                         \
#    --core-opt-dataset      "OpenFF Optimization Set 1"                            \
#    --core-opt-dataset      "SMIRNOFF Coverage Set 1"                              \
#    --core-opt-dataset      "OpenFF Aniline Para Opt v1.0"                         \
#    --core-opt-dataset      "OpenFF Protein Capped 1-mers 3-mers Optimization Dataset v1.0" \
#    --iodine-opt-dataset    "OpenFF Gen2 Optimization Dataset Protomers v1.0"      \
#    --iodine-opt-dataset    "OpenFF Iodine Chemistry Optimization Dataset v1.0"    \
#    --initial-force-field   "initial-force-field.offxml"                           \
#    --opt-records-to-remove "../optimization-records-to-remove.dat"                \
#    --max-opt-conformers    12                                                     \
#    --output-dataset-path   "training-datasets/optimization-training-dataset.json" \
#    --output-smirks-path    "training-datasets/optimization-training-smirks.json"  \
#    --verbose

#python ../2-curate-training-datasets.py download-torsiondrive                       \
#    --core-td-dataset        "OpenFF Gen 2 Torsion Set 1 Roche 2"                   \
#    --core-td-dataset        "OpenFF Gen 2 Torsion Set 2 Coverage 2"                \
#    --core-td-dataset        "OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2"      \
#    --core-td-dataset        "OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2"  \
#    --core-td-dataset        "OpenFF Gen 2 Torsion Set 5 Bayer 2"                   \
#    --core-td-dataset        "OpenFF Gen 2 Torsion Set 6 supplemental 2"            \
#    --protein-td-dataset     "OpenFF Protein Dipeptide 2-D TorsionDrive v2.1"       \
#    --protein-td-dataset     "OpenFF Protein Capped 1-mer Sidechains v1.2"          \
#    --aux-td-dataset         "SMIRNOFF Coverage Torsion Set 1"                      \
#    --aux-td-dataset         "OpenFF Group1 Torsions"                               \
#    --aux-td-dataset         "OpenFF Group1 Torsions 2"                             \
#    --aux-td-dataset         "OpenFF Group1 Torsions 3"                             \
#    --aux-td-dataset         "Pfizer discrepancy torsion dataset 1"                 \
#    --aux-td-dataset         "OpenFF Gen3 Torsion Set v1.0"                         \
#    --aux-td-dataset         "OpenFF Amide Torsion Set v1.0"                        \
#    --aux-td-dataset         "OpenFF WBO Conjugated Series v1.0"                    \
#    --aux-td-dataset         "OpenFF DANCE 1 eMolecules t142 v1.0"                  \
#    --initial-force-field    "initial-force-field.offxml"                           \
#    --explicit-ring-torsions "../explicit-ring-torsions.dat"                        \
#    --td-records-to-remove   "../torsiondrive-records-to-remove.dat"                \
#    --additional-td-records  "../additional-torsiondrive-records.json"              \
#    --cap-size               5                                                      \
#    --cap-method             "pick_heavy"                                           \
#    --n-processes            8                                                      \
#    --output-dataset-path    "training-datasets/torsiondrive-training-dataset.json" \
#    --output-smirks-path     "training-datasets/torsiondrive-training-smirks.json"  \
#    --verbose

#python ../3-get-msm-parameters.py                                                 \
#    --initial-force-field  "initial-force-field.offxml"                           \
#    --output-force-field   "msm-force-field.offxml"                               \
#    --optimization-dataset "training-datasets/optimization-training-dataset.json" \
#    --working-directory    "modified-seminario-method"                            \
#    --verbose

python ../4-create-forcebalance-inputs.py                                              \
    --tag                       "default-weights/tmp-forcebalance"                     \
    --force-field-path          "tmp-msm-force-field.offxml"                           \
    --optimization-dataset-path "training-datasets/optimization-training-dataset.json" \
    --torsiondrive-dataset-path "training-datasets/torsiondrive-training-dataset.json" \
    --valence-smirks-path       "training-datasets/optimization-training-smirks.json"  \
    --torsion-smirks-path       "training-datasets/torsiondrive-training-smirks.json"  \
    --smarts-to-exclude         "../smarts-to-exclude.dat"                             \
    --smiles-to-exclude         "../smiles-to-exclude.dat"                             \
    --protein-record-ids-path   "../protein-record-ids.dat"                            \
    --port                      55125                                                  \
    --verbose

