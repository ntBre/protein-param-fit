#!/bin/bash

#python 1-get-initial-force-field.py                                                  \
#    --amber-nonbonded-path       "NBAmber.dat"                                       \
#    --input-force-field          "openff_unconstrained-2.1.0.offxml"                 \
#    --library-charge-force-field "../library-charges/protein-library-charges.offxml" \
#    --output-force-field         "initial-force-field.offxml"

#rsync -azv ../null-model/training-datasets .

#python ../3-get-msm-parameters.py                                                 \
#    --initial-force-field  "initial-force-field.offxml"                           \
#    --output-force-field   "msm-force-field.offxml"                               \
#    --optimization-dataset "training-datasets/optimization-training-dataset.json" \
#    --working-directory    "modified-seminario-method"                            \
#    --verbose

python ../4-create-forcebalance-inputs.py                                              \
    --tag                       "abinitio/forcebalance"                                \
    --force-field-path          "msm-force-field.offxml"                               \
    --optimization-dataset-path "training-datasets/optimization-training-dataset.json" \
    --torsiondrive-dataset-path "training-datasets/torsiondrive-training-dataset.json" \
    --valence-smirks-path       "training-datasets/optimization-training-smirks.json"  \
    --torsion-smirks-path       "training-datasets/torsiondrive-training-smirks.json"  \
    --smarts-to-exclude         "../smarts-to-exclude.dat"                             \
    --smiles-to-exclude         "../smiles-to-exclude.dat"                             \
    --protein-record-ids-path   "../protein-record-ids.dat"                            \
    --opt-geo-weight            0.005                                                  \
    --port                      55125                                                  \
    --torsiondrive-weight       0.01                                                   \
    --torsiondrive-target-type  "AbInitio"                                             \
    --verbose

