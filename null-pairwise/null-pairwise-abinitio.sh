#!/bin/bash

rsync -azv ../null-nagl/initial-force-field.offxml .

rsync -azv ../null-nagl/training-datasets .

rsync -azv ../null-nagl/msm-force-field.offxml .

mkdir -p abinitio/forcebalance
rsync -azv ../null-nagl/abinitio/forcebalance/forcefield abinitio/forcebalance/
rsync -azv ../null-nagl/abinitio/forcebalance/targets.tar.gz abinitio/forcebalance/
rsync -azv ../null-nagl/abinitio/forcebalance/optimize.in abinitio/forcebalance/

sed -i 's/AbInitio_SMIRNOFF/AbInitioPairwise_SMIRNOFF/g' abinitio/forcebalance/optimize.in

