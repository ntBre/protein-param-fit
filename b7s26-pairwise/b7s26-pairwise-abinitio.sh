#!/bin/bash

rsync -azv ../b7s26-nagl/initial-force-field.offxml .

rsync -azv ../b7s26-nagl/training-datasets .

rsync -azv ../b7s26-nagl/msm-force-field.offxml .

mkdir -p abinitio/forcebalance
rsync -azv ../b7s26-nagl/abinitio/forcebalance/forcefield abinitio/forcebalance/
rsync -azv ../b7s26-nagl/abinitio/forcebalance/targets.tar.gz abinitio/forcebalance/
rsync -azv ../b7s26-nagl/abinitio/forcebalance/optimize.in abinitio/forcebalance/

sed -i 's/AbInitio_SMIRNOFF/AbInitioPairwise_SMIRNOFF/g' abinitio/forcebalance/optimize.in

