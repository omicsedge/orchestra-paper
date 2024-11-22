#!/usr/bin/env bash
set -ex

mkdir -p /results/simulation

cd /code/simulation
python simulate.py \
    -sc 1 -ec 22 \
    -sp /data/toy_example/Source_panel.vcf.gz \
    -sm /data/toy_example/SampleTable.forTraining.txt \
    -v "example-0.01" \
    -t "random" \
    -nt 2 \
    -o /results/simulation
