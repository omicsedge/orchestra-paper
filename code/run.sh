#!/usr/bin/env bash
set -ex

echo "This script will be running about 3h (chr 5 to 22)."

mkdir -p /results/simulation
mkdir -p /results/training
mkdir -p /results/inference
version="example-0.01"

cd /code/simulation
python simulation.py \
    -sc 5 -ec 22 \
    -sp /data/toy_example/Source_panel.vcf.gz \
    -sm /data/toy_example/SampleTable.forTraining.txt \
    -v $version \
    -t "random" \
    -nt 2 \
    -o /results/simulation

echo "✓ Completed simulation"

cd /code/training
for chr in "5-6" "7 8" "9 10" "11 12" "13 14" "15 16" "17 18" "19 22"; do
    start_chr=$(echo $chr | cut -d' ' -f1)
    end_chr=$(echo $chr | cut -d' ' -f2)

    python training.py \
        -sd /results/simulation \
        -sc $start_chr \
        -ec $end_chr \
        -ws 600 \
        -l 3 \
        -o /results/training \
        -v $version \
        -e 100

    echo "✓ Completed training chromosomes $start_chr-$end_chr"
done

find /results/training/ -type d -regex '.*/m[0-9]+\.[0-9]+' -exec rm -rf {} +

cd /code/inference
python inference.py \
    -p /data/toy_example/Admixed_Mexicans.target_panel.vcf.gz \
    -o /results/inference \
    -m /results/training/$version

 echo "✓ Completed inference"
