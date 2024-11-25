#!/usr/bin/env bash
set -ex

echo "This script will be running about 3h (chr 5 to 22)."

EXPERIMENT_NAME="example-0.01"

python /main.py simulation \
    -sc 5 -ec 22 \
    -sp /data/toy_example/Source_panel.vcf.gz \
    -sm /data/toy_example/SampleTable.forTraining.txt \
    -v $EXPERIMENT_NAME \
    -t "random" \
    -nt 2 \
    -o /results/simulation

echo "✓ Completed simulation"

for chr in "5 6" "7 8" "9 10" "11 12" "13 14" "15 16" "17 18" "19 22"; do
    start_chr=$(echo $chr | cut -d' ' -f1)
    end_chr=$(echo $chr | cut -d' ' -f2)

    python /main.py training \
        -sd /results/simulation \
        -sc $start_chr \
        -ec $end_chr \
        -ws 600 \
        -l 3 \
        -o /results/training \
        -v $EXPERIMENT_NAME \
        -e 100

    rm -rf /results/training/m$start_chr.$end_chr

    echo "✓ Completed training chromosomes $start_chr-$end_chr"
done

python /main.py inference \
    -p /data/toy_example/Admixed_Mexicans.target_panel.vcf.gz \
    -o /results/inference \
    -m /results/training/$EXPERIMENT_NAME

 echo "✓ Completed inference"
