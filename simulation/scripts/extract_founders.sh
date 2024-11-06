#!/bin/bash

abort () {
    aws stepfunctions send-task-failure --task-token "$TASK_TOKEN" \
        --error $? --cause "$BASH_COMMAND failed Extract Founders for $COHORT cohort and chr:$CHR"
    exit 1
}
trap 'abort' ERR

# send heartbeat to state machine
while true; do
    aws stepfunctions send-task-heartbeat --task-token "$TASK_TOKEN" || kill -- -$$ && sleep "$HEARTBEAT_INTERVAL"
done &

echo "Start Extract Founders for $COHORT cohort and chr:$CHR"

# INTERNAL FOLDER STRUCTURE
data_folder="/app/data"
results_folder="/app/results"
dataset_path="$data_folder/$COHORT"

# CREATE DATA FOLDER
mkdir -p $data_folder/fasta_hg38
mkdir -p $data_folder/hapmap_recombination_rate_hg38
mkdir -p $results_folder
mkdir -p $dataset_path

echo "Download data from s3"
aws s3 cp $S3_ORIGINAL/chr$CHR.vcf.gz $data_folder/source/ --quiet
aws s3 cp $S3_SIMULATED/$COHORT/phenotype-data/sample_map.tsv $data_folder/ --quiet
aws s3 cp $S3_SIMULATED/population_map.tsv $data_folder/ --quiet

echo "Function: extract gtype info"
samples=$(cut -f1 $data_folder/sample_map.tsv)

bcftools view -v snps -m2 -M2 -s $(echo $samples | tr ' ' '\n' | tail -n+2 | tr '\n' ',' | sed -e 's/,$//g') $data_folder/source/chr$CHR.vcf.gz -Ov -o $dataset_path/founders.vcf

cd $results_folder/
dataset_vcf="${COHORT}_set.vcf"
grep '#' $dataset_path/founders.vcf > header.txt
grep -v '#' $dataset_path/founders.vcf | awk '$8="."' | awk '$7="PASS"' | tr ' ' '\t' > body.txt
cat header.txt body.txt > $dataset_vcf
rm header.txt body.txt

# Recheck only bi-allelic variants (= remove multi-allelic variants)
# this step is probably not necessary!
# There were some issues with Alex's VCFs before and this step was needed
echo "Start filter_variants.py"
grep -v '#' $dataset_vcf | cut -f2 | sort | uniq -c | awk '$1>1' | awk '{print $2}' > remove_variants.txt
python /app/filter_variants.py --vcf-file $dataset_vcf --remove-variants-file remove_variants.txt --output-file aux && mv aux $dataset_vcf
echo "Done filter_variants.py"

# Chr format
sed -i "s/^chr//g" $dataset_vcf
sed -i "s/contig=<ID=chr$CHR>/contig=<ID=$CHR>/g" $dataset_vcf

echo "After filter_variants: $(cat $dataset_vcf | awk '{if ($1 !~/^#/) print}' | wc -l)"

awk 'BEGIN {FS=OFS="\t"} {
    for (i=10; i<=NF; i++) {  # Start from the 10th field where genotype information begins
        gsub("/","|",$i);
    }
    print;
}' $dataset_vcf > aux && mv aux $dataset_vcf

# Remove variants with unknown genotypes ('./.')
grep -v './.' $dataset_vcf > aux && mv aux $dataset_vcf

echo "After grep: $(cat $dataset_vcf | awk '{if ($1 !~/^#/) print}' | wc -l)"

########################################################
## GENERATE FOUNDER LABELS                            ##
########################################################
echo "Start generate_founder_labels.py"
python /app/generate_founder_labels.py \
    --vcf-file $results_folder/$dataset_vcf \
    --sample-map-file $data_folder/sample_map.tsv \
    --output-file $results_folder/Labels.generation0.txt
echo "Done generate_founder_labels.py"

echo "Finished Extract Founders for $COHORT and chr$CHR"
aws s3 cp $results_folder/$dataset_vcf $S3_SIMULATED/$COHORT/intermediate-data/simulation/chr$CHR/founders.vcf

bgzip $results_folder/$dataset_vcf
gzip $results_folder/Labels.generation0.txt
aws s3 cp $results_folder/$dataset_vcf.gz $S3_SIMULATED/$COHORT/genotype-data/master/chromosomes/chr$CHR/gen0/chr$CHR.vcf.gz
aws s3 cp $results_folder/Labels.generation0.txt.gz $S3_SIMULATED/$COHORT/phenotype-data/chromosomes/chr$CHR/gen0/ancestry.txt.gz

echo "Start conversion to x-array (zarr)"
python /app/to_xarray.py \
    --data-file $results_folder/$dataset_vcf.gz \
    --label-file $results_folder/Labels.generation0.txt.gz \
    --population-map $data_folder/population_map.tsv \
    --output-path $results_folder/chr$CHR.zarr

cd $results_folder && tar -czvf chr$CHR.zarr.tar.gz chr$CHR.zarr

echo "Finished conversion"
aws s3 cp $results_folder/chr$CHR.zarr.tar.gz $S3_SIMULATED/$COHORT/zarr-files/chromosomes/chr$CHR/gen0/
