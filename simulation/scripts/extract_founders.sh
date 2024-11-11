#!/bin/bash
set -e

# Input arguments
COHORT="${1,,}"
SAMPLE_MAP=$2
DATASET_CHR_PATH=$3
OUTPUT_DIR=$4
CHR=$5

temp_dir=$(mktemp -d)
cd $temp_dir

echo "*** Start Extract Founders for $COHORT cohort and chr:$CHR ***"
echo "Working directory: $temp_dir"


samples=$(cut -f1 $SAMPLE_MAP)

echo "Function: extract gtype info"
bcftools view -v snps -m2 -M2 -s $(echo $samples | tr ' ' '\n' | tail -n+2 | tr '\n' ',' | sed -e 's/,$//g') $DATASET_CHR_PATH -Ov -o founders.vcf

grep '#' founders.vcf > dataset.vcf
grep -v '#' founders.vcf | awk '$8="."' | awk '$7="PASS"' | tr ' ' '\t' >> dataset.vcf


# Recheck only bi-allelic variants (= remove multi-allelic variants)
# this step is probably not necessary!
# There were some issues with Alex's VCFs before and this step was needed
echo "Start filter_variants.py"
grep -v '#' dataset.vcf | cut -f2 | sort | uniq -c | awk '$1>1' | awk '{print $2}' > remove_variants.txt
python /app/filter_variants.py --vcf-file dataset.vcf --remove-variants-file remove_variants.txt --output-file aux && mv aux dataset.vcf
echo "Done filter_variants.py"

# Chr format
sed -i "s/^chr//g" dataset.vcf
sed -i "s/contig=<ID=chr$CHR>/contig=<ID=$CHR>/g" dataset.vcf

echo "After filter_variants: $(cat dataset.vcf | awk '{if ($1 !~/^#/) print}' | wc -l)"

awk 'BEGIN {FS=OFS="\t"} {
    for (i=10; i<=NF; i++) {  # Start from the 10th field where genotype information begins
        gsub("/","|",$i);
    }
    print;
}' dataset.vcf > aux && mv aux dataset.vcf

# Remove variants with unknown genotypes ('./.')
grep -v './.' dataset.vcf > aux && mv aux dataset.vcf

echo "After grep: $(cat dataset.vcf | awk '{if ($1 !~/^#/) print}' | wc -l)"

########################################################
## GENERATE FOUNDER LABELS                            ##
########################################################
echo "Start generate_founder_labels.py"
python /app/generate_founder_labels.py \
    --vcf-file dataset.vcf \
    --sample-map-file $SAMPLE_MAP \
    --output-file ancestry.txt.gz # Labels.generation0.txt.gz
echo "Done generate_founder_labels.py"
mkdir -p $OUTPUT_DIR/$COHORT/phenotype-data/chromosomes/chr$CHR/gen0
cp ancestry.txt.gz $OUTPUT_DIR/$COHORT/phenotype-data/chromosomes/chr$CHR/gen0/

echo "Finished Extract Founders for $COHORT and chr$CHR"
mkdir -p $OUTPUT_DIR/$COHORT/genotype-data/chr$CHR/gen0/
bgzip dataset.vcf && cp dataset.vcf.gz $OUTPUT_DIR/$COHORT/genotype-data/chr$CHR/gen0/chr$CHR.vcf.gz

echo "Start conversion to x-array (zarr)"
python /app/to_xarray.py \
    --data-file dataset.vcf.gz \
    --label-file ancestry.txt.gz \
    --population-map $OUTPUT_DIR/population_map.tsv \
    --output-path chr$CHR.zarr

echo "Finished conversion"
mkdir -p $OUTPUT_DIR/$COHORT/zarr-files/chr$CHR/gen0
tar -czf chr$CHR.zarr.tar.gz chr$CHR.zarr
cp chr$CHR.zarr.tar.gz $OUTPUT_DIR/$COHORT/zarr-files/chr$CHR/gen0/

rm -rf $temp_dir