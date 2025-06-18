#!/bin/bash
set -e
# Input arguments
generations_run=6
COHORT="${1,,}"
SAMPLE_MAP=$2
DATASET_CHR_PATH=$3
OUTPUT_DIR=$4
CHR=$5
SIMULATION_TYPE=$6
N_TIMES=$7

temp_dir=$(mktemp -d)
cd $temp_dir

echo "*** Start Generate Pedigree for $COHORT and default chr:$CHR ***"
echo "Working directory: $temp_dir"


########################################################
## FIRST PART: RUN SIMULATIONS FOR Generating Pedigree #
## (STEP1) CHOOSE SAMPLES AND PARAMETERS FOR ANALYSIS ##
########################################################

# Pathways
# Change base_folder and slim_soft paths!
# Be careful! We need SLiM version 3.7!

# Script is run on only one dataset (test or train)
# Test/train samples (files with sample IDs)
samples=$(cut -f1 $SAMPLE_MAP)

echo "Number of samples: $(wc -l < $SAMPLE_MAP)"

##################################################
## (STEP2) CREATE VCF WITH ALL SELECTED SAMPLES ##
##################################################

# Obtain VCF with appropriate samples and chr:
# (i)   Extract test/train samples from source VCF
# (ii)  Take only biallelic variants and remove monomorphic
# (iii) Format changes: column INFO -number 8- to "." (needed for SLiM)
bcftools view -v snps -m2 -M2 -s $(echo $samples | tr ' ' '\n' | tail -n+2 | tr '\n' ',' | sed -e 's/,$//g') $DATASET_CHR_PATH -Ov -o founders.vcf


echo "Create $(pwd)/${COHORT}_set.vcf"

grep '#' founders.vcf > header.txt
grep -v '#' founders.vcf | awk '$8="."' | awk '$7="PASS"' | tr ' ' '\t' > body.txt
cat header.txt body.txt > ${COHORT}_set.vcf

echo "Optional step to remove multi-allelic variants"
# Recheck only bi-allelic variants (= remove multi-allelic variants)
# NOTE: this step is probably not necessary! There were some issues with
grep -v '#' ${COHORT}_set.vcf | cut -f2 | sort | uniq -c | awk '$1>1' | awk '{print $2}' > remove_variants.txt
python /code/simulation/app/filter_variants.py --vcf-file ${COHORT}_set.vcf --remove-variants-file remove_variants.txt --output-file aux
mv aux ${COHORT}_set.vcf

echo "Removed variants number: $(cat  remove_variants.txt | wc -l)"

echo "Remove chr prefix in chrom column if exist, update contig header."
sed -i "s/^chr//g" ${COHORT}_set.vcf
sed -i "s/contig=<ID=chr$CHR>/contig=<ID=$CHR>/g" ${COHORT}_set.vcf

# Remove variants with unknown genotypes ('./.')
grep -v './.' ${COHORT}_set.vcf > aux && mv aux ${COHORT}_set.vcf

######################
## (STEP3) RUN SLIM ##
######################

# Create file with (super)pop info per sample
individuals=$(grep '#C' ${COHORT}_set.vcf | tr '\t' '\n' | tail -n+10)
for ind in $individuals; do
    awk -v ind=$ind '$1==ind' $SAMPLE_MAP | cut -f1,3,2 | awk '{ print $1 "\t" $3 "\t" $2}'
done > samples.keep


if [ $SIMULATION_TYPE == "real" ]; then
    # Input files:
    # (i)   Split VCF into different super populations 
    # (ii)  Remove snps with unknown genotypes - already done!
    # (iii) Remove useless gtype-associated metrics (dosages) for SLiM
    # (iv)  Obtain number of individuals that are going to be simulated
    sed -i "s/\(##contig=<ID=$CHR>\)/\1\n##INFO=<ID=ARTIFICIAL,Number=1,Type=String,Description='To track recombination'>/g" ${COHORT}_set.vcf
    for superpop in $(cut -f2 samples.keep | sort | uniq)
    do
        awk -v superpop=$superpop '$2==superpop' samples.keep | cut -f1 > samples.$superpop.keep
        bcftools view -s $(cat samples.$superpop.keep | tr '\n' ',' | sed -e 's/,$//g') ${COHORT}_set.vcf -Ov -o data.$superpop.vcf

        cat data.$superpop.vcf | sed -e 's,:[:,.0-9]*\t,\t,g' | sed -e 's,:[:,.0-9]*$,,g' | grep -v './.' > aux && mv aux data.$superpop.vcf 
        
        sample_size=$(cut -f2 samples.keep | grep $superpop | wc -l)
        echo $(($sample_size*2)) > samplesize.tosimulate.train.$superpop.txt
        echo $(($sample_size*8)) > samplesize.tosimulate.test.$superpop.txt
    done
else
    cat ${COHORT}_set.vcf | sed -e 's,:[:,.0-9]*\t,\t,g' | sed -e 's,:[:,.0-9]*$,,g' | grep -v './.' > data.vcf
    sample_size=$(cut -f2 samples.keep | wc -l)
    echo $(($sample_size*$N_TIMES)) > samplesize.tosimulate.train.txt
    echo $(($sample_size*4*$N_TIMES)) > samplesize.tosimulate.test.txt
fi


# Prepare SLiM recipe: use appropriate parameters
# (i)   chr + pathways (fasta hg38, hapmap, results folder)
# (ii)  sample sizes: initial and simulated
cp /code/simulation/scripts/SLiM.initial_chr.$SIMULATION_TYPE.recipe SLiM.torun

sed -i "s,PARAMETER_CHR,$CHR,g" SLiM.torun

if [ "$SIMULATION_TYPE" == "real" ]; then
    sed -i "s,PARAMETER_EUR_SAMPLESIZE,$(wc -l samples.EUR.keep | awk '{print $1}'),g" SLiM.torun
    sed -i "s,PARAMETER_EAS_SAMPLESIZE,$(wc -l samples.EAS.keep | awk '{print $1}'),g" SLiM.torun
    sed -i "s,PARAMETER_SAS_SAMPLESIZE,$(wc -l samples.SAS.keep | awk '{print $1}'),g" SLiM.torun
    sed -i "s,PARAMETER_AFR_SAMPLESIZE,$(wc -l samples.AFR.keep | awk '{print $1}'),g" SLiM.torun
    sed -i "s,PARAMETER_AMR_SAMPLESIZE,$(wc -l samples.AMR.keep | awk '{print $1}'),g" SLiM.torun
else
    sed -i "s,PARAMETER_SAMPLESIZE,$(wc -l samples.keep | awk '{print $1}'),g" SLiM.torun
fi

if [[ $COHORT == "test" ]]; then
    if [ "$SIMULATION_TYPE" == "real" ]; then
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_EUR,$(cat samplesize.tosimulate.test.EUR.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_EAS,$(cat samplesize.tosimulate.test.EAS.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_SAS,$(cat samplesize.tosimulate.test.SAS.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_AFR,$(cat samplesize.tosimulate.test.AFR.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_AMR,$(cat samplesize.tosimulate.test.AMR.txt),g" SLiM.torun
    else
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS,$(cat samplesize.tosimulate.test.txt),g" SLiM.torun
    fi

    sed -i "s,DIVISION_PARAMETER,4,g" SLiM.torun  
else
    if [ "$SIMULATION_TYPE" == "real" ]; then
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_EUR,$(cat samplesize.tosimulate.train.EUR.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_EAS,$(cat samplesize.tosimulate.train.EAS.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_SAS,$(cat samplesize.tosimulate.train.SAS.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_AFR,$(cat samplesize.tosimulate.train.AFR.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_AMR,$(cat samplesize.tosimulate.train.AMR.txt),g" SLiM.torun
    else
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS,$(cat samplesize.tosimulate.train.txt),g" SLiM.torun
    fi

    sed -i "s,DIVISION_PARAMETER,1,g" SLiM.torun
fi


## Prepare SLiM recipe: one script per generation
for generation_round in $(seq 1 $generations_run); do 
    generations_output=$(($generation_round+1))
    cp SLiM.torun SLiM.torun_gen$generation_round
    sed -i "s,PARAMETER_GENERATIONS_N,$generations_output,g" SLiM.torun_gen$generation_round
done

# Check if the required fasta file exists
FASTA_FILE="/data/fasta/$CHR.fa"
if [ ! -f "$FASTA_FILE" ]; then
    echo "ERROR: Required fasta file $FASTA_FILE does not exist!"
    echo "Please ensure the fasta file for chromosome $CHR is available at $FASTA_FILE"
    exit 1
fi

echo '########################################'
echo "Run SLiM script for $COHORT cohort"
echo '########################################'

for gen in $(seq 1 $generations_run); do
    echo "Generation: $gen"
    slim SLiM.torun_gen$gen >/dev/null 2>&1
    mv output.trees output.generation_$gen.trees
    mv mating.txt mating.generation_$gen.txt
    mv death.txt death.generation_$gen.txt
done

echo "Copying results"
OUTPUT_DIR=$OUTPUT_DIR/$COHORT/phenotype-data/pedigrees/
mkdir -p $OUTPUT_DIR && mv death* mating* $OUTPUT_DIR

rm -rf $temp_dir
