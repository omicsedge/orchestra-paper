#!/bin/bash

echo "Start Generate Pedigree for $COHORT and default chr:$CHR"

# HARDCODED ARGUMENTS
generations_run=6
COHORT="${1,,}"
SAMPLE_MAP=$2
DATASET_CHR_PATH=$3
OUTPUT_DIR=$4

# INTERNAL FOLDER STRUCTURE

results_folder="$OUTPUT_DIR/results"
dataset_path="$OUTPUT_DIR/$COHORT"

# CREATE DATA FOLDER
mkdir -p $OUTPUT_DIR/fasta_hg38
mkdir -p $OUTPUT_DIR/hapmap_recombination_rate_hg38
mkdir -p $results_folder $dataset_path

echo "Download data from s3"
aws s3 cp $S3_BUCKET/reference-files/fasta/hg38/Homo_sapiens.GRCh38.dna.chromosome.$CHR.fa $OUTPUT_DIR/fasta_hg38 --quiet
aws s3 cp $S3_BUCKET/reference-files/fasta/hg38/Homo_sapiens.GRCh38.dna.chromosome.$CHR.fa.fai $OUTPUT_DIR/fasta_hg38 --quiet
aws s3 cp $S3_BUCKET/reference-files/recombination_map/hapmap/hg38/hapmap_recomb_hg38_chr$CHR.txt $OUTPUT_DIR/hapmap_recombination_rate_hg38/ --quiet


########################################################
## FIRST PART: RUN SIMULATIONS FOR Generating Pedigree #
########################################################

########################################################
## (STEP1) CHOOSE SAMPLES AND PARAMETERS FOR ANALYSIS ##
########################################################

# Pathways
# Change base_folder and slim_soft paths!
# Be careful! We need SLiM version 3.7!

# Script is run on only one dataset (test or train)

# Test and train samples (files with sample IDs)
samples=$(cut -f1 $SAMPLE_MAP)

echo "Number of samples: $samples"

##################################################
## (STEP2) CREATE VCF WITH ALL SELECTED SAMPLES ##
##################################################

# Obtain VCF with appropriate samples and chr:
# (i)   Extract test/train samples from source VCF
# (ii)  Take only biallelic variants and remove monomorphic
# (iii) Format changes: column INFO -number 8- to "." (needed for SLiM)
bcftools view -v snps -m2 -M2 -s $(echo $samples | tr ' ' '\n' | tail -n+2 | tr '\n' ',' | sed -e 's/,$//g') $DATASET_CHR_PATH -Ov -o $dataset_path/founders.vcf

cd $results_folder
dataset_vcf="${COHORT}_set.vcf"
echo "Create $results_folder/$dataset_vcf"

grep '#' $dataset_path/founders.vcf > header.txt
grep -v '#' $dataset_path/founders.vcf | awk '$8="."' | awk '$7="PASS"' | tr ' ' '\t' > body.txt
cat header.txt body.txt > $dataset_vcf
rm header.txt body.txt

echo "Optional step to remove multi-allelic variants"
# Recheck only bi-allelic variants (= remove multi-allelic variants)
# NOTE: this step is probably not necessary! There were some issues with
# Alex's VCFs before and this step was needed
grep -v '#' $dataset_vcf | cut -f2 | sort | uniq -c | awk '$1>1' | awk '{print $2}' > remove_variants.txt
python /app/filter_variants.py --vcf-file $dataset_vcf --remove-variants-file remove_variants.txt --output-file aux
mv aux $dataset_vcf

echo "Removed variants number: $(cat  remove_variants.txt | wc -l)"

echo "Remove chr prefix in chrom column if exist, update contig header."
sed -i "s/^chr//g" $dataset_vcf
sed -i "s/contig=<ID=chr$CHR>/contig=<ID=$CHR>/g" $dataset_vcf

# Remove variants with unknown genotypes ('./.')
grep -v './.' $dataset_vcf > aux && mv aux $dataset_vcf

######################
## (STEP3) RUN SLIM ##
######################

# Create file with (super)pop info per sample
individuals=$(grep '#C' ${COHORT}_set.vcf | tr '\t' '\n' | tail -n+10)
for ind in $individuals; do
    awk -v ind=$ind '$1==ind' $SAMPLE_MAP | cut -f1,3,2 | awk '{ print $1 "\t" $3 "\t" $2}'
done > aux.txt && mv aux.txt samples.keep


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
cd "/app/scripts" || exit
cp SLiM.initial_chr.$SIMULATION_TYPE.recipe SLiM.torun

sed -i "s,PARAMETER_CHR,$CHR,g" SLiM.torun
sed -i "s,PARAMETER_FASTA_HG38_PATHWAY,$data_folder/fasta_hg38,g" SLiM.torun
sed -i "s,PARAMETER_HAPMAP_PATHWAY,$data_folder/hapmap_recombination_rate_hg38,g" SLiM.torun
sed -i "s,PARAMETER_RESULTS_PATHWAY,$results_folder,g" SLiM.torun

if [ "$SIMULATION_TYPE" == "real" ]; then
    sed -i "s,PARAMETER_EUR_SAMPLESIZE,$(wc -l $results_folder/samples.EUR.keep | awk '{print $1}'),g" SLiM.torun
    sed -i "s,PARAMETER_EAS_SAMPLESIZE,$(wc -l $results_folder/samples.EAS.keep | awk '{print $1}'),g" SLiM.torun
    sed -i "s,PARAMETER_SAS_SAMPLESIZE,$(wc -l $results_folder/samples.SAS.keep | awk '{print $1}'),g" SLiM.torun
    sed -i "s,PARAMETER_AFR_SAMPLESIZE,$(wc -l $results_folder/samples.AFR.keep | awk '{print $1}'),g" SLiM.torun
    sed -i "s,PARAMETER_AMR_SAMPLESIZE,$(wc -l $results_folder/samples.AMR.keep | awk '{print $1}'),g" SLiM.torun
else
    sed -i "s,PARAMETER_SAMPLESIZE,$(wc -l $results_folder/samples.keep | awk '{print $1}'),g" SLiM.torun
fi

if [[ $COHORT == "test" ]]; then
    if [ "$SIMULATION_TYPE" == "real" ]; then
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_EUR,$(cat $results_folder/samplesize.tosimulate.test.EUR.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_EAS,$(cat $results_folder/samplesize.tosimulate.test.EAS.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_SAS,$(cat $results_folder/samplesize.tosimulate.test.SAS.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_AFR,$(cat $results_folder/samplesize.tosimulate.test.AFR.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_AMR,$(cat $results_folder/samplesize.tosimulate.test.AMR.txt),g" SLiM.torun
    else
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS,$(cat $results_folder/samplesize.tosimulate.test.txt),g" SLiM.torun
    fi

    sed -i "s,DIVISION_PARAMETER,4,g" SLiM.torun  
else
    if [ "$SIMULATION_TYPE" == "real" ]; then
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_EUR,$(cat $results_folder/samplesize.tosimulate.train.EUR.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_EAS,$(cat $results_folder/samplesize.tosimulate.train.EAS.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_SAS,$(cat $results_folder/samplesize.tosimulate.train.SAS.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_AFR,$(cat $results_folder/samplesize.tosimulate.train.AFR.txt),g" SLiM.torun
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS_AMR,$(cat $results_folder/samplesize.tosimulate.train.AMR.txt),g" SLiM.torun
    else
        sed -i "s,PARAMETER_NUMBER_SIMULATED_INDIVIDUALS,$(cat $results_folder/samplesize.tosimulate.train.txt),g" SLiM.torun
    fi

    sed -i "s,DIVISION_PARAMETER,1,g" SLiM.torun
fi

rm $results_folder/samplesize.tosimulate*


## Prepare SLiM recipe: one script per generation
for generation_round in $(seq 1 $generations_run); do 
    generations_output=$(($generation_round+1))
    cp SLiM.torun SLiM.torun_gen$generation_round
    sed -i "s,PARAMETER_GENERATIONS_N,$generations_output,g" SLiM.torun_gen$generation_round
done

echo '########################################'
echo "Run SLiM script for $COHORT cohort"
echo '########################################'

for gen in $(seq 1 $generations_run); do
    echo "Generation: $gen"
    slim SLiM.torun_gen$gen >/dev/null 2>&1
    mv $results_folder/output.trees $results_folder/output.generation_$gen.trees
    mv $results_folder/mating.txt $results_folder/mating.generation_$gen.txt
    mv $results_folder/death.txt $results_folder/death.generation_$gen.txt
done
rm SLiM.torun*

## Remove some files
echo "Copying results"
# Move matings and deaths recorded pedigrees to the tmp folder
mkdir -p tmp && mv $results_folder/death* $results_folder/mating*  $results_folder/pheno/pedigrees/

