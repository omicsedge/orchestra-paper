#!/bin/bash
set -e

abort () {
    aws stepfunctions send-task-failure --task-token "$TASK_TOKEN" \
        --error $? --cause "$BASH_COMMAND failed Simulation for chr:$CHR and gen:$GEN"
    exit 1
}
trap 'abort' ERR

# send heartbeat to state machine
while true; do
    aws stepfunctions send-task-heartbeat --task-token "$TASK_TOKEN" || kill -- -$$ && sleep "$HEARTBEAT_INTERVAL"
done &

echo "Start Simulation for chr:$CHR and gen:$GEN"

# FETCH CLI ARGUMENTS
chr=$CHR
generation=$GEN
type=$SIMULATION_TYPE
dataset=$COHORT

# INTERNAL FOLDER STRUCTURE
data_folder="/app/data"
results_folder="/app/results"
scripts_folder="/app/scripts"

echo "Download data from s3"
aws s3 cp $S3_SIMULATED/$COHORT/phenotype-data/sample_map.tsv $data_folder/ --quiet
aws s3 cp $S3_SIMULATED/$COHORT/intermediate-data/simulation/chr${chr}/founders.vcf $results_folder/chr$chr/$dataset/${dataset}_set.vcf --quiet
aws s3 cp $S3_BUCKET/reference-files/fasta/hg38/Homo_sapiens.GRCh38.dna.chromosome.$chr.fa $data_folder/fasta_hg38/ --quiet
aws s3 cp $S3_BUCKET/reference-files/fasta/hg38/Homo_sapiens.GRCh38.dna.chromosome.$chr.fa.fai $data_folder/fasta_hg38/ --quiet
aws s3 cp $S3_BUCKET/reference-files/recombination_map/hapmap/hg38/hapmap_recomb_hg38_chr$chr.txt $data_folder/hapmap_recombination_rate_hg38/ --quiet
aws s3 cp $S3_SIMULATED/$COHORT/phenotype-data/pedigrees/death.generation_${generation}.txt $results_folder/ --quiet
aws s3 cp $S3_SIMULATED/$COHORT/phenotype-data/pedigrees/mating.generation_${generation}.txt $results_folder/ --quiet
aws s3 cp $S3_SIMULATED/population_map.tsv $data_folder/ --quiet

echo "Start Simulation for chr:$CHR and gen:$GEN"

echo "Function: SLiMulations following pedigree obtained in the first part"

# Create file with (super)pop info per sample
cd $results_folder/chr$chr/$dataset/
individuals=$(grep '#C' ${COHORT}_set.vcf | tr '\t' '\n' | tail -n+10)

awk -v inds="$individuals" '
BEGIN {
    FS=OFS="\t"
    split(inds, ind_array, " ")
}

FNR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == "sample_id") idx1 = i
        if ($i == "superpopulation_code") idx2 = i
        if ($i == "level_3_code") idx3 = i
    }
}

FNR > 1 {
    for (i in ind_array) {
        if (ind_array[i] == $idx1) {
            print $idx1, $idx2, $idx3
        }
    }
}
' $data_folder/sample_map.tsv > samples.keep

if [ $type == 'real' ]; then
    # Input files:
    # (i)   Don't split VCF into different super populations - just super population files in the correct order! Don't need to specify migration with different superpops - already in mating files
    # (ii)  Remove snps with unknown genotypes - already done!
    # (iii) Remove useless gtype-associated metrics (dosages) for SLiM
    sed -i "s/\(##contig=<ID=$chr>\)/\1\n##INFO=<ID=ARTIFICIAL,Number=1,Type=String,Description='To track recombination'>/g" ${COHORT}_set.vcf
    rm -f samples.ALL.keep    
    for superpop in $(echo "EUR EAS SAS AFR AMR" | tr ' ' '\n')
    do
        awk -v superpop=$superpop '$2==superpop' samples.keep | cut -f1 >> samples.ALL.keep
    done
    bcftools view -s $(cat samples.ALL.keep | tr '\n' ',' | sed -e 's/,$//g') ${COHORT}_set.vcf -Ov -o data.ALL.vcf
    cat data.ALL.vcf | sed -e 's,:[:,.0-9]*\t,\t,g' | sed -e 's,:[:,.0-9]*$,,g' | grep -v './.' > aux && mv aux data.ALL.vcf
    echo "$(wc -l samples.ALL.keep)"
else
    cat ${COHORT}_set.vcf | sed -e 's,:[:,.0-9]*\t,\t,g' | sed -e 's,:[:,.0-9]*$,,g' | grep -v './.' > data.ALL.vcf
    echo "$(wc -l samples.keep)"
fi


# Prepare SLiM recipe: use appropriate parameters
# (i)   chr + pathways (fasta hg38, hapmap, results folder)
# (ii)  use reproduce_pedigree recipe
cd $scripts_folder
cp SLiM.reproduce_pedigree.recipe SLiM.torun

## Prepare SLiM recipe: one script per generation
generations_output=$(($generation+1))
sed -i "s,PARAMETER_CHR,$chr,g" SLiM.torun
sed -i "s,PARAMETER_FASTA_HG38_PATHWAY,$data_folder/fasta_hg38,g" SLiM.torun
sed -i "s,PARAMETER_HAPMAP_PATHWAY,$data_folder/hapmap_recombination_rate_hg38,g" SLiM.torun
sed -i "s,PARAMETER_RESULTS_PATHWAY,$results_folder/chr$chr/$dataset,g" SLiM.torun
sed -i "s,PARAMETER_MATING_DEATH_PATHWAY,$results_folder,g" SLiM.torun
sed -i "s,PARAMETER_ALL_SAMPLESIZE,$(wc -l $results_folder/chr$chr/$dataset/samples.keep | awk '{print $1}'),g" SLiM.torun
sed -i "s,PARAMETER_GENERATIONS_MATING_DEATH,$generation,g" SLiM.torun
sed -i "s,PARAMETER_GENERATIONS_N,$generations_output,g" SLiM.torun


echo "Run slim script command. $(slim -v)"
slim SLiM.torun

cd $results_folder/chr$chr/$dataset/

# samples keep file
if [ $type == 'real' ]; then
    keep_file=samples.ALL.keep
else
    keep_file=samples.keep
fi

# Build sample order file
awk 'OFS="\t"{ print NR*2-2,$1,$3; print NR*2-1,$1,$3 }' $keep_file > samples.order


## Run extraction - extract info from tree sequences (Labels + Gtypes)
## Be careful! High memory usage
echo "Function: to extract info from tree sequence: $COHORT"
keep_n=$(cat $keep_file | wc -l)
python3 $scripts_folder/extract_tree_info.py \
        --results-path $results_folder/chr$chr/$dataset \
        --gen $generation \
        --extract-individuals $(echo $(($keep_n*$N_TIMES)))


echo "Function: to obtain final Gtypes"

## (i) Update VCF with proper info (correct alternative gtypes and position)
## (ii) VCF format (add missing column) + individuals names
sed -i 's,|[2-9],|1,g' trees.generation_${generation}.vcf
sed -i 's,[2-9]|,1|,g' trees.generation_${generation}.vcf 

(grep '#' trees.generation_${generation}.vcf && grep -v '#' trees.generation_${generation}.vcf | awk '$2=$2+1' | tr ' ' '\t') > aux && mv aux trees.generation_${generation}.vcf

sed -i "s,^1,$chr,g" trees.generation_${generation}.vcf
sed -i 's,PASS,1000\tPASS,g' trees.generation_${generation}.vcf
sed -i "s,tsk_\([0-9]*\),i\1_${dataset}_generation${generation},g" trees.generation_${generation}.vcf

## (iii) Use initial gtype file to recover those monomorphic positions that were removed by SLiM + store variant info to add to labels
# Recover monomorphic variants from initial vcf + store info from polymorphic variants for Labels
snps_initial_vcf=$(grep -v '#' ${dataset}_set.vcf | cut -f2)
grep -v '##' trees.generation_${generation}.vcf | cut -f1,2 > polymorphic_vts.generation_${generation}.vcf
snps_slim_vcf=$(cut -f2 polymorphic_vts.generation_${generation}.vcf | tail -n+2)
echo $snps_initial_vcf $snps_slim_vcf | tr ' ' '\n' | sort | uniq -u > monomorphic_snps.txt

if [ "$(cat monomorphic_snps.txt | wc -l)" != 0 ];
then
    # Add genotypes from monomorphic variants to simulated genomes + sort to obtain proper vcf
    rm -f add_monomorphic_vts.vcf
    awk 'FNR==NR{a[$1]=$1; next}; $2 in a {print $0;}' monomorphic_snps.txt ${dataset}_set.vcf | cut -f11 > column.monomorphic_gtypes
    number_of_fake_gtypes=$(grep '#C' trees.generation_${generation}.vcf | tr '\t' '\n' | tail -n+10 | wc -l)

    echo "rm(list=ls())" > command.R
    echo "column<-read.table('$results_folder/chr$chr/$dataset/column.monomorphic_gtypes')" >> command.R
    echo "colnames(column)<-'repeat'" >> command.R
    echo "number<-$number_of_fake_gtypes" >> command.R
    echo "df<-column[,rep('repeat',number)]" >> command.R
    echo "write.table(df,file='$results_folder/chr$chr/$dataset/add_monomorphic_vts.txt',sep='\t', quote=F, row.names = F, col.names = F)" >> command.R
    Rscript command.R
    rm command.R

    paste <(cat monomorphic_snps.txt | sed "s/^/$chr\t/g" | awk '$3="."' | awk '$4=1000' | awk '$5="."' | awk '$6=1000' | awk '$7="PASS"' | awk '$8="."' | awk '$9="GT"' | tr ' ' '\t') <(cat add_monomorphic_vts.txt) > add_monomorphic_vts.vcf
    rm column.monomorphic_gtypes monomorphic_snps.txt add_monomorphic_vts.txt

    cat trees.generation_${generation}.vcf add_monomorphic_vts.vcf > aux && mv aux trees.generation_${generation}.vcf
    mv add_monomorphic_vts.vcf monomorphic_vts.generation_${generation}.vcf
    vcf-sort trees.generation_${generation}.vcf > aux && mv aux trees.generation_${generation}.vcf
else
    mv monomorphic_snps.txt monomorphic_vts.generation_${generation}.vcf
fi

############
## LABELS: #
############
# Create file with (super)pop info per sample
## (i) Join haplotypes in diploid genomes + add monomorphic variants that were removed by SLiM
echo "Start generate_labels function"
samples_n=$(cat samples.keep | wc -l)

python $scripts_folder/generate_labels.py \
      --results-path $results_folder/chr$chr/$dataset \
      --cohort $COHORT \
      --gen $generation \
      --initial-genomes $(echo $(($samples_n*2))) \
      --total-genomes $(echo $(($samples_n*2*3))) \
      --samples-order $results_folder/chr$chr/$dataset/samples.order


bgzip $results_folder/chr$chr/$dataset/trees.generation_${generation}.vcf
gzip $results_folder/chr$chr/$dataset/Labels.generation${generation}.txt

echo "Finished Simulation for chr:$CHR and gen:$GEN"
aws s3 cp $results_folder/chr$chr/$dataset/trees.generation_${GEN}.vcf.gz $S3_SIMULATED/$COHORT/genotype-data/master/chromosomes/chr$CHR/gen$GEN/chr$chr.vcf.gz
aws s3 cp $results_folder/chr$chr/$dataset/Labels.generation${GEN}.txt.gz $S3_SIMULATED/$COHORT/phenotype-data/chromosomes/chr$CHR/gen$GEN/ancestry.txt.gz

echo "Start conversion to x-array (zarr)"
python /app/to_xarray.py \
    --data-file $results_folder/chr$chr/$dataset/trees.generation_${generation}.vcf.gz \
    --label-file $results_folder/chr$chr/$dataset/Labels.generation${generation}.txt.gz \
    --population-map $data_folder/population_map.tsv \
    --output-path $results_folder/chr$chr/$dataset/chr$CHR.zarr

cd $results_folder/chr$chr/$dataset/ && tar -czvf chr$CHR.zarr.tar.gz chr$CHR.zarr

echo "Finished conversion"
aws s3 cp $results_folder/chr$chr/$dataset/chr$CHR.zarr.tar.gz $S3_SIMULATED/$COHORT/zarr-files/chromosomes/chr$CHR/gen$GEN/
