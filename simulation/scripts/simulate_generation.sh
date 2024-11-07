#!/bin/bash
set -e

# Input arguments
COHORT="${1,,}"
SAMPLE_MAP=$2
DATASET_CHR_PATH=$3
OUTPUT_DIR=$4
CHR=$5
GEN=$6
SIMULATION_TYPE=$7
N_TIMES=$8

temp_dir=$(mktemp -d)
cd $temp_dir

echo "*** Start Simulation for cohort:$COHORT, chr:$CHR and gen:$GEN ***"
echo "Working directory: $temp_dir"


# INTERNAL FOLDER STRUCTURE
cp $OUTPUT_DIR/$COHORT/genotype-data/master/chromosomes/chr$CHR/gen0/chr$CHR.vcf.gz dataset.vcf.gz
cp $OUTPUT_DIR/$COHORT/phenotype-data/pedigrees/death.generation_${GEN}.txt death.generation.txt
cp $OUTPUT_DIR/$COHORT/phenotype-data/pedigrees/mating.generation_${GEN}.txt mating.generation.txt
bgzip -d dataset.vcf.gz

echo "Function: SLiMulations following pedigree obtained in the first part"

# Create file with (super)pop info per sample
individuals=$(grep '#C' dataset.vcf | tr '\t' '\n' | tail -n+10)

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
' $SAMPLE_MAP > samples.keep

if [ "$SIMULATION_TYPE" = "real" ]; then
    # Input files:
    # (i)   Don't split VCF into different super populations - just super population files in the correct order! Don't need to specify migration with different superpops - already in mating files
    # (ii)  Remove snps with unknown genotypes - already done!
    # (iii) Remove useless gtype-associated metrics (dosages) for SLiM
    sed -i "s/\(##contig=<ID=$CHR>\)/\1\n##INFO=<ID=ARTIFICIAL,Number=1,Type=String,Description='To track recombination'>/g" dataset.vcf
    rm -f samples.ALL.keep    
    for superpop in "EUR" "EAS" "SAS" "AFR" "AMR"; do
        awk -v superpop=$superpop '$2==superpop' samples.keep | cut -f1 >> samples.ALL.keep
    done
    bcftools view -s $(cat samples.ALL.keep | tr '\n' ',' | sed -e 's/,$//g') dataset.vcf -Ov -o data.ALL.vcf
    cat data.ALL.vcf | sed -e 's,:[:,.0-9]*\t,\t,g' | sed -e 's,:[:,.0-9]*$,,g' | grep -v './.' > aux && mv aux data.ALL.vcf
    echo "$(wc -l samples.ALL.keep)"
else
    cat dataset.vcf | sed -e 's,:[:,.0-9]*\t,\t,g' | sed -e 's,:[:,.0-9]*$,,g' | grep -v './.' > data.ALL.vcf
    echo "$(wc -l samples.keep)"
fi


# Prepare SLiM recipe: use appropriate parameters
# (i)   chr + pathways (fasta hg38, hapmap, results folder)
# (ii)  use reproduce_pedigree recipe
cp /scripts/SLiM.reproduce_pedigree.recipe SLiM.torun

## Prepare SLiM recipe: one script per generation
sed -i "s,PARAMETER_CHR,$CHR,g" SLiM.torun
sed -i "s,PARAMETER_ALL_SAMPLESIZE,$(wc -l samples.keep | awk '{print $1}'),g" SLiM.torun
sed -i "s,PARAMETER_GENERATIONS_N,$(($GEN+1)),g" SLiM.torun

echo "Run slim script command. $(slim -v)"
slim SLiM.torun # >/dev/null 2>&1

# samples keep file
if [ "$SIMULATION_TYPE" = "real" ]; then
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
python3 /scripts/extract_tree_info.py \
        --input-file output.trees \
        --output-file trees.generation.vcf \
        --extract-individuals $(echo $(($keep_n*$N_TIMES)))


echo "Function: to obtain final Gtypes"

## (i) Update VCF with proper info (correct alternative gtypes and position)
## (ii) VCF format (add missing column) + individuals names
sed -i 's,|[2-9],|1,g' trees.generation.vcf
sed -i 's,[2-9]|,1|,g' trees.generation.vcf 

(grep '#' trees.generation.vcf && grep -v '#' trees.generation.vcf | awk '$2=$2+1' | tr ' ' '\t') > aux && mv aux trees.generation.vcf

sed -i "s,^1,$CHR,g" trees.generation.vcf
sed -i 's,PASS,1000\tPASS,g' trees.generation.vcf
sed -i "s,tsk_\([0-9]*\),i\1_${COHORT}_generation${GEN},g" trees.generation.vcf

## (iii) Use initial gtype file to recover those monomorphic positions that were removed by SLiM + store variant info to add to labels
# Recover monomorphic variants from initial vcf + store info from polymorphic variants for Labels
snps_initial_vcf=$(grep -v '#' dataset.vcf | cut -f2)
grep -v '##' trees.generation.vcf | cut -f1,2 > polymorphic_vts.generation.vcf
snps_slim_vcf=$(cut -f2 polymorphic_vts.generation.vcf | tail -n+2)
echo $snps_initial_vcf $snps_slim_vcf | tr ' ' '\n' | sort | uniq -u > monomorphic_snps.txt

if [ "$(cat monomorphic_snps.txt | wc -l)" != 0 ];
then
    # Add genotypes from monomorphic variants to simulated genomes + sort to obtain proper vcf
    rm -f add_monomorphic_vts.vcf
    awk 'FNR==NR{a[$1]=$1; next}; $2 in a {print $0;}' monomorphic_snps.txt dataset.vcf | cut -f11 > column.monomorphic_gtypes
    number_of_fake_gtypes=$(grep '#C' trees.generation.vcf | tr '\t' '\n' | tail -n+10 | wc -l)

    echo "rm(list=ls())" > command.R
    echo "column<-read.table('column.monomorphic_gtypes')" >> command.R
    echo "colnames(column)<-'repeat'" >> command.R
    echo "number<-$number_of_fake_gtypes" >> command.R
    echo "df<-column[,rep('repeat',number)]" >> command.R
    echo "write.table(df,file='add_monomorphic_vts.txt',sep='\t', quote=F, row.names = F, col.names = F)" >> command.R
    Rscript command.R
    rm command.R

    paste <(cat monomorphic_snps.txt | sed "s/^/$CHR\t/g" | awk '$3="."' | awk '$4=1000' | awk '$5="."' | awk '$6=1000' | awk '$7="PASS"' | awk '$8="."' | awk '$9="GT"' | tr ' ' '\t') <(cat add_monomorphic_vts.txt) > add_monomorphic_vts.vcf
    rm column.monomorphic_gtypes monomorphic_snps.txt add_monomorphic_vts.txt

    cat trees.generation.vcf add_monomorphic_vts.vcf > aux && mv aux trees.generation.vcf
    mv add_monomorphic_vts.vcf monomorphic_vts.generation.vcf
    vcf-sort trees.generation.vcf > aux && mv aux trees.generation.vcf
else
    mv monomorphic_snps.txt monomorphic_vts.generation.vcf
fi

############
## LABELS: #
############
# Create file with (super)pop info per sample
## (i) Join haplotypes in diploid genomes + add monomorphic variants that were removed by SLiM
echo "Start generate_labels function"
samples_n=$(cat samples.keep | wc -l)

python /scripts/generate_labels.py \
     --output-tree-file output.trees \
     --monomorphic-vts-file monomorphic_vts.generation.vcf \
     --polymorphic-vts-file polymorphic_vts.generation.vcf \
     --samples-order-file samples.order \
     --output-labels-file Labels.generation.txt.gz \
     --cohort $COHORT \
     --gen $GEN \
     --initial-genomes $(echo $(($samples_n*2))) \
     --total-genomes $(echo $(($samples_n*2*3))) \


echo "Finished Simulation for chr:$CHR and gen:$GEN"
mkdir -p $OUTPUT_DIR/$COHORT/genotype-data/master/chromosomes/chr$CHR/gen$GEN/
bgzip trees.generation.vcf && cp trees.generation.vcf.gz $OUTPUT_DIR/$COHORT/genotype-data/master/chromosomes/chr$CHR/gen$GEN/chr$CHR.vcf.gz
mkdir -p $OUTPUT_DIR/$COHORT/phenotype-data/chromosomes/chr$CHR/gen$GEN/
cp Labels.generation.txt.gz $OUTPUT_DIR/$COHORT/phenotype-data/chromosomes/chr$CHR/gen$GEN/ancestry.txt.gz

echo "Start conversion to x-array (zarr)"
python /scripts/to_xarray.py \
    --data-file trees.generation.vcf.gz \
    --label-file Labels.generation.txt.gz \
    --population-map $OUTPUT_DIR/population_map.tsv \
    --output-path chr$CHR.zarr

echo "Finished conversion"
tar -czf chr$CHR.zarr.tar.gz chr$CHR.zarr
mkdir -p $OUTPUT_DIR/$COHORT/zarr-files/chromosomes/chr$CHR/gen$GEN/
mv chr$CHR.zarr.tar.gz $OUTPUT_DIR/$COHORT/zarr-files/chromosomes/chr$CHR/gen$GEN/

rm -rf $temp_dir
