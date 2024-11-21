############################################################
# IMPORTANT. PLEASE READ BEFORE EXECUTING THIS TOY EXAMPLE #
############################################################
#
####
# This step is merely a guide that outlines the process for obtaining and preparing publicly available genotype data from the 1000 Genomes Project, the Human Genome Diversity Project, 
# and the Simons Genome Diversity Project for local ancestry inference analysis using Orchestra. While this protocol guides you to assemble the same reference and source population panels, 
# we excluded the direct retrieval of SGDP data due to the size and download demands of these datasets. Instead, we have preloaded these panels for your use, and it is advisable to utilize 
# them as they include SGDP samples. Please proceed accordingly.
#
#
## Target admixed dataset:
# - Admixed Mexicans in Los Angeles, which was collected as part of the 1000 Genomes Project
#
## Reference panel comprising ancestral source populations:
# - Spanish (Iberians from 1000 Genomes Project)
# - Sub-Saharan Africans (Yoruba from 1000 Genomes Project)
# - Native Americans (Mayan and Pima populations from the Human Genome Diversity Project, and Mixe and Zapotec populations from the Simons Genome Diversity Project)
# 
# 
####
# Analyses: 
# (1) LAI analysis and LAD score:
# We have uploaded datasets that are ready for use with Orchestra for LAI analysis
# This dataset includes the target admixed population, the reference panel with ancestral source populations, and the SNP set featured in our article, which includes ancestry-associated and batch-corrected variants. 
# There is no need to retrieve or impute any data: you can retrieve those panels, merge source panels into a reference panel and proceed directly to Step2, immediately utilizing the provided panels for your LAI analysis with Orchestra.
# If you prefer to prepare the panels yourself, we have provided this Step 1 as a guide.
#
# (2) Fadm score and selection signal:
# If you wish to conduct a comprehensive analysis and also explore selection signals from adaptive admixture, calculating Fadm statistics in the process, you will need to access the complete genotype information. 
# You can obtain this dataset either through your own methods or by following this protocol completely.
# In this regard, this dataset should contain all SNPs (after applying a Minor Allele Frequency threshold to exclude low-frequency or rare variants) for frequency calculation. 
# We have also uploaded the SNP frequencies, including those observed in admixed Mexicans and those expected based on the ancestral source populations and global admixture proportions, to facilitate the calculation of Fadm statistics.
# Thus, you can skip the steps focused on obtaining this data and just proceed directly to the key sections.
#
#
####
# Recommended protocol:
# 1) No need to run Step1. Just run 'bcftools merge' (bcftools merge Source_NativeAmericans.reference_panel.vcf.gz Source_Europeans.reference_panel.vcf.gz Source_Africans.reference_panel.vcf.gz -Oz -o Source_panel.vcf.gz)
# to merge the source panels provided in GitHub/CodeOcean into a single reference panel that is ready for use with Orchestra. Redo folder structure pathways as this Step1 guide.
# - Target genotype data: Admixed_Mexicans.target_panel.vcf.gz
# - Sample ID-population membership map: SampleTable.forTrainings
# - Reference genotype data: Source_Africans.reference_panel.vcf.gz, Source_Europeans.reference_panel.vcf.gz and Source_NativeAmericans.reference_panel.vcf.gz (these should be merged: Source_panel.vcf.gz)
#
# 2) Run Orchestra with target (Admixed_Mexicans.target_panel.vcf.gz) and reference (Source_panel.vcf.gz) panels
# 3) Get LAD score (take inferred local ancestry results from Orchesra in output_inference/summary_results/smooth_samples.chr*.tsv.gz)
# 4) Get Fadm score. No need to run Step4 completely. Simply utilize the provided SNP frequency data located in the /reference_files/observed_and_expected_frequencies_for_Fadm folder. 
# Start from the 'Harmonize frequencies across same set of SNPs (source panels)' section and proceed to calculate the Fadm score using the expected and observed SNP frequencies provided.
#
# 5) Get final Manhattan plot and discover the selection signals
##########################################################################



##############################
###### Prepare SNP sets ######
##############################

### Software used ###
#
# bcftools v1.21 for VCF manipulation
#
# Picard tools for liftover (version 3.0.0; https://broadinstitute.github.io/picard/). It requires Java (we used version 17.0.10)
# UCSC chain files are needed to guide this operation, along with the reference sequence (in fasta format) for the target genome build (http://hgdownload.soe.ucsc.edu/downloads.html)
#
# Beagle for genotype imputation (v5.4; http://faculty.washington.edu/browning/beagle/beagle.html)
# And resources: HapMap genetic map, 1000 Genomes Project phase 3 reference panel 
#


### Pathways (change accordingly) ###
# Ensure you use the same analysis_folder consistently throughout all script steps in the toy example.
picard_path=""
folder_with_liftover_files=""

beagle_path=""
folder_with_imputation_resources=""

analysis_folder=""


### Get SNP sets into VCF format ###
## Download our ancestry-associated and batch-corrected SNP set from GitHub: QC_ancestry_associated_SNP_set.hg38.keep
## (this SNP set was for the custom-35pop reference panel used as input for Orchestra in our article Lerga-Jaso et al.)

cd $analysis_folder

echo '##fileformat=VCFv4.3
#CHROM POS ID REF ALT QUAL FILTER INFO' | tr ' ' '\t' > QC_ancestry_associated_SNP_set.hg38.vcf
cat QC_ancestry_associated_SNP_set.hg38.keep | awk '{print $1,$2,".",$3,$4,".",".","."}' | sed 's/^/chr/g' | tr ' ' '\t' >> QC_ancestry_associated_SNP_set.hg38.vcf
bgzip QC_ancestry_associated_SNP_set.hg38.vcf
tabix -p vcf QC_ancestry_associated_SNP_set.hg38.vcf.gz


### Liftover QC SNP set to hg19 ###
zcat QC_ancestry_associated_SNP_set.hg38.vcf.gz | sed 's/^chr//g' > QC_ancestry_associated_SNP_set.TO_LIFT.vcf

java -Xms20G -jar $picard_path/picard.jar LiftoverVcf \
    I=QC_ancestry_associated_SNP_set.TO_LIFT.vcf \
    O=lifted_over.QC_ancestry_associated_SNP_set.hg19.vcf \
    CHAIN=$folder_with_liftover_files/hg38ToHg19.over.chain \
    REJECT=rejected_variants.pseudo_VCF.vcf \
    R=$folder_with_liftover_files/hg19.fa RECOVER_SWAPPED_REF_ALT=true

rm QC_ancestry_associated_SNP_set.TO_LIFT.vcf rejected_variants.pseudo_VCF.vcf lifted_over.QC_ancestry_associated_SNP_set.hg19.vcf.idx

mv lifted_over.QC_ancestry_associated_SNP_set.hg19.vcf QC_ancestry_associated_SNP_set.hg19.vcf
bgzip QC_ancestry_associated_SNP_set.hg19.vcf
tabix -p vcf QC_ancestry_associated_SNP_set.hg19.vcf.gz



#################################
###### Prepare sample sets ######
#################################

### Get sample sets for analysis ###
# Source data (1KGP)
echo 'HG01500,HG01503,HG01504,HG01506,HG01507,HG01509,HG01510,HG01512,HG01513,HG01515,HG01516,HG01518,HG01519,HG01521,HG01522,HG01524,HG01525,HG01527,HG01528,HG01530,HG01531,HG01536,HG01602,HG01603,
HG01605,HG01606,HG01607,HG01608,HG01610,HG01612,HG01613,HG01615,HG01617,HG01618,HG01619,HG01620,HG01623,HG01624,HG01625,HG01631,HG01632,HG01670,HG01673,HG01675,HG01676,HG01678,HG01679,HG01680,HG01682,
HG01684,HG01685,HG01686,HG01694,HG01697,HG01699,HG01700,HG01702,HG01705,HG01707,HG01709,HG01746,HG01747,HG01756,HG01757,HG01761,HG01762,HG01765,HG01766,HG01767,HG01768,HG01770,HG01771,HG01773,HG01775,
HG01776,HG01777,HG01779,HG01781,HG01783,HG01784,HG02220,HG02221,HG02223,HG02224,HG02231,HG02232,HG02233,HG02235,HG02236,HG02238' | tr ',' '\n' > samples_1kgp.EUR.source.keep

echo 'NA18486,NA18487,NA18488,NA18489,NA18498,NA18499,NA18501,NA18502,NA18504,
NA18505,NA18507,NA18508,NA18510,NA18511,NA18516,NA18517,NA18519,NA18520,NA18522,NA18523,NA18852,NA18853,NA18855,NA18856,NA18858,NA18859,NA18861,NA18862,NA18864,NA18865,NA18867,NA18868,NA18870,NA18871,
NA18873,NA18874,NA18877,NA18878,NA18879,NA18881,NA18907,NA18908,NA18909,NA18910,NA18912,NA18913,NA18915,NA18916,NA18917,NA18923,NA18924,NA18933,NA18934,NA19092,NA19093,NA19095,NA19096,NA19098,NA19099,
NA19101,NA19102,NA19107,NA19108,NA19113,NA19114,NA19116,NA19117,NA19119,NA19121,NA19122,NA19127,NA19128,NA19130,NA19131,NA19137,NA19138,NA19140,NA19141,NA19143,NA19144,NA19146,NA19147,NA19149,NA19150,
NA19152,NA19153,NA19159,NA19160,NA19171,NA19172,NA19175,NA19176,NA19184,NA19185,NA19189,NA19190,NA19197,NA19198,NA19200,NA19201,NA19203,NA19204,NA19206,NA19207,NA19209,NA19210,NA19213,NA19214,NA19222,
NA19223,NA19225,NA19226,NA19235,NA19236,NA19239,NA19247,NA19248,NA19256,NA19257' | tr ',' '\n' > samples_1kgp.AFR.source.keep

cat samples_1kgp.EUR.source.keep samples_1kgp.AFR.source.keep > samples_1kgp.source.keep

# Source data (HGDP)
echo 'HGDP00854,HGDP00855,HGDP00856,HGDP00857,HGDP00858,HGDP00859,HGDP00860,HGDP00861,HGDP00862,HGDP00863,HGDP00864,HGDP00865,HGDP00868,HGDP00869,HGDP00870,HGDP00871,HGDP00872,HGDP00873,
HGDP00875,HGDP00876,HGDP00877,HGDP01037,HGDP01041,HGDP01043,HGDP01044,HGDP01050,HGDP01053,HGDP01055,HGDP01056,HGDP01057,HGDP01058,HGDP01059,HGDP01060' | tr ',' '\n' > samples_hgdp.NAM.source.keep 

# Source data (SGDP)
echo 'S_Mixe-2,B_Mixe-1,S_Mixe-3,S_Zapotec-1' | tr ',' '\n' > samples_sgdp.NAM.source.keep 

# Admixed population (1KGP)
echo 'NA19648,NA19649,NA19651,NA19652,NA19654,NA19655,NA19657,NA19658,NA19660,NA19661,NA19663,NA19664,NA19669,NA19670,NA19676,NA19678,NA19679,NA19681,NA19682,NA19684,NA19716,NA19717,NA19719,
NA19720,NA19722,NA19723,NA19725,NA19726,NA19728,NA19729,NA19731,NA19732,NA19734,NA19735,NA19740,NA19741,NA19746,NA19747,NA19749,NA19750,NA19752,NA19755,NA19756,NA19758,NA19759,NA19761,NA19762,
NA19764,NA19770,NA19771,NA19773,NA19774,NA19776,NA19777,NA19779,NA19780,NA19782,NA19783,NA19785,NA19786,NA19788,NA19789,NA19792,NA19794,NA19795' | tr ',' '\n' > samples_1kgp.admixed_Mexicans.keep



###################################################
###### Download and prepare 1KGP phase3 data ######
###################################################

### Download 1KGP data ###
download_1kgp() {
    chr=$1

    wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
    wget ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi


    ### Filter by SNP set (get SNPs > 0.005 MAF for Fadm statistic calculation) ###
    bcftools view -S samples_1kgp.admixed_Mexicans.keep -q 0.005:minor CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz -Ov -o Panel_MAF_filtered.chr${chr}.admixed_Mexicans.vcf.gz
    bcftools view -S samples_1kgp.source.keep -q 0.005:minor CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz -Ov -o Panel_MAF_filtered.chr${chr}.1kgp.vcf.gz

    ### Filter by SNP set (get ancestry-associate SNPs for Orchestra and LAD statistic calculation) ###
    bcftools isec -c none -p QC_filtered.chr${chr} -n=2 -w1 CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz QC_ancestry_associated_SNP_set.hg38.vcf.gz -Oz
    mv QC_filtered.chr${chr}/0000.vcf.gz ./chr${chr}.QC_SNP_set.vcf.gz
    mv QC_filtered.chr${chr}/0000.vcf.gz.tbi ./chr${chr}.QC_SNP_set.vcf.gz.tbi
    rm QC_filtered.chr${chr}/README.txt QC_filtered.chr${chr}/sites.txt
    rm -d QC_filtered.chr${chr}


    ### Extract samples ###
    bcftools view -S samples_1kgp.EUR.source.keep -Oz -o Source.chr${chr}.1kgp.EUR.vcf.gz chr${chr}.QC_SNP_set.vcf.gz
    bcftools view -S samples_1kgp.AFR.source.keep -Oz -o Source.chr${chr}.1kgp.AFR.vcf.gz chr${chr}.QC_SNP_set.vcf.gz
    bcftools view -S samples_1kgp.admixed_Mexicans.keep -Oz -o Admixed_Mexicans.chr${chr}.1kgp.vcf.gz chr${chr}.QC_SNP_set.vcf.gz

    rm chr${chr}.QC_SNP_set.vcf.gz chr${chr}.QC_SNP_set.vcf.gz.tbi CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
}

export -f download_1kgp
parallel --jobs 22 download_1kgp ::: {1..22}



############################################
###### Download and prepare HGDP data ######
############################################

### Download HGDP samples ####
for sample_HGDP in $(cat samples_hgdp.NAM.source.keep | tr ',' '\n'); do
    echo "#### ${sample_HGDP} ####"    
    wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/gVCFs/${sample_HGDP}.hgdp_wgs.20190516.vcf.gz
    wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/gVCFs/${sample_HGDP}.hgdp_wgs.20190516.vcf.gz.tbi
done


### Merge samples into one VCF ###
cp $folder_with_liftover_files/hg38.fa ./hg38.fasta
sed -i 's/^>/>chr/g' hg38.fasta
samtools faidx hg38.fasta
samtools dict hg38.fasta > hg38.dict

bcftools merge --gvcf hg38.fasta HGDP00854.hgdp_wgs.20190516.vcf.gz \
    HGDP00855.hgdp_wgs.20190516.vcf.gz \
    HGDP00856.hgdp_wgs.20190516.vcf.gz \
    HGDP00857.hgdp_wgs.20190516.vcf.gz \
    HGDP00858.hgdp_wgs.20190516.vcf.gz \
    HGDP00859.hgdp_wgs.20190516.vcf.gz \
    HGDP00860.hgdp_wgs.20190516.vcf.gz \
    HGDP00861.hgdp_wgs.20190516.vcf.gz \
    HGDP00862.hgdp_wgs.20190516.vcf.gz \
    HGDP00863.hgdp_wgs.20190516.vcf.gz \
    HGDP00864.hgdp_wgs.20190516.vcf.gz \
    HGDP00865.hgdp_wgs.20190516.vcf.gz \
    HGDP00868.hgdp_wgs.20190516.vcf.gz \
    HGDP00869.hgdp_wgs.20190516.vcf.gz \
    HGDP00870.hgdp_wgs.20190516.vcf.gz \
    HGDP00871.hgdp_wgs.20190516.vcf.gz \
    HGDP00872.hgdp_wgs.20190516.vcf.gz \
    HGDP00873.hgdp_wgs.20190516.vcf.gz \
    HGDP00875.hgdp_wgs.20190516.vcf.gz \
    HGDP00876.hgdp_wgs.20190516.vcf.gz \
    HGDP00877.hgdp_wgs.20190516.vcf.gz \
    HGDP01037.hgdp_wgs.20190516.vcf.gz \
    HGDP01041.hgdp_wgs.20190516.vcf.gz \
    HGDP01043.hgdp_wgs.20190516.vcf.gz \
    HGDP01044.hgdp_wgs.20190516.vcf.gz \
    HGDP01050.hgdp_wgs.20190516.vcf.gz \
    HGDP01053.hgdp_wgs.20190516.vcf.gz \
    HGDP01055.hgdp_wgs.20190516.vcf.gz \
    HGDP01056.hgdp_wgs.20190516.vcf.gz \
    HGDP01057.hgdp_wgs.20190516.vcf.gz \
    HGDP01058.hgdp_wgs.20190516.vcf.gz \
    HGDP01059.hgdp_wgs.20190516.vcf.gz \
    HGDP01060.hgdp_wgs.20190516.vcf.gz \
    -o HGDPcohort.g.vcf.gz -Oz

bcftools convert --gvcf2vcf HGDPcohort.g.vcf.gz --fasta-ref hg38.fasta -Oz -o HGDPcohort.vcf.gz 
tabix -p vcf HGDPcohort.vcf.gz

rm HGDP*.hgdp_wgs.*.vcf.gz HGDP*.hgdp_wgs.*.vcf.gz.tbi HGDPcohort.g.vcf.gz HGDPcohort.g.vcf.gz.tbi
rm hg38.dict hg38.fasta hg38.fasta.fai


### Imputation (Beagle) ###
# Split HGDP dataset into chromosomes
echo 'chr1 1,chr2 2,chr3 3,chr4 4,chr5 5,chr6 6,chr7 7,chr8 8,chr9 9,chr10 10,chr11 11,chr12 12,chr13 13,chr14 14,chr15 15,chr16 16,chr17 17,chr18 18,chr19 19,chr20 20,chr21 21,chr22 22' | tr ',' '\n' > rename_chrs.txt
cat rename_chrs.txt | awk '$3=$1' | awk '{print $2,$3}' > rename_chrs.reverse.txt
bcftools index -s HGDPcohort.vcf.gz | cut -f 1 | while read C; do bcftools view -O z -o HGDPcohort.split.${C}.vcf.gz HGDPcohort.vcf.gz "${C}" ; done
rm HGDPcohort.vcf.gz HGDPcohort.vcf.gz.tbi

for chr in $(seq 1 22); do

    echo '#######'
    echo "chr$chr"
    echo '#######'

    # Liftover (hg38 -> hg19)
    bcftools annotate --rename-chrs rename_chrs.txt HGDPcohort.split.chr${chr}.vcf.gz -Oz -o HGDPcohort.nochr_annot.chr${chr}.vcf.gz

    java -Xms20G -jar $picard_path/picard.jar LiftoverVcf \
        I=HGDPcohort.nochr_annot.chr${chr}.vcf.gz \
        O=HGDPcohort.liftover_hg19.chr${chr}.vcf.gz \
        CHAIN=$folder_with_liftover_files/hg38ToHg19.over.chain \
        REJECT=rejected_variants.pseudo_VCF.vcf \
        R=$folder_with_liftover_files/hg19.fa RECOVER_SWAPPED_REF_ALT=true

    rm rejected_variants.pseudo_VCF.vcf

    # Imputation (using 1KGP phase3 b37 as reference)
    java -Xmx40g -jar $beagle_path/beagle.22Jul22.46e.jar \
        chrom=${chr} \
        gt=HGDPcohort.liftover_hg19.chr${chr}.vcf.gz \
        ref=$folder_with_imputation_resources/chr${chr}.1kg.phase3.v5a.b37.bref3 \
        out=HGDPcohort.imputed.chr${chr} \
        map=$folder_with_imputation_resources/plink.GRCh37.map

    # Liftover (hg19 -> hg38)
    java -Xms20G -jar $picard_path/picard.jar LiftoverVcf \
        I=HGDPcohort.imputed.chr${chr}.vcf.gz \
        O=HGDPcohort.hg38.chr${chr}.vcf.gz \
        CHAIN=$folder_with_liftover_files/hg19ToHg38.over.chain \
        REJECT=rejected_variants.pseudo_VCF.vcf \
        R=$folder_with_liftover_files/hg38.fa RECOVER_SWAPPED_REF_ALT=true

    bcftools annotate --rename-chrs rename_chrs.reverse.txt HGDPcohort.hg38.chr${chr}.vcf.gz -Oz -o HGDPcohort.hg38.chr_annot.chr${chr}.vcf.gz
    mv HGDPcohort.hg38.chr_annot.chr${chr}.vcf.gz HGDPcohort.hg38.chr${chr}.vcf.gz
    tabix -p vcf -f HGDPcohort.hg38.chr${chr}.vcf.gz
    
    rm rejected_variants.pseudo_VCF.vcf HGDPcohort.imputed.chr${chr}.vcf.gz HGDPcohort.nochr_annot.chr${chr}.vcf.gz HGDPcohort.liftover_hg19.chr${chr}.vcf.gz HGDPcohort.liftover_hg19.chr${chr}.vcf.gz.tbi HGDPcohort.imputed.chr${chr}.log


    ### Filter by SNP set (get ancestry-associate SNPs for Orchestra and LAD statistic calculation) ###
    bcftools isec -c none -p QC_filtered.chr${chr} -n=2 -w1 HGDPcohort.hg38.chr${chr}.vcf.gz QC_ancestry_associated_SNP_set.hg38.vcf.gz -Oz
    mv QC_filtered.chr${chr}/0000.vcf.gz ./Source.chr${chr}.HGDP.NAM.vcf.gz
    mv QC_filtered.chr${chr}/0000.vcf.gz.tbi ./Source.chr${chr}.HGDP.NAM.vcf.gz.tbi
    rm QC_filtered.chr${chr}/README.txt QC_filtered.chr${chr}/sites.txt
    rm -d QC_filtered.chr${chr}

    ### Filter by SNP set (get SNPs for Fadm statistic calculation: WE USE SAME SET AS 1KGP SOURCE PANELS) ###
    tabix -p vcf -f Panel_MAF_filtered.chr${chr}.1kgp.vcf.gz
    bcftools isec -c none -p MAF_filtered.chr${chr} -n=2 -w1 HGDPcohort.hg38.chr${chr}.vcf.gz Panel_MAF_filtered.chr${chr}.1kgp.vcf.gz -Oz
    mv MAF_filtered.chr${chr}/0000.vcf.gz ./Panel_MAF_filtered.chr${chr}.HGDP.vcf.gz
    mv MAF_filtered.chr${chr}/0000.vcf.gz.tbi ./Panel_MAF_filtered.chr${chr}.HGDP.vcf.gz.tbi
    rm MAF_filtered.chr${chr}/README.txt MAF_filtered.chr${chr}/sites.txt
    rm -d MAF_filtered.chr${chr}

    rm HGDPcohort.hg38.chr${chr}.vcf.gz HGDPcohort.hg38.chr${chr}.vcf.gz.tbi HGDPcohort.split.chr${chr}.vcf.gz

done

rm rename_chrs.txt rename_chrs.reverse.txt


################################
###### Download SGDP data ######
################################
#
# We did not include the code for obtaining SGDP data to keep the setup simple and quick. 
# However, the uploaded genotype data includes additional SGDP individuals. Zapotec and Mixe individuals should be added to the reference panel using the same methods previously outlined:
# 
# - Mixe: mixe0002 (S_Mixe-2), mixe0007 (B_Mixe-1), mixe0042 (S_Mixe-3)
# - Zapotec: zapo0098 (S_Zapotec-1)
# 
# Access the source data here: https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/phased_data2021
#



##################################################
###### Finalize reference and target panels ######
##################################################

### Organize panels ###
cd $analysis_folder
mkdir -p Fadm_estimation
mkdir -p Target_and_reference_panels
mkdir -p Target_and_reference_panels/subsets_and_samples/

mv Panel_MAF_filtered.chr* Fadm_estimation/


### Concatenate chromosomes (target and reference panels) ###
bcftools concat Admixed_Mexicans.chr1.1kgp.vcf.gz Admixed_Mexicans.chr2.1kgp.vcf.gz Admixed_Mexicans.chr3.1kgp.vcf.gz Admixed_Mexicans.chr4.1kgp.vcf.gz Admixed_Mexicans.chr5.1kgp.vcf.gz Admixed_Mexicans.chr6.1kgp.vcf.gz Admixed_Mexicans.chr7.1kgp.vcf.gz Admixed_Mexicans.chr8.1kgp.vcf.gz Admixed_Mexicans.chr9.1kgp.vcf.gz Admixed_Mexicans.chr10.1kgp.vcf.gz Admixed_Mexicans.chr11.1kgp.vcf.gz Admixed_Mexicans.chr12.1kgp.vcf.gz Admixed_Mexicans.chr13.1kgp.vcf.gz Admixed_Mexicans.chr14.1kgp.vcf.gz Admixed_Mexicans.chr15.1kgp.vcf.gz Admixed_Mexicans.chr16.1kgp.vcf.gz Admixed_Mexicans.chr17.1kgp.vcf.gz Admixed_Mexicans.chr18.1kgp.vcf.gz Admixed_Mexicans.chr19.1kgp.vcf.gz Admixed_Mexicans.chr20.1kgp.vcf.gz Admixed_Mexicans.chr21.1kgp.vcf.gz Admixed_Mexicans.chr22.1kgp.vcf.gz -Oz -o Admixed_Mexicans.1kgp.vcf.gz

bcftools concat Source.chr1.1kgp.AFR.vcf.gz Source.chr2.1kgp.AFR.vcf.gz Source.chr3.1kgp.AFR.vcf.gz Source.chr4.1kgp.AFR.vcf.gz Source.chr5.1kgp.AFR.vcf.gz Source.chr6.1kgp.AFR.vcf.gz Source.chr7.1kgp.AFR.vcf.gz Source.chr8.1kgp.AFR.vcf.gz Source.chr9.1kgp.AFR.vcf.gz Source.chr10.1kgp.AFR.vcf.gz Source.chr11.1kgp.AFR.vcf.gz Source.chr12.1kgp.AFR.vcf.gz Source.chr13.1kgp.AFR.vcf.gz Source.chr14.1kgp.AFR.vcf.gz Source.chr15.1kgp.AFR.vcf.gz Source.chr16.1kgp.AFR.vcf.gz Source.chr17.1kgp.AFR.vcf.gz Source.chr18.1kgp.AFR.vcf.gz Source.chr19.1kgp.AFR.vcf.gz Source.chr20.1kgp.AFR.vcf.gz Source.chr21.1kgp.AFR.vcf.gz Source.chr22.1kgp.AFR.vcf.gz -Oz -o Source.AFR.1kgp.vcf.gz

bcftools concat Source.chr1.1kgp.EUR.vcf.gz Source.chr2.1kgp.EUR.vcf.gz Source.chr3.1kgp.EUR.vcf.gz Source.chr4.1kgp.EUR.vcf.gz Source.chr5.1kgp.EUR.vcf.gz Source.chr6.1kgp.EUR.vcf.gz Source.chr7.1kgp.EUR.vcf.gz Source.chr8.1kgp.EUR.vcf.gz Source.chr9.1kgp.EUR.vcf.gz Source.chr10.1kgp.EUR.vcf.gz Source.chr11.1kgp.EUR.vcf.gz Source.chr12.1kgp.EUR.vcf.gz Source.chr13.1kgp.EUR.vcf.gz Source.chr14.1kgp.EUR.vcf.gz Source.chr15.1kgp.EUR.vcf.gz Source.chr16.1kgp.EUR.vcf.gz Source.chr17.1kgp.EUR.vcf.gz Source.chr18.1kgp.EUR.vcf.gz Source.chr19.1kgp.EUR.vcf.gz Source.chr20.1kgp.EUR.vcf.gz Source.chr21.1kgp.EUR.vcf.gz Source.chr22.1kgp.EUR.vcf.gz -Oz -o Source.EUR.1kgp.vcf.gz

bcftools concat -a Source.chr1.HGDP.NAM.vcf.gz Source.chr2.HGDP.NAM.vcf.gz Source.chr3.HGDP.NAM.vcf.gz Source.chr4.HGDP.NAM.vcf.gz Source.chr5.HGDP.NAM.vcf.gz Source.chr6.HGDP.NAM.vcf.gz Source.chr7.HGDP.NAM.vcf.gz Source.chr8.HGDP.NAM.vcf.gz Source.chr9.HGDP.NAM.vcf.gz Source.chr10.HGDP.NAM.vcf.gz Source.chr11.HGDP.NAM.vcf.gz Source.chr12.HGDP.NAM.vcf.gz Source.chr13.HGDP.NAM.vcf.gz Source.chr14.HGDP.NAM.vcf.gz Source.chr15.HGDP.NAM.vcf.gz Source.chr16.HGDP.NAM.vcf.gz Source.chr17.HGDP.NAM.vcf.gz Source.chr18.HGDP.NAM.vcf.gz Source.chr19.HGDP.NAM.vcf.gz Source.chr20.HGDP.NAM.vcf.gz Source.chr21.HGDP.NAM.vcf.gz Source.chr22.HGDP.NAM.vcf.gz -Oz -o Source.NAM.hgdp.vcf.gz

rm Admixed_Mexicans.chr*.1kgp.vcf.gz
rm Source.chr*.1kgp.AFR.vcf.gz
rm Source.chr*.1kgp.EUR.vcf.gz
rm Source.chr*.HGDP.NAM.vcf.gz*


### Merge ancestral reference population into one source panel with common SNPs ###
bcftools index Source.EUR.1kgp.vcf.gz
bcftools index Source.AFR.1kgp.vcf.gz
bcftools index Source.NAM.hgdp.vcf.gz

bcftools isec -n=3 -c all -w1 Source.EUR.1kgp.vcf.gz Source.AFR.1kgp.vcf.gz Source.NAM.hgdp.vcf.gz -Oz -p CommonSNPs_SourcePanel

bcftools view -T CommonSNPs_SourcePanel/0000.vcf.gz Source.EUR.1kgp.vcf.gz -Oz -o EUR.common.vcf.gz
bcftools view -T CommonSNPs_SourcePanel/0000.vcf.gz Source.AFR.1kgp.vcf.gz -Oz -o AFR.common.vcf.gz
bcftools view -T CommonSNPs_SourcePanel/0000.vcf.gz Source.NAM.hgdp.vcf.gz -Oz -o NAM.common.vcf.gz

bcftools index EUR.common.vcf.gz
bcftools index AFR.common.vcf.gz
bcftools index NAM.common.vcf.gz

bcftools merge EUR.common.vcf.gz AFR.common.vcf.gz NAM.common.vcf.gz -Oz -o Source_panel.vcf.gz
rm EUR.common.vcf.gz AFR.common.vcf.gz NAM.common.vcf.gz EUR.common.vcf.gz.csi AFR.common.vcf.gz.csi NAM.common.vcf.gz.csi
rm CommonSNPs_SourcePanel/*
rm -d CommonSNPs_SourcePanel


### Filter target and reference panels by common SNPs ###
bcftools annotate -x ^FORMAT/GT Source_panel.vcf.gz -Oz -o Source_panel.cleaned.vcf.gz
mv Source_panel.cleaned.vcf.gz Source_panel.vcf.gz

tabix -p vcf Source_panel.vcf.gz
tabix -p vcf Admixed_Mexicans.1kgp.vcf.gz

bcftools isec -n=2 -c all -w1 Source_panel.vcf.gz Admixed_Mexicans.1kgp.vcf.gz -Oz -p CommonSNPs_bothPanels

bcftools view -T CommonSNPs_bothPanels/0000.vcf.gz Source_panel.vcf.gz -Oz -o Source.common.vcf.gz
bcftools view -T CommonSNPs_bothPanels/0000.vcf.gz Admixed_Mexicans.1kgp.vcf.gz -Oz -o Admixed_Mexicans.common.vcf.gz
mv Source.common.vcf.gz Source_panel.vcf.gz
mv Admixed_Mexicans.common.vcf.gz Admixed_Mexicans.1kgp.vcf.gz

rm *.tbi
rm CommonSNPs_bothPanels/*
rm -d CommonSNPs_bothPanels


### Filter duplicated entries ###
zgrep -v '#' Source_panel.vcf.gz | cut -f1,2 | sort | uniq -d > duplicated_vts.remove

bcftools view -T ^duplicated_vts.remove Admixed_Mexicans.1kgp.vcf.gz -o Admixed_Mexicans.f.vcf.gz
bcftools view -T ^duplicated_vts.remove Source_panel.vcf.gz -o Source_panel.f.vcf.gz
mv Admixed_Mexicans.f.vcf.gz Admixed_Mexicans.1kgp.vcf.gz
mv Source_panel.f.vcf.gz Source_panel.vcf.gz
rm duplicated_vts.remove


### Sample metadata for training (reference panel) ###
echo 'SampleID,level3,level1,Population,Dataset
HG01500,Iberian,EUR,Iberian,1kg30x
HG01503,Iberian,EUR,Iberian,1kg30x
HG01504,Iberian,EUR,Iberian,1kg30x
HG01506,Iberian,EUR,Iberian,1kg30x
HG01507,Iberian,EUR,Iberian,1kg30x
HG01509,Iberian,EUR,Iberian,1kg30x
HG01510,Iberian,EUR,Iberian,1kg30x
HG01512,Iberian,EUR,Iberian,1kg30x
HG01513,Iberian,EUR,Iberian,1kg30x
HG01515,Iberian,EUR,Iberian,1kg30x
HG01516,Iberian,EUR,Iberian,1kg30x
HG01518,Iberian,EUR,Iberian,1kg30x
HG01519,Iberian,EUR,Iberian,1kg30x
HG01521,Iberian,EUR,Iberian,1kg30x
HG01522,Iberian,EUR,Iberian,1kg30x
HG01524,Iberian,EUR,Iberian,1kg30x
HG01525,Iberian,EUR,Iberian,1kg30x
HG01527,Iberian,EUR,Iberian,1kg30x
HG01528,Iberian,EUR,Iberian,1kg30x
HG01530,Iberian,EUR,Iberian,1kg30x
HG01531,Iberian,EUR,Iberian,1kg30x
HG01536,Iberian,EUR,Iberian,1kg30x
HG01602,Iberian,EUR,Iberian,1kg30x
HG01603,Iberian,EUR,Iberian,1kg30x
HG01605,Iberian,EUR,Iberian,1kg30x
HG01606,Iberian,EUR,Iberian,1kg30x
HG01607,Iberian,EUR,Iberian,1kg30x
HG01608,Iberian,EUR,Iberian,1kg30x
HG01610,Iberian,EUR,Iberian,1kg30x
HG01612,Iberian,EUR,Iberian,1kg30x
HG01613,Iberian,EUR,Iberian,1kg30x
HG01615,Iberian,EUR,Iberian,1kg30x
HG01617,Iberian,EUR,Iberian,1kg30x
HG01618,Iberian,EUR,Iberian,1kg30x
HG01619,Iberian,EUR,Iberian,1kg30x
HG01620,Iberian,EUR,Iberian,1kg30x
HG01623,Iberian,EUR,Iberian,1kg30x
HG01624,Iberian,EUR,Iberian,1kg30x
HG01625,Iberian,EUR,Iberian,1kg30x
HG01631,Iberian,EUR,Iberian,1kg30x
HG01632,Iberian,EUR,Iberian,1kg30x
HG01670,Iberian,EUR,Iberian,1kg30x
HG01673,Iberian,EUR,Iberian,1kg30x
HG01675,Iberian,EUR,Iberian,1kg30x
HG01676,Iberian,EUR,Iberian,1kg30x
HG01678,Iberian,EUR,Iberian,1kg30x
HG01679,Iberian,EUR,Iberian,1kg30x
HG01680,Iberian,EUR,Iberian,1kg30x
HG01682,Iberian,EUR,Iberian,1kg30x
HG01684,Iberian,EUR,Iberian,1kg30x
HG01685,Iberian,EUR,Iberian,1kg30x
HG01686,Iberian,EUR,Iberian,1kg30x
HG01694,Iberian,EUR,Iberian,1kg30x
HG01697,Iberian,EUR,Iberian,1kg30x
HG01699,Iberian,EUR,Iberian,1kg30x
HG01700,Iberian,EUR,Iberian,1kg30x
HG01702,Iberian,EUR,Iberian,1kg30x
HG01705,Iberian,EUR,Iberian,1kg30x
HG01707,Iberian,EUR,Iberian,1kg30x
HG01709,Iberian,EUR,Iberian,1kg30x
HG01746,Iberian,EUR,Iberian,1kg30x
HG01747,Iberian,EUR,Iberian,1kg30x
HG01756,Iberian,EUR,Iberian,1kg30x
HG01757,Iberian,EUR,Iberian,1kg30x
HG01761,Iberian,EUR,Iberian,1kg30x
HG01762,Iberian,EUR,Iberian,1kg30x
HG01765,Iberian,EUR,Iberian,1kg30x
HG01766,Iberian,EUR,Iberian,1kg30x
HG01767,Iberian,EUR,Iberian,1kg30x
HG01768,Iberian,EUR,Iberian,1kg30x
HG01770,Iberian,EUR,Iberian,1kg30x
HG01771,Iberian,EUR,Iberian,1kg30x
HG01773,Iberian,EUR,Iberian,1kg30x
HG01775,Iberian,EUR,Iberian,1kg30x
HG01776,Iberian,EUR,Iberian,1kg30x
HG01777,Iberian,EUR,Iberian,1kg30x
HG01779,Iberian,EUR,Iberian,1kg30x
HG01781,Iberian,EUR,Iberian,1kg30x
HG01783,Iberian,EUR,Iberian,1kg30x
HG01784,Iberian,EUR,Iberian,1kg30x
HG02220,Iberian,EUR,Iberian,1kg30x
HG02221,Iberian,EUR,Iberian,1kg30x
HG02223,Iberian,EUR,Iberian,1kg30x
HG02224,Iberian,EUR,Iberian,1kg30x
HG02231,Iberian,EUR,Iberian,1kg30x
HG02232,Iberian,EUR,Iberian,1kg30x
HG02233,Iberian,EUR,Iberian,1kg30x
HG02235,Iberian,EUR,Iberian,1kg30x
HG02236,Iberian,EUR,Iberian,1kg30x
HG02238,Iberian,EUR,Iberian,1kg30x
NA18486,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18487,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18488,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18489,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18498,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18499,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18501,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18502,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18504,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18505,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18507,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18508,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18510,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18511,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18516,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18517,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18519,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18520,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18522,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18523,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18852,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18853,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18855,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18856,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18858,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18859,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18861,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18862,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18864,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18865,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18867,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18868,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18870,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18871,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18873,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18874,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18877,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18878,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18879,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18881,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18907,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18908,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18909,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18910,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18912,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18913,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18915,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18916,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18917,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18923,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18924,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18933,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA18934,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19092,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19093,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19095,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19096,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19098,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19099,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19101,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19102,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19107,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19108,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19113,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19114,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19116,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19117,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19119,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19121,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19122,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19127,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19128,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19130,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19131,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19137,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19138,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19140,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19141,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19143,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19144,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19146,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19147,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19149,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19150,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19152,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19153,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19159,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19160,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19171,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19172,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19175,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19176,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19184,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19185,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19189,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19190,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19197,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19198,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19200,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19201,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19203,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19204,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19206,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19207,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19209,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19210,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19213,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19214,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19222,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19223,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19225,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19226,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19235,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19236,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19239,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19247,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19248,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19256,SubSaharanAfrican,AFR,Yoruba,1kg30x
NA19257,SubSaharanAfrican,AFR,Yoruba,1kg30x
HGDP00854,NativeAmericans,NAM,Maya,HGDP
HGDP00855,NativeAmericans,NAM,Maya,HGDP
HGDP00856,NativeAmericans,NAM,Maya,HGDP
HGDP00857,NativeAmericans,NAM,Maya,HGDP
HGDP00858,NativeAmericans,NAM,Maya,HGDP
HGDP00859,NativeAmericans,NAM,Maya,HGDP
HGDP00860,NativeAmericans,NAM,Maya,HGDP
HGDP00861,NativeAmericans,NAM,Maya,HGDP
HGDP00862,NativeAmericans,NAM,Maya,HGDP
HGDP00863,NativeAmericans,NAM,Maya,HGDP
HGDP00864,NativeAmericans,NAM,Maya,HGDP
HGDP00865,NativeAmericans,NAM,Maya,HGDP
HGDP00868,NativeAmericans,NAM,Maya,HGDP
HGDP00869,NativeAmericans,NAM,Maya,HGDP
HGDP00870,NativeAmericans,NAM,Maya,HGDP
HGDP00871,NativeAmericans,NAM,Maya,HGDP
HGDP00872,NativeAmericans,NAM,Maya,HGDP
HGDP00873,NativeAmericans,NAM,Maya,HGDP
HGDP00875,NativeAmericans,NAM,Maya,HGDP
HGDP00876,NativeAmericans,NAM,Maya,HGDP
HGDP00877,NativeAmericans,NAM,Maya,HGDP
HGDP01037,NativeAmericans,NAM,Pima,HGDP
HGDP01041,NativeAmericans,NAM,Pima,HGDP
HGDP01043,NativeAmericans,NAM,Pima,HGDP
HGDP01044,NativeAmericans,NAM,Pima,HGDP
HGDP01050,NativeAmericans,NAM,Pima,HGDP
HGDP01053,NativeAmericans,NAM,Pima,HGDP
HGDP01055,NativeAmericans,NAM,Pima,HGDP
HGDP01056,NativeAmericans,NAM,Pima,HGDP
HGDP01057,NativeAmericans,NAM,Pima,HGDP
HGDP01058,NativeAmericans,NAM,Pima,HGDP
HGDP01059,NativeAmericans,NAM,Pima,HGDP
HGDP01060,NativeAmericans,NAM,Pima,HGDP' > SampleTable.forTraining.txt

mv Admixed_Mexicans.1kgp.vcf.gz Target_and_reference_panels/Admixed_Mexicans.vcf.gz
mv SampleTable.forTraining.txt Target_and_reference_panels/
mv Source_panel.vcf.gz Target_and_reference_panels/

mv Source.NAM.hgdp.vcf.gz* Target_and_reference_panels/subsets_and_samples/
mv Source.AFR.1kgp.vcf.gz* Target_and_reference_panels/subsets_and_samples/
mv Source.EUR.1kgp.vcf.gz* Target_and_reference_panels/subsets_and_samples/
mv samples_1kgp.admixed_Mexicans.keep Target_and_reference_panels/subsets_and_samples/
mv samples_1kgp.source.keep Target_and_reference_panels/subsets_and_samples/
mv samples_hgdp.NAM.source.keep Target_and_reference_panels/subsets_and_samples/
mv samples_1kgp.AFR.source.keep Target_and_reference_panels/subsets_and_samples/
mv samples_1kgp.EUR.source.keep Target_and_reference_panels/subsets_and_samples/

cd Target_and_reference_panels/


### Remove chr annotation from VCFs ###
echo 'chr1 1,chr2 2,chr3 3,chr4 4,chr5 5,chr6 6,chr7 7,chr8 8,chr9 9,chr10 10,chr11 11,chr12 12,chr13 13,chr14 14,chr15 15,chr16 16,chr17 17,chr18 18,chr19 19,chr20 20,chr21 21,chr22 22' | tr ',' '\n' > rename_chrs.txt
bcftools annotate --rename-chrs rename_chrs.txt Source_panel.vcf.gz -Oz -o Source_panel.nochr.vcf.gz
bcftools annotate --rename-chrs rename_chrs.txt Admixed_Mexicans.vcf.gz -Oz -o Admixed_Mexicans.nochr.vcf.gz
mv Source_panel.nochr.vcf.gz Source_panel.vcf.gz
mv Admixed_Mexicans.nochr.vcf.gz Admixed_Mexicans.vcf.gz
rm rename_chrs.txt


#### Description of files for LAI with Orchestra:
# Target panel:                Admixed_Mexicans.vcf.gz
# Reference panel:             Source_panel.vcf.gz
# Sample population map:       SampleTable.forTraining.txt


