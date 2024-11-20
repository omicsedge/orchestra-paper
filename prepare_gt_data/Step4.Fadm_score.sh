
#######################
## PRELIMINARY STEPS ##
#######################

### Software ###
#
# PLINK v2.00 for SNP frequency calculation
#


### Pathways (change accordingly) ###
plink2_path=""
analysis_folder=""
Fadm_statistic_calc_folder="$analysis_folder/Fadm_estimation"
samples_folder="$analysis_folder/Target_and_reference_panels/subsets_and_samples"


### Genotype data stored (MAF > 0.005 filtered SNP set) ###
# Check files:
# - Panel_MAF_filtered.chr${chr}.admixed_Mexicans.vcf.gz ==> Admixed Mexicans set
# - Panel_MAF_filtered.chr${chr}.1kgp.vcf.gz ==> AFR and EUR sets
# - Panel_MAF_filtered.chr${chr}.HGDP.vcf.gz ==> NAM set


### Separate AFR and EUR sets from 1KGP ###
cd $Fadm_statistic_calc_folder

run_commands() {
    chr=$1
    samples_path=$2

    bcftools view -S $samples_path/samples_1kgp.EUR.source.keep Panel_MAF_filtered.chr${chr}.1kgp.vcf.gz -Oz -o Panel_MAF_filtered.chr${chr}.EUR.vcf.gz
    bcftools view -S $samples_path/samples_1kgp.AFR.source.keep Panel_MAF_filtered.chr${chr}.1kgp.vcf.gz -Oz -o Panel_MAF_filtered.chr${chr}.AFR.vcf.gz
}
export -f run_commands

parallel --jobs 4 run_commands ::: {1..22} ::: $samples_folder


### Re-label VCFs ###
for chr in $(seq 22 -1 1); do
    mv Panel_MAF_filtered.chr${chr}.HGDP.vcf.gz Panel_MAF_filtered.chr${chr}.NAM.vcf.gz
done


### Rename SNP IDs ###
for chr in $(seq 22 -1 1); do
    bcftools annotate -x ID Panel_MAF_filtered.chr${chr}.NAM.vcf.gz | bcftools annotate --set-id +'%CHROM:%POS:%REF:%ALT' -Oz -o Panel_MAF_filtered.chr${chr}.NAM.annotated.vcf.gz
done



################################################################################################################################################
####################################### OBTAIN OBSERVED SNP FREQUENCIES IN THE ADMIXED TARGET POPULATION #######################################
################################################################################################################################################

### Calculate observed frequencies ###
cd $Fadm_statistic_calc_folder
mkdir -p Observed_freq && cd Observed_freq

run_commands() {
    chr=$1
    plink_path=$2

    $plink_path/plink2 --vcf ../Panel_MAF_filtered.chr${chr}.admixed_Mexicans.vcf.gz --freq --out chr${chr}.admixed_Mexicans.freq --threads 5
    rm chr${chr}.admixed_Mexicans.freq.log
}
export -f run_commands

parallel --jobs 4 run_commands ::: {1..22} ::: $plink2_path



##################################################################################################################################################
####################################### OBTAIN EXPECTED SNP FREQUENCIES USING ANCESTRAL SOURCE POPULATIONS #######################################
##################################################################################################################################################

### Average estimated adxmiture proportions ###
# Check file: Admixture_proportions.txt


############################################################
## Observed SNP frequency per ancestral source population ##
############################################################

### Observed SNP frequency: Iberians from 1KGP ###
cd $Fadm_statistic_calc_folder

mkdir -p Expected_freq/Source_pop.SNP_freq
cd Expected_freq/Source_pop.SNP_freq

get_freq() {
    chr=$1
    plink_path=$2

    $plink_path/plink2 --vcf ../../Panel_MAF_filtered.chr${chr}.EUR.vcf.gz --freq --out chr${chr}.EUR.freq --threads 50 --vcf-half-call m
}
export -f get_freq
parallel --jobs 22 get_freq ::: {1..22} ::: $plink2_path


### Observed SNP frequency: Yoruba from 1KGP ###
get_freq() {
    chr=$1
    plink_path=$2

    $plink_path/plink2 --vcf ../../Panel_MAF_filtered.chr${chr}.AFR.vcf.gz --freq --out chr${chr}.AFR.freq --threads 50 --vcf-half-call m
}
export -f get_freq
parallel --jobs 22 get_freq ::: {1..22} ::: $plink2_path


### Observed SNP frequency: Native Americans from 1KGP ###
get_freq() {
    chr=$1
    plink_path=$2

    $plink_path/plink2 --vcf ../../Panel_MAF_filtered.chr${chr}.NAM.annotated.vcf.gz --freq --out chr${chr}.NAM.freq --threads 50 --vcf-half-call m
    sed -i 's/chr//g' chr${chr}.NAM.freq.afreq
}
export -f get_freq
parallel --jobs 22 get_freq ::: {1..22} ::: $plink2_path

rm *.log



###################################################################
## Harmonize frequencies across same set of SNPs (source panels) ##
###################################################################

### Take common SNPs and same format ###
cd $Fadm_statistic_calc_folder/Expected_freq/Source_pop.SNP_freq

for chr in $(seq 22 -1 1); do
    echo "# chr$chr #"

    # Join by common positions
    grep -v "^#" chr${chr}.AFR.freq.afreq | cut -f2-6 > AFR_data
    grep -v "^#" chr${chr}.EUR.freq.afreq | cut -f2-6 > EUR_data
    grep -v "^#" chr${chr}.NAM.freq.afreq | cut -f2-6 > NAM_data

    sort -k1,1 AFR_data > AFR_sorted
    sort -k1,1 EUR_data > EUR_sorted
    sort -k1,1 NAM_data > NAM_sorted
    join -1 1 -2 1 AFR_sorted EUR_sorted > temp1
    join -1 1 -2 1 temp1 NAM_sorted > joined_data
    rm temp1

    # Adjust frequency using same allele as reference
    # If alleles do not match, remove that row 
    awk '{
        ID = $1
        AFR_REF = $2; AFR_ALT = $3; AFR_FREQ = $4; AFR_OBS = $5
        EUR_REF = $6; EUR_ALT = $7; EUR_FREQ = $8; EUR_OBS = $9
        NAM_REF = $10; NAM_ALT = $11; NAM_FREQ = $12; NAM_OBS = $13

        # Check allele matching
        if (AFR_REF == EUR_REF && AFR_ALT == EUR_ALT) {
            # Same orientation
            EUR_FREQ_adj = EUR_FREQ
        } else if (AFR_REF == EUR_ALT && AFR_ALT == EUR_REF) {
            # Alleles are swapped
            EUR_FREQ_adj = 1 - EUR_FREQ
        } else {
            # Alleles do not match
            next
        }

        if (AFR_REF == NAM_REF && AFR_ALT == NAM_ALT) {
            # Same orientation
            NAM_FREQ_adj = NAM_FREQ
        } else if (AFR_REF == NAM_ALT && AFR_ALT == NAM_REF) {
            # Alleles are swapped
            NAM_FREQ_adj = 1 - NAM_FREQ
        } else {
            # Alleles do not match
            next
        }

        # Output the SNP with adjusted frequencies
        print ID, AFR_REF, AFR_ALT, AFR_FREQ, AFR_OBS, EUR_FREQ_adj, EUR_OBS, NAM_FREQ_adj, NAM_OBS
    }' joined_data > common_snps.txt

    # Extract final frequencies using common SNPs
    awk '{print $1, $2, $3, $4, $5}' common_snps.txt | awk -v chr=$chr '$6=chr' | awk '{print $6, $1,$2,$3,$4,$5}' | tr ' ' '\t' > chr${chr}.AFR.afreq
    awk '{print $1, $2, $3, $6, $7}' common_snps.txt | awk -v chr=$chr '$6=chr' | awk '{print $6, $1,$2,$3,$4,$5}' | tr ' ' '\t' > chr${chr}.EUR.afreq
    awk '{print $1, $2, $3, $8, $9}' common_snps.txt | awk -v chr=$chr '$6=chr' | awk '{print $6, $1,$2,$3,$4,$5}' | tr ' ' '\t' > chr${chr}.NAM.afreq

    rm common_snps.txt joined_data

done
rm chr*freq.afreq
rm AFR_sorted NAM_sorted EUR_sorted
rm AFR_data NAM_data EUR_data



##############################
## Get expected frequencies ##
##############################

### Exclude multi-allelic variants ###
cd $Fadm_statistic_calc_folder/Expected_freq/Source_pop.SNP_freq
grep '[AGCT],[AGCT]' * | cut -f2  | sort | uniq > exclude.SNPs

exclude=$(cat exclude.SNPs | tr '\n' ',' | sed 's/,$/\n/g' | sed 's;,;\\|;g')
if [ "$exclude" != "" ]; then
    for chr in $(seq 22 -1 1); do
        grep -v "$exclude" chr${chr}.AFR.afreq > aux && mv aux chr${chr}.AFR.afreq
        grep -v "$exclude" chr${chr}.EUR.afreq > aux && mv aux chr${chr}.EUR.afreq
        grep -v "$exclude" chr${chr}.NAM.afreq > aux && mv aux chr${chr}.NAM.afreq
    done
fi
rm exclude.SNPs


### Get expected frequency using SNP frequency from ancestral source populations ###
cd $Fadm_statistic_calc_folder/Expected_freq

echo "
for(chr in c(1:22)){
    print('###############')
    print(paste0('Chr: ',chr))
    print('###############')

    ### Load admixture proportions from Step3 (LAD score) ###
    pop.freqs <- read.table('../../Admixture_proportions.txt', sep='\t', header=T) 


    ### Create a list of data frames for each population file ###
    pop.files <- unique(paste0('Source_pop.SNP_freq/chr',chr,'.',pop.freqs\$Ancestry,'.afreq'))
    pop.list_data <- lapply(pop.files, function(x) {
      tmp <- (read.table(x, header=F, sep='\t', comment.char='#'))
      colnames(tmp) <- c('chr','snp','ref','alt','freq','count')
      tmp\$pop <- gsub('.*\\\.','',gsub('.afreq','',x))
      return(tmp)
    })


    ### Combine population data into one data frame ###
    pop.data <- do.call(rbind, pop.list_data)


    ### Merge with population frequencies ###
    pop.data.with_freq_allchrs <- merge(pop.data, pop.freqs, by.y=c('Ancestry'), by.x=c('pop'))

     
    ### Calculate expected frequency per SNP and allele ###
    pop.data.with_freq_allchrs\$exp_freq <- pop.data.with_freq_allchrs\$freq * pop.data.with_freq_allchrs\$Percentage


    ### Summarize by SNP and allele ###
    result.with_freq_allchrs <- aggregate(exp_freq ~ snp + ref + alt, data=pop.data.with_freq_allchrs, FUN=sum)


    ### Obtain max and min freq in source populations to remove SNPs with frequencies outside of the range ###
    library(dplyr)
    principal_pops = pop.freqs\$Ancestry

    aux.allchrs <- pop.data[which(pop.data\$pop %in% principal_pops),]

    min_freq.allchrs <- aggregate(freq ~ snp, data = aux.allchrs, FUN = min)
    max_freq.allchrs <- aggregate(freq ~ snp, data = aux.allchrs, FUN = max)
    min_max_freq.allchrs <- cbind(min_freq.allchrs, max_freq.allchrs)
    min_max_freq.allchrs <- min_max_freq.allchrs[, c(1,2,4)]
    colnames(min_max_freq.allchrs) <- c('SNP', 'Min', 'Max')
    rownames(min_max_freq.allchrs) <- min_max_freq.allchrs\$SNP


    ### Rename columns ###
    colnames(result.with_freq_allchrs) <- c('ID', 'REF', 'ALT', 'exp.freq')
    rownames(result.with_freq_allchrs) <- as.character(result.with_freq_allchrs\$ID)
    result.with_freq_allchrs <- result.with_freq_allchrs[pop.list_data[[1]]\$snp,]
    result.with_freq_allchrs\$chr <- chr

    if(chr == 1){
        final_df.with_freq_allchrs <- result.with_freq_allchrs    
        final_mm.with_freq_allchrs <- min_max_freq.allchrs[rownames(result.with_freq_allchrs),]
    }else{
        final_df.with_freq_allchrs <- rbind(final_df.with_freq_allchrs, result.with_freq_allchrs)
        final_mm.with_freq_allchrs <- rbind(final_mm.with_freq_allchrs,
                                            min_max_freq.allchrs[rownames(result.with_freq_allchrs),])
    }
}
write.table(final_df.with_freq_allchrs, file='allchrs.expected_freq.txt', col.names=T, row.names=F, sep='\t', quote=F)
write.table(final_mm.with_freq_allchrs, file='minmax_freqs.allchrs.txt', col.names=T, row.names=F, sep='\t', quote=F)
" > Rcommand.R
Rscript Rcommand.R
rm Rcommand.R *.log


# Description:
## - allchrs.expected_freq.txt file
# This file includes the expected frequency for each SNP, derived from the reference populations that serve as ancestral sources. 
# SNP frequencies from source populations are weighted according to the global (genome-wide) admixture proportions determined by Orchestra for the admixed population.
#
## - minmax_freqs.allchrs.txt file
# This file presents the minimum and maximum frequencies observed in the source populations for each SNP.
# These values represent the potential range within which the frequencies might vary in the admixed population.
#



################################################################################################################################################
####################################### OBTAIN FADM STATISTIC USING OBSERVED VS EXPECTED SNP FREQUENCIES #######################################
################################################################################################################################################

#############################################
## Prepare observed and expected frequency ##
#############################################

### Merge all chromosome (observed frequency) ###
cd $Fadm_statistic_calc_folder/
mkdir -p Fadm 

cd Observed_freq
cat chr1.admixed_Mexicans.freq.afreq chr2.admixed_Mexicans.freq.afreq chr3.admixed_Mexicans.freq.afreq chr4.admixed_Mexicans.freq.afreq chr5.admixed_Mexicans.freq.afreq chr6.admixed_Mexicans.freq.afreq chr7.admixed_Mexicans.freq.afreq chr8.admixed_Mexicans.freq.afreq chr9.admixed_Mexicans.freq.afreq chr10.admixed_Mexicans.freq.afreq chr11.admixed_Mexicans.freq.afreq chr12.admixed_Mexicans.freq.afreq chr13.admixed_Mexicans.freq.afreq chr14.admixed_Mexicans.freq.afreq chr15.admixed_Mexicans.freq.afreq chr16.admixed_Mexicans.freq.afreq chr17.admixed_Mexicans.freq.afreq chr18.admixed_Mexicans.freq.afreq chr19.admixed_Mexicans.freq.afreq chr20.admixed_Mexicans.freq.afreq chr21.admixed_Mexicans.freq.afreq chr22.admixed_Mexicans.freq.afreq > ../Fadm/observed_freq.txt

cd $Fadm_statistic_calc_folder/Fadm/
grep -v '#CHROM' observed_freq.txt > aux && mv aux observed_freq.txt


### Annotate SNP position (hg38) ###
cd $Fadm_statistic_calc_folder/Fadm
awk 'BEGIN{FS=OFS="\t"} NR==1{print $1, "POS", $2, $3, $4, $5, $6} NR>1{split($2, a, ":"); print $1, a[2], $2, $3, $4, $5, $6}' observed_freq.txt > aux && mv aux observed_freq.txt


### Get expected frequency ###
cp $Fadm_statistic_calc_folder/Expected_freq/allchrs.expected_freq.txt $Fadm_statistic_calc_folder/Fadm/expected_freq.txt


## Use the same SNP set (observed and expected from both reference and target panels)
sort -k1,1 expected_freq.txt > sorted_expected.txt
sort -k3,3 observed_freq.txt > sorted_observed.txt

join -1 1 -2 3 -o 1.1 1.2 1.3 1.4 1.5 2.1 2.2 2.3 2.4 2.5 2.6 2.7 -t $'\t' sorted_expected.txt sorted_observed.txt > final_joined_output.txt

cat final_joined_output.txt | cut -f1-5 | tr ' ' '\t' > expected_freq.common_snps.txt
cat final_joined_output.txt | cut -f6-12 | tr ' ' '\t' > observed_freq.common_snps.txt

rm final_joined_output.txt sorted_expected.txt sorted_observed.txt



##############################
## Calculate Fadm statistic ##
##############################

### R script to get Fadm statistics ###
cd $Fadm_statistic_calc_folder/Fadm
echo "
    rm(list=ls())

    ### Load observed and expected frequencies ###
    print('Loading data')
    expected <- read.table('./expected_freq.common_snps.txt', header=F)
    colnames(expected) <- c('ID', 'REF', 'ALT', 'FREQ', 'CHR')

    observed <- read.table('./observed_freq.common_snps.txt', header = F, comment.char='')
    colnames(observed) <- c('CHR', 'POS', 'ID', 'REF', 'ALT', 'FREQ', 'COUNTS')


    ### Remove duplicated variants ###
    duplicates <- c(observed[duplicated(observed\$ID),]\$ID, 
                    expected[duplicated(expected\$ID),]\$ID)
    observed <- observed[!(observed\$ID %in% duplicates),]
    expected <- expected[!(expected\$ID %in% duplicates), ]


    ### Get SNPs with same ref-alt orientation ###
    if( all(observed\$ID == expected\$ID) ){
        print('First check - fine')
    }else{
        print('Different SNP IDs')
    }

    if( all(observed\$REF == expected\$REF) ){
        print('Second check - fine')
    }else{
        filter <- observed\$REF == expected\$REF
        observed <- observed[filter, ]
        expected <- expected[filter, ]
        expected_chrpacks <- expected_chrpacks[filter, ]
    }

    observed\$FREQ <- as.numeric(observed\$FREQ)
    expected\$FREQ <- as.numeric(expected\$FREQ)


    ### Get Fadm statistic ###
    numerador <- ((observed\$FREQ - expected\$FREQ) ** 2) + (((1 - observed\$FREQ) - (1 - expected\$FREQ)) ** 2)
    denominador <- 2 * (1 - ((expected\$FREQ ** 2) + ((1 - expected\$FREQ) ** 2)))
    Fadm <- numerador / denominador
    observed\$FADM <- Fadm
    colnames(observed)[6] <- 'FREQ_OBSERVED'
    colnames(expected)[4] <- 'FREQ_EXPECTED'
    expected <- cbind(observed, expected[4])
    final_df <- expected

    print('Ready to write...')
    write.table(final_df, file='Fadm_score.txt', col.names=T,row.names=F,sep='\t', quote=F)
" > Rcommand.R
Rscript Rcommand.R
rm Rcommand.R



########################################
## Additional quality control filters ##
########################################

### Filter SNPs with observed frequencies outside of expected range ###
# Add frequency range
cd $Fadm_statistic_calc_folder/Fadm
cut -f1 expected_freq.common_snps.txt > common.snps
awk 'FNR==NR{a[$1]=$1; next}; $1 in a {print $0;}' common.snps $Fadm_statistic_calc_folder/Expected_freq/minmax_freqs.allchrs.txt > minmax_freqs.common_snps.txt
rm common.snps

awk 'BEGIN { FS=OFS="\t" }
     NR==FNR { min_freq[$1]=$2; max_freq[$1]=$3; next }
     FNR==1 { print $0, "MIN_FREQ", "MAX_FREQ"; next }
     { print $0, min_freq[$3], max_freq[$3] }' minmax_freqs.common_snps.txt Fadm_score.txt > merged_output.txt
mv merged_output.txt Fadm_score.txt

# Filter by range
head -n+1 Fadm_score.txt > Fadm_score.FILTERED.txt
awk '$6 > $10' Fadm_score.txt | awk '$6 < $11' >> Fadm_score.FILTERED.txt


### Remove indels and complementary alleles (G/C, A/T) ###
head -n+1 Fadm_score.FILTERED.txt > Fadm_score.FILTERED2.txt
cat Fadm_score.FILTERED.txt | awk '($4 == "A") || ($4 == "G") || ($4 == "C") || ($4 == "T")' | awk '($5 == "A") || ($5 == "G") || ($5 == "C") || ($5 == "T")' | tr ' ' '\t' >> Fadm_score.FILTERED2.txt
mv Fadm_score.FILTERED2.txt Fadm_score.FILTERED.txt

head -n+1 Fadm_score.FILTERED.txt > Fadm_score.FILTERED2.txt
cat Fadm_score.FILTERED.txt | awk '!(($4 == "A") && ($5 == "T"))' | awk '!(($4 == "T") && ($5 == "A"))' | awk '!(($4 == "G") && ($5 == "C"))' | awk '!(($4 == "C") && ($5 == "G"))' | tr ' ' '\t' >> Fadm_score.FILTERED2.txt
mv Fadm_score.FILTERED2.txt Fadm_score.FILTERED.txt


### Filter by low frequency (expected and observed freq) ###
head -n+1 Fadm_score.FILTERED.txt > Fadm_score.FILTERED2.txt
awk '$9 < 0.99' Fadm_score.FILTERED.txt | awk '$9 > 0.01' | awk '$6 > 0.01' | awk '$6 < 0.99' >> Fadm_score.FILTERED2.txt
mv Fadm_score.FILTERED2.txt Fadm_score.FILTERED.txt


### Filter by low frequency (in specific source populations: min and max freqs) ###
head -n+1 Fadm_score.FILTERED.txt > Fadm_score.FILTERED2.txt
awk '$10 > 0.01' Fadm_score.FILTERED.txt | awk '$11 < 0.99' >> Fadm_score.FILTERED2.txt
mv Fadm_score.FILTERED2.txt Fadm_score.FILTERED.txt


### Filter by potential weird artefacts (weird MAF distance > 40%, potential allele switch or wrongly called marker?) ###
head -n+1 Fadm_score.FILTERED.txt > Fadm_score.FILTERED2.txt
cat Fadm_score.FILTERED.txt | awk '!($6 >= $9 + 0.4)' | awk '!($6 <= $9 - 0.4)' >> Fadm_score.FILTERED2.txt

mv Fadm_score.FILTERED2.txt Fadm_score.FILTERED.txt
mv Fadm_score.FILTERED.txt Fadm_score.SNP_selection.txt



##########
## Plot ##
##########
cd $Fadm_statistic_calc_folder/Fadm
echo "

    ### Packages ###
    library(ggplot2)
    library(dplyr)
    library(data.table)


    ###### Load data ######
    freq_deviation_file <- 'Fadm_score.SNP_selection.txt'
    FADM.allchrs <- fread(freq_deviation_file, header = T)


    ##############
    # Fadm SCORE #
    ##############

    ###### Prepare data for plot ######
    FADM.allchrs\$POS <- FADM.allchrs\$POS / 1000000
    FADM.allchrs\$CHR <- as.character(FADM.allchrs\$CHR)
    FADM.allchrs\$CHR <- factor(FADM.allchrs\$CHR, levels = as.character(c(1:22)))
    data_cum <- FADM.allchrs %>% 
    group_by(CHR) %>% 
    summarise(max_bp = max(POS)) %>% 
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(CHR, bp_add)

    FADM2 <- FADM.allchrs %>% 
    inner_join(data_cum, by = 'CHR') %>% 
    mutate(POS_cum = POS + bp_add)

    axis_set <- FADM2 %>% 
    group_by(CHR) %>% 
    summarize(center = mean(POS_cum))

    FADM2 <- na.omit(FADM2)
    FADM2 <- FADM2[FADM2\$FADM != Inf, ]
    FADM.allchrs <- FADM2


    ###### Fadm score plot ######
    FADM.allchrs_plot <- FADM.allchrs
    FADM.allchrs_plot <- FADM.allchrs_plot[FADM.allchrs_plot\$FADM > 0.10 ,]
    fadm_plot.allchrs <- ggplot(FADM.allchrs_plot, aes(x = POS_cum, y = (FADM), col=CHR))+
        geom_point()+
        scale_color_manual(values = rep(c('palevioletred4','palevioletred2'),11), guide='none')+
        theme_classic()+
        scale_x_continuous(label = axis_set\$CHR, breaks = axis_set\$center)+
        ylab(paste0('Fadm'))


    ###### Save plots ######
    ggsave(filename=paste0('../../Fadm.png'), fadm_plot.allchrs)

" > Rcommand.R
Rscript Rcommand.R
rm Rcommand.R




