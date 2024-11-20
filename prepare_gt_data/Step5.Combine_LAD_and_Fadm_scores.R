rm(list=ls())

### Packages ###
library(ggplot2)
library(dplyr)
library(data.table)


### Pathways ###
dataset <- 'Admixed_Mexicans'

analysis_folder <- ''
lad_score_folder <- paste0(analysis_folder, '/LAI_results')
fadm_score_folder <- paste0(analysis_folder, '/Fadm_estimation/Fadm')

setwd(analysis_folder)



#####################################################################################################################
######################## Calculate empirical P value per SNP for Fadm and LAD (individually) ########################
#####################################################################################################################

##################################
## Load LAD and Fadm statistics ##
##################################
FADM <- read.table(paste0(fadm_score_folder, '/Fadm_score.SNP_selection.txt'), header = T, sep = '\t')
LAD <- read.table(paste0(lad_score_folder, '/Selection_results.Source_populations.txt'), header = T, sep = '\t')


##############################################
## Get empirical P value for Fadm statistic ##
##############################################

###### Estimate Fadm empirical P value ######
FADM <- FADM[order(FADM$FADM, decreasing = T), ]
FADM$rank <- rank(-FADM$FADM, ties.method = "min")
FADM$empirical_p <- FADM$rank / nrow(FADM)


###### Regions data ######
# Ancestry regions defined for our custom panel (uploaded in GitHub). This can be defined for each SNP set and window size used when running Orchestra
windows_regions <- read.table("./Ancestry_regions.hg38.txt", header=F)

# Extend first and last windows per chromosome to cover the whole genome
colnames(windows_regions) <- c("Chr","Start","Stop")
for(chr in c(1:22)){
  lengths <- data.frame(Chr=c(1:22),
                        length=c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622,
                                 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468))
  
  windows_regions[windows_regions$Chr ==paste0("chr", chr),][1,2] <- 0
  windows_regions[windows_regions$Chr ==paste0("chr", chr),][nrow(windows_regions[windows_regions$Chr ==paste0("chr", chr),]),3] <- lengths[lengths$Chr == chr,2]
}

# Close gaps between windows
regions_new <- windows_regions
for(chr in c(1:22)){
  gap <- regions_new[regions_new$Chr == paste0("chr", chr), ]$Start[-1] - regions_new[regions_new$Chr == paste0("chr", chr), ]$Stop[-nrow(regions_new[regions_new$Chr == paste0("chr", chr), ])]
  regions_new2 <- data.frame(
    Chr = regions_new[regions_new$Chr == paste0("chr", chr), ]$Chr,
    Start = regions_new[regions_new$Chr == paste0("chr", chr), ]$Start - c(0,gap/2),
    Stop = regions_new[regions_new$Chr == paste0("chr", chr), ]$Stop + c(gap/2,0)
  )
  regions_new[regions_new$Chr == paste0("chr", chr), ] <- regions_new2
}
regions_new$Chr <- gsub("chr", "", regions_new$Chr)
regions_new$Window <- c(1:nrow(regions_new))


###### Assign windows to SNPs (Fadm) ######
regions_new$Chr <- as.numeric(regions_new$Chr)
regions_new$Start <- as.numeric(regions_new$Start)
regions_new$Stop <- as.numeric(regions_new$Stop)

FADM$window <- 0
for(chr in c(1:22)){
  breakpoints <- c(regions_new[regions_new$Chr == chr,]$Start[1], regions_new[regions_new$Chr == chr,]$Stop)
  potential_windows <- cut(FADM[FADM$CHR == chr,]$POS, breakpoints, labels = FALSE)
  potential_windows <- potential_windows + regions_new[regions_new$Chr == chr,4][1] - 1
  FADM[FADM$CHR == chr,]$window <- potential_windows
}


#############################################
## Get empirical P value for LAD statistic ##
#############################################

###### Check big imbalances in extreme portions of chromosomes ######
# Convert each chromosome into the chromosome pack
aux_regions <- unique(LAD[,c(2,1)])
aux_regions <- aux_regions[order(aux_regions$window),]
aux_regions <- as.data.table(aux_regions)
aux_regions$chr <- as.character(aux_regions$chr)

chr_lookup <- c(
  "chr1" = "chr1_2", "chr2" = "chr1_2",
  "chr3" = "chr3_4", "chr4" = "chr3_4",
  "chr5" = "chr5_6", "chr6" = "chr5_6",
  "chr7" = "chr7_8", "chr8" = "chr7_8",
  "chr9" = "chr9_10", "chr10" = "chr9_10",
  "chr11" = "chr11_12", "chr12" = "chr11_12",
  "chr13" = "chr13_14", "chr14" = "chr13_14",
  "chr15" = "chr15_16", "chr16" = "chr15_16",
  "chr17" = "chr17_18", "chr18" = "chr17_18",
  "chr19" = "chr19_22", "chr20" = "chr19_22", 
  "chr21" = "chr19_22", "chr22" = "chr19_22"
)
aux_regions[, window_pack := chr_lookup[chr]]
aux_regions$chr <- as.numeric(gsub("chr", "", aux_regions$chr))

# Get chromosome tips (5 regions for each start and end chromosome tips)
setorder(aux_regions, chr, window)
result <- aux_regions[, .SD[c(1:5, (.N-4):.N)], by = window_pack]

# Check if LAD has to be adjusted at chromosome tips
chr_tips.LAD <- LAD[which(LAD$window %in% result$window),]$LAD
chr_in.LAD <- LAD[!(LAD$window %in% result$window),]$LAD
pval_diff_in_LAD <- t.test(abs(chr_tips.LAD), abs(chr_in.LAD))$p.val

if(pval_diff_in_LAD < 0.05){
  sd.tips <- sd(chr_tips.LAD)
  sd.in <- sd(chr_in.LAD)
  
  ratio_adj <- sd.in / sd.tips
  LAD[which(LAD$window %in% result$window),]$LAD <- LAD[which(LAD$window %in% result$window),]$LAD * ratio_adj
}


###### Get empirical P value for LAD ######
LAD$LAD_abs <- abs(LAD$LAD)

colnames(FADM) <- c("CHR", "POS", "SNP_ID", "REF", "ALT",
                    "OBSERVED_FREQ", "OBS_COUNTS", "FADM", "EXPECTED_FREQ", "MIN_FREQ", "MAX_FREQ",
                    "RANK_FADM", "EMPIRICAL_P_FADM", "WINDOW_FOR_LAD")
FADM$LAD_ABS <- LAD[FADM$WINDOW_FOR_LAD, ]$LAD_abs

CombineStats.df <- FADM
CombineStats.df <- CombineStats.df[order(CombineStats.df$LAD_ABS, decreasing = T), ]
CombineStats.df$RANK_LAD <- rank(-CombineStats.df$LAD_ABS, ties.method = "min")
CombineStats.df$EMPIRICAL_P_LAD <- CombineStats.df$RANK_LAD / nrow(CombineStats.df)


###### Prepare analysis dataset ######
CombineStats.df <- CombineStats.df[,c(1,2,3,4,5,6,9,8,10,11,12,13,14,15,16,17)]

CombineStats.df$POS <- CombineStats.df$POS / 1000000
CombineStats.df$CHR <- as.character(CombineStats.df$CHR)
CombineStats.df$CHR <- factor(CombineStats.df$CHR, levels = as.character(c(1:22)))

data_cum <- CombineStats.df %>% 
  group_by(CHR) %>% 
  summarise(max_bp = max(POS)) %>% 
  mutate(BPADD = lag(cumsum(max_bp), default = 0)) %>% 
  select(CHR, BPADD)

CombineStats.df_intermediate <- CombineStats.df %>% 
  inner_join(data_cum, by = "CHR") %>% 
  mutate(POS_CUM = POS + BPADD)

axis_set <- CombineStats.df_intermediate %>% 
  group_by(CHR) %>% 
  summarize(center = mean(POS_CUM))

CombineStats.df_intermediate <- na.omit(CombineStats.df_intermediate)
CombineStats.df_intermediate <- CombineStats.df_intermediate[CombineStats.df_intermediate$FADM != Inf, ]
CombineStats.df <- CombineStats.df_intermediate


###################################################################################################################################
######################## Combine LAD and Fadm statistics and get combined P value per SNP (manhattan plot) ########################
###################################################################################################################################

###### Calculate Chi-squared distribution and estimate combined P values using Fadm and LAD statistics ######
# We approximate significance threshold = 0.05 / (1992 * 2)
significance_threshold = 10 ** (-6)

CombineStats.df$chi_squared <- (-2) * (log(CombineStats.df$EMPIRICAL_P_FADM) + log(CombineStats.df$EMPIRICAL_P_LAD))
CombineStats.df$combined_Pval <- pchisq(CombineStats.df$chi_squared, df = 4, lower.tail = FALSE)

axis_set <- CombineStats.df %>% 
  group_by(CHR) %>% 
  summarize(center = mean(POS_CUM))

###### Position of selection signal for adaptive admixed found in Cuadros-Espinoza (2022) for Mexicans ######  
CuadrosEspinoza_signal_position <- median(CombineStats.df[CombineStats.df$WINDOW_FOR_LAD == 759,]$POS_CUM)


###### Plot ######  
CombineStats.df_plot <- CombineStats.df[CombineStats.df$combined_Pval < 0.05, ]
CombineStats.df_plot$CHR <- as.character(CombineStats.df_plot$CHR)
CombineStats.df_plot$CHR <- factor(CombineStats.df_plot$CHR, levels = as.character(c(1:22)))

CombineStats.df_plot$significant <- "no"
CombineStats.df_plot[CombineStats.df_plot$combined_Pval < significance_threshold, ]$significant <- "yes"

ManhattanPlot <- ggplot(CombineStats.df_plot, aes(x = POS_CUM, y = -log10(combined_Pval), col=CHR, size = significant))+
  geom_vline(xintercept = CuadrosEspinoza_signal_position, col="darkseagreen2", linewidth=5)+      ## Signals per population in CuadrosEspinoza (2022)
  geom_point()+
  geom_hline(yintercept = c(-log10(significance_threshold)), 
             linetype = 2, col="coral")+    
  scale_color_manual(values = rep(c("#385b9d","#70a0e0"),11), guide="none")+
  theme_classic()+
  theme(text = element_text(size = 16))+
  scale_size_manual(values = c(1,2.5), guide = "none")+
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center)+
  ylab(paste0("-log10(P value)"))+
  xlab("Chromosome")
ManhattanPlot


###### Save plots ######  
setwd(analysis_folder)
ggsave(filename=paste0("GWsignals_selection.", dataset, ".png"), ManhattanPlot)


