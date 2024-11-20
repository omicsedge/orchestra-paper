rm(list=ls())

### Packages ###
library(ggplot2)
library(cowplot)
library(magrittr)
library(data.table)



###################################
## PREPARE ORCHESTRA LAI RESULTS ## 
###################################

### Pathways (change accordingly) ###
analysis_folder <- '/home/jonlerga/Escritorio/SelfDecode/Testing/AncestryPaper/NATURE_COMMUNICATIONS/REVIEWERS_COMMENTS/ToyExample/Analysis_folder/'
orchestra_results_folder <- paste0(analysis_folder, '/LAI_results')


### Get LAI results from Orchestra ###
setwd(orchestra_results_folder)

smooth_res_files <- dir() %>%
  .[!grepl('\\.txt$|\\.png$', .)] %>%
  gsub('smooth_samples\\.chr', '', .) %>%
  gsub('\\.tsv\\.gz$', '', .) %>%
  .[order(as.numeric(gsub('\\..*', '', .)))]


### Merge chromosomes ###
for(chrs in smooth_res_files) {
  file_path <- sprintf("smooth_samples.chr%s.tsv.gz", chrs)
  dt <- fread(file_path, skip = 1)
  
  dt <- dt[, .(V2, V3, V4, V5, V6, ID_Combined = paste(V4, V5, sep = "_"))]
  
  if(chrs == '1.2'){
    fwrite(dt, "LAIpredictions.txt", append = FALSE, col.names = FALSE, sep = "\t") 
  }else{
    fwrite(dt, "LAIpredictions.txt", append = TRUE, col.names = FALSE, sep = "\t")
  }
}



####################################
## ANCESTRY PERCENTAGE PER WINDOW ## 
####################################

### Load LAI data ###
LAIpred <- read.table(paste0("./LAIpredictions.txt"), header = T)
colnames(LAIpred) <- c("SampleID", "Hap", "Chr", "Window", "Ancestry", "Window_code")


### Filter sample IDs if needed ###
LAIpred <- LAIpred[grep('NA', LAIpred$SampleID),] # 1KGP samples


### Average estimated genome-wide ancestry proportions ###
global_perc <- data.frame( Ancestry = names(table(LAIpred$Ancestry)),
                  Count = as.numeric(table(LAIpred$Ancestry)) )
global_perc$Percentage <- global_perc$Count / sum(global_perc$Count)
write.table(global_perc, file = "../Admixture_proportions.txt", row.names = F, col.names = T, sep = "\t", quote = F)


### Select ancestry label to get percentage per window ###
# table(LAIpred$Ancestry)
ancestry_label <- 'EUR' # Here, European ancestry for Mexicans, but it could be something else for other admixed populations

windows_all <- as.character(unique(LAIpred$Window_code))
window_counter <- 1
for(window in windows_all){
  
  # Get data for window w (LAI info from all individuals)
  rows_aux <- LAIpred[LAIpred$Window_code == window,]
  chr <- unique(rows_aux$Chr)
  
  # Ancestry percentage in window w
  perc_anc <- nrow(rows_aux[which(rows_aux$Ancestry %in% ancestry_label), ]) / nrow(rows_aux) * 100

  # Store results
  if(window_counter == 1){
    df_selection_test <- data.frame(window = window_counter,
                          chr = paste0("chr",chr),
                          Perc = perc_anc)
    
  }else{
    df_selection_test <- rbind(df_selection_test,
                     data.frame(window = window_counter,
                                chr = paste0("chr",chr),
                                Perc = perc_anc))
  }
  window_counter <- window_counter + 1
}
df_selection_test$Source_population <- ancestry_label



##############################
## LAD STATISTIC PER WINDOW ## 
##############################
#
# alpha_w_p: admixture proportion from source population p for a given window w 
# alpha_mean_p: estimated genome-wide admixture proportion (we use alpha_mean_p per chromosome pack)
#

### Define chromosome packs used in Orchestra ###
df_selection_test$LAD <- -999 # Dummy value to be replaced
df_chr_pack <- data.frame(pack = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,10,10),
                          chrs = c(1:22))


### LAD score ###
for( chr_pack in unique(df_chr_pack$pack) ){
  
  # Chromosomes in chromosome pack
  chrs <- df_chr_pack[df_chr_pack$pack == chr_pack,]$chrs
  chrs <- paste0("chr",chrs)
  
  # Alpha_mean_p (= average estimated admixture proportion across all analyzed windows)
  alpha_mean_p <- mean(df_selection_test[which(df_selection_test$chr %in% chrs),]$Perc)
  
  # Alpha_w_p (= estimated admixture proportion for a specific window)
  alpha_w_p <- df_selection_test[which(df_selection_test$chr %in% chrs),]$Perc
  
  # LAD score
  df_selection_test[which(df_selection_test$chr %in% chrs),]$LAD <- alpha_w_p - alpha_mean_p

}


### Save LAD scores ###
write.table(df_selection_test,
            file = paste0("Selection_results.Source_populations.txt"),
            row.names = F, col.names = T, sep = "\t", quote = F)



##########
## PLOT ##
##########

### Load data ###
LAD_info <- read.table('Selection_results.Source_populations.txt', header=T)


### Ancestry label and LAD standard deviation ###
dataset <- 'Mexicans (1KGP)'
ancestry_label <- as.character(unique(LAD_info$Source_population))
sd <- sd(LAD_info$LAD)

LAD_info$chr <- factor(LAD_info$chr, levels = as.character(paste0("chr", c(1:22))))


### Adaptive signal in Cuadros-Espinoza (2022) ###
Plot.Percentage <- ggplot(LAD_info, aes(x = window, y = Perc, col=chr))+
  geom_vline(xintercept = 759, col="darkseagreen2", linewidth=2.5)

Plot.LAD <- ggplot(LAD_info, aes(x = window, y = LAD, col = chr))+
  geom_vline(xintercept = 759, col="darkseagreen2", linewidth=2.5)


### Plots ###

# Percentage
Plot.Percentage <- Plot.Percentage +
  geom_point()+
  scale_color_manual(values = rep(c("palevioletred4","palevioletred2"),11), guide="none")+
  theme_classic()+
  ylab(paste0("% ancestry\n(", ancestry_label, ")"))+
  ggtitle(paste0('|| ',dataset, ' ||\nPercentage'))

# LAD score
Plot.LAD <- Plot.LAD +
  geom_hline(yintercept = c(0, 3 * sd, -3 * sd), linetype=c(1,2,2), col="grey80", linewidth=1)+
  geom_point()+
  scale_color_manual(values = rep(c("palevioletred4","palevioletred2"),11), guide="none")+
  theme_classic()+
  ylab(paste0("LAD statistic\n(", ancestry_label, ")"))+
  ggtitle('LAD score')

# Final plot
Plot.Final <- plot_grid(Plot.Percentage, Plot.LAD, nrow = , ncol = 1)
print(Plot.Final)


### Save plots ###
ggsave(filename=paste0("../LADscore.png"), Plot.Final)



