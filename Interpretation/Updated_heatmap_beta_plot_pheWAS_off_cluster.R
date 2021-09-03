###################################################################################################################

### DO A CROSS-NEURO PHENOTYPE PLOT FOR BONFERRONI SIGNIFICANCE 

###################################################################################################################

# Read in bonferroni files and get the list of SeqIds which were significant across all neurological markers 

## First make a list of seqIds which pass across all phenotypes 

### APOE 
apoe <- read.csv("Y:/Danni/Cluster_Files/pheWAS_scaled_310521/results_scaled_annotated/APOE_JOINT_annotated_covs_bon.csv")
list1 <- apoe[which(apoe$Status == "pass"),]
list1 <- list1[c(1)] # 11 proteins for APOE 

setwd("Y:/Danni/Cluster_Files/phewas/")

### IMAGING 
brain <- read.csv("brain_accel_result_annotated_bon.csv")
brain <- brain[which(brain$Status == "pass"),]
list1 <- rbind(list1, brain[c(1)])

GGM <- read.csv("Global_GM_Volume_result_annotated_bon.csv")
GGM <- GGM[which(GGM$Status == "pass"),]
list1 <- rbind(list1, GGM[c(1)])

GFA <- read.csv("gFA_result_annotated_bon.csv")
GFA <- GFA[which(GFA$Status == "pass"),]
list1 <- rbind(list1, GFA[c(1)])

GMD <- read.csv("gMD_result_annotated_bon.csv")
GMD <- GMD[which(GMD$Status == "pass"),]
list1 <- rbind(list1, GMD[c(1)])

# CER <- read.csv("Cerebrum_WM_Volume_result_annotated_bon.csv")
# CER <- CER[which(CER$Status == "pass"),]
# list1 <- rbind(list1, CER[c(1)])

FAZ <- read.csv("Y:/Danni/Cluster_Files/pheWAS_WMHV_110821/WMHV_annotated_FDR.csv")
FAZ <- FAZ[which(FAZ$Status == "pass"),]
list1 <- rbind(list1, FAZ[c(1)])

WBV <- read.csv("WBV_With_Ventricles_result_annotated_bon.csv")
WBV <- WBV[which(WBV$Status == "pass"),]
list1 <- rbind(list1, WBV[c(1)])

setwd("Y:/Danni/Cluster_Files/pheWAS_scaled_310521/results_scaled_annotated_cog/")

### COG SCORES
g <- read.csv("g_result_annotated_bon.csv")
g <- g[which(g$Status == "pass"),]
list1 <- rbind(list1, g[c(1)])

gf <- read.csv("gf_result_annotated_bon.csv")
gf <- gf[which(gf$Status == "pass"),]
list1 <- rbind(list1, gf[c(1)])

digit_symbol <- read.csv("digit_symbol_result_annotated_bon.csv")
digit_symbol <- digit_symbol[which(digit_symbol$Status == "pass"),]
list1 <- rbind(list1, digit_symbol[c(1)])

LM <- read.csv("LM_result_annotated_bon.csv")
LM <- LM[which(LM$Status == "pass"),]
list1 <- rbind(list1, LM[c(1)])

mr_correct <- read.csv("mr_correct_result_annotated_bon.csv")
mr_correct <- mr_correct[which(mr_correct$Status == "pass"),]
list1 <- rbind(list1, mr_correct[c(1)])

verbal_total <- read.csv("verbal_total_result_annotated_bon.csv")
verbal_total <- verbal_total[which(verbal_total$Status == "pass"),]
list1 <- rbind(list1, verbal_total[c(1)])

vocabulary <- read.csv("vocabulary_result_annotated_bon.csv")
vocabulary <- vocabulary[which(vocabulary$Status == "pass"),]
list1 <- rbind(list1, vocabulary[c(1)])


## Next, read the files in and chop down to key info for plotting 
apoe <- read.csv("Y:/Danni/Cluster_Files/pheWAS_scaled_310521/results_scaled_annotated/APOE_JOINT_annotated_covs_bon.csv")

setwd("Y:/Danni/Cluster_Files/phewas/")
brain <- read.csv("brain_accel_result_annotated_bon.csv")
GGM <- read.csv("Global_GM_Volume_result_annotated_bon.csv")
GFA <- read.csv("gFA_result_annotated_bon.csv")
GMD <- read.csv("gMD_result_annotated_bon.csv")
# CER <- read.csv("Cerebrum_WM_Volume_result_annotated_bon.csv")
FAZ <- read.csv("Y:/Danni/Cluster_Files/pheWAS_WMHV_110821/WMHV_annotated_FDR.csv")
WBV <- read.csv("WBV_With_Ventricles_result_annotated_bon.csv")

setwd("Y:/Danni/Cluster_Files/pheWAS_scaled_310521/results_scaled_annotated_cog/")
g <- read.csv("g_result_annotated_bon.csv")
gf <- read.csv("gf_result_annotated_bon.csv")
digit_symbol <- read.csv("digit_symbol_result_annotated_bon.csv")
LM <- read.csv("LM_result_annotated_bon.csv")
mr_correct <- read.csv("mr_correct_result_annotated_bon.csv")
verbal_total <- read.csv("verbal_total_result_annotated_bon.csv")
vocabulary <- read.csv("vocabulary_result_annotated_bon.csv")

# Reverse polarity of the gfa and gmd variable results for beta column only 
GFA$beta <- GFA$beta*(-1)
GMD$beta <- GMD$beta*(-1)

g <- g[c(1,9,3,5,7)]
g$phenotype <- "General Cognitive Ability"
gf <- gf[c(1,9,3,5,7)]
gf$phenotype <- "General Fluid Cognitive Ability"
digit_symbol <- digit_symbol[c(1,9,3,5,7)]
digit_symbol$phenotype <- "Processing Speed"
LM <- LM[c(1,9,3,5,7)]
LM$phenotype <- "Logical Memory"
mr_correct <- mr_correct[c(1,9,3,5,7)]
mr_correct$phenotype <- "Non-Verbal Reasoning"
verbal_total <- verbal_total[c(1,9,3,5,7)]
verbal_total$phenotype <- "Verbal Reasoning"
vocabulary <- vocabulary[c(1,9,3,5,7)]
vocabulary$phenotype <- "Vocabulary"

bind <- rbind(verbal_total, vocabulary)
bind <- rbind(bind, LM)
bind <- rbind(bind, mr_correct)
bind <- rbind(bind, digit_symbol)
bind <- rbind(bind, gf)
bind <- rbind(bind, g)

apoe <- apoe[c(1,9,3,5,7)]
apoe$phenotype <- "APOE"

bind <- rbind(bind, apoe)

brain <- brain[c(1,9,3,5,7)]
brain$phenotype <- "Relative Brain Age"
GGM <- GGM[c(1,9,3,5,7)]
GGM$phenotype <- "Global Grey Matter Volume"
GFA <- GFA[c(1,9,3,5,7)]
GFA$phenotype <- "General Fractional Anisotropy"
GMD <- GMD[c(1,9,3,5,7)]
GMD$phenotype <- "General Mean Diffusivity"
# CER <- CER[c(1,9,3,5,7)]
# CER$phenotype <- "Cerebrum White Matter Volume"
FAZ <- FAZ[c(1,2,3,5,7)]
FAZ$phenotype <- "White Matter Hyperintensity Volume"
WBV <- WBV[c(1,9,3,5,7)]
WBV$phenotype <- "Whole Brain Volume"

bind <- rbind(bind, brain)
bind <- rbind(bind, GGM)
bind <- rbind(bind, GFA)
bind <- rbind(bind, GMD)
# bind <- rbind(bind, CER)
bind <- rbind(bind, FAZ)
bind <- rbind(bind, WBV)

# Order by P value to see most significant associations 
bind <- bind[order(bind$Pcalc),]

# Set to another variable for common verison 
bind2 <- bind

# Index the top 100 associations by significance of p value as cutoff for * insertion
bind$stars <- cut(bind$Pcalc, breaks=c(-Inf, 0.00001296099, 0.01, 0.05, Inf), label=c("*", "", "", ""))

# Subset to the SeqIds in the top 100 assocs
bind <- bind[which(bind$SeqId %in% bind$SeqId[1:100]),]

# # Restrict to the proteins of interest across all phenotypes from above 
# bind <- bind[which(bind$SeqId %in% list1$SeqId),]

# # Add stars for significance (bonferroni)
# bind$stars <- cut(bind$Pcalc, breaks=c(-Inf, 0.0000180638, 0.01, 0.05, Inf), label=c("*", "", "", ""))

# Make genes unique 
names(bind)[2] <- "gene"

length(unique(bind$SeqId)) # 66 proteins 

# Plot 
library(ggplot2)


pdf("Y:/Danni/Cluster_Files/heatmap_plot_betas_joint_250821_top_100_assocs_GFA_flipped.pdf", width = 14.5, height = 5)
ggplot(aes(x=gene, y=phenotype, fill=beta), data=bind) +
 geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
 geom_text(aes(label=stars), color="black", size=5) +
theme(axis.text.x = element_text(size=10, angle=45, hjust = 0.9),
          axis.text.y = element_text(size=12, angle=0)) + xlab("") + ylab("") +
scale_y_discrete(limits = c("APOE", "Relative Brain Age", "Global Grey Matter Volume", "General Fractional Anisotropy",
  "General Mean Diffusivity", "White Matter Hyperintensity Volume", "Whole Brain Volume",
  "General Cognitive Ability", "General Fluid Cognitive Ability", "Logical Memory", "Processing Speed", "Non-Verbal Reasoning", "Vocabulary",
  "Verbal Reasoning")) + guides(fill=guide_legend(title="Beta effect"))
dev.off()


###################################################################################################

### Do version with common proteins across imaging and cognitive 

###################################################################################################

# Read in the common associations table 
library(readxl)
common <- read_excel("Y:/Danni/stradl_markers/Rev_Co_author_circulation/old/Table_common_markers.xlsx")
common <- as.data.frame(common)

# Read in GMD and gDA
file <- read_excel("Y:/Danni/stradl_markers/Rev_Co_author_circulation/Suppl_tables/gFA_gMD.xlsx")
file <- as.data.frame(file)

# Reverse signs for all 3 files
file[,3] <- file[,3]*(-1)
file[,9] <- file[,9]*(-1)

# Save out for suppl files 
setwd("Y:/Danni/stradl_markers/Rev_Co_author_circulation/Suppl_tables/")
write.csv(file, "GFS_GMD_flipped.csv", row.names = F)

# Subset bind file to the 26 proteins with the common associations of inetrest 
bind3 <- bind2[which(bind$SeqId %in% common$SeqId),]

# Make genes unique 
names(bind3)[2] <- "gene"

length(unique(bind3$SeqId)) # 26 proteins 

# Plot 
library(ggplot2)


pdf("Y:/Danni/Cluster_Files/heatmap_plot_betas_joint_250821_common_26_assocs_GFA_flipped_unannotated.pdf", width = 14.5, height = 5)
ggplot(aes(x=gene, y=phenotype, fill=beta), data=bind3) +
 geom_tile() + scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") +
scale_y_discrete(limits = c("APOE", "Relative Brain Age", "Global Grey Matter Volume", "General Fractional Anisotropy",
  "General Mean Diffusivity", "White Matter Hyperintensity Volume", "Whole Brain Volume",
  "General Cognitive Ability", "General Fluid Cognitive Ability", "Logical Memory", "Processing Speed", "Non-Verbal Reasoning", "Vocabulary",
  "Verbal Reasoning")) +
theme(axis.text.x = element_text(size=10, angle=45, hjust = 0.9),
          axis.text.y = element_text(size=12, angle=0)) + xlab("") + ylab("")  + guides(fill=guide_legend(title="Beta effect"))
dev.off()




###################################################################################################

### Plot by haplotype for top proteins for APOE 

###################################################################################################

# low="#2C7BB6", mid="white", high="#D7191C"

library(tidyverse)
# Load in results files for APOE
APOE <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated/APOE_JOINT_annotated_rank_inverse_Pcalc_bon.csv")
# APOE_covs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated/APOE_JOINT_COVS_annotated_rank_inverse_Pcalc_bon.csv")

# Get the names of the top proteins passing FDR in both cases
APOE <- APOE %>% filter(Status == "pass")

data <- APOE[c(1,9)] # get just seqID and gene names for proteins 

data <- data %>% unique() # 11 unique protein levels associated with either e2 or e4 allele status 

# Read in the protein data 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/file_for_heatmaps_cog.csv", check.names = F)

# Plot for each of the proteins in relation to the main protein measurements and haplotypes 
library(ggpubr)
library(ggplot2)

pred <- data

plot_list <- list()
plot_list2 <- list()

### Individual plots using code above as basis, so that we can patch work them together 

# 3  10082-251           NEFL #3 - blue
i <- 3
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  d1 <- prot[,"APOE"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d4 <- cbind(d1,d2)
  names(d4) <- c("APOE", "prot")
  d4$APOE <- as.factor(d4$APOE)
  d4 <- na.omit(d4)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

  q = ggplot(d4, aes(x=APOE , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

  plot_list[[i]] <- p
  plot_list2[[i]] <- q

   pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/", protein_name, "_", seqid_name, "_plot_sig_proteins_by_hapltotype.pdf"))
  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))
  dev.off()

NEFL <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))

# 6   10046-55          BIRC2 #6 - blue 
i <- 6
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  d1 <- prot[,"APOE"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d4 <- cbind(d1,d2)
  names(d4) <- c("APOE", "prot")
  d4$APOE <- as.factor(d4$APOE)
  d4 <- na.omit(d4)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

  q = ggplot(d4, aes(x=APOE , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

  plot_list[[i]] <- p
  plot_list2[[i]] <- q

   pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/", protein_name, "_", seqid_name, "_plot_sig_proteins_by_hapltotype.pdf"))
  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))
  dev.off()

BIRC2 <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))

# 7    19158-1            PAF #7 - blue
i <- 7
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  d1 <- prot[,"APOE"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d4 <- cbind(d1,d2)
  names(d4) <- c("APOE", "prot")
  d4$APOE <- as.factor(d4$APOE)
  d4 <- na.omit(d4)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

  q = ggplot(d4, aes(x=APOE , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

  plot_list[[i]] <- p
  plot_list2[[i]] <- q

   pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/", protein_name, "_", seqid_name, "_plot_sig_proteins_by_hapltotype.pdf"))
  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))
  dev.off()

PAF <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))

# 11   2797-56           APOB #11 - red
i <- 11
seqid_name <- as.character(data[i,1])
  protein_name <- as.character(data[i,2])

  d1 <- prot[,"apoe"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("apoe", "prot")
  d3 <- na.omit(d3)

  d1 <- prot[,"APOE"] %>% as.data.frame()
  d2 <- prot[,seqid_name] %>% as.data.frame()
  d4 <- cbind(d1,d2)
  names(d4) <- c("APOE", "prot")
  d4$APOE <- as.factor(d4$APOE)
  d4 <- na.omit(d4)

  p = ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

  q = ggplot(d4, aes(x=APOE , y=prot)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal() + ggtitle(protein_name)

  plot_list[[i]] <- p
  plot_list2[[i]] <- q

   pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/", protein_name, "_", seqid_name, "_plot_sig_proteins_by_hapltotype.pdf"))
  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='lightskyblue3', color = 'lightskyblue3')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))
  dev.off()

APOB <-  ggplot(d3, aes(x=apoe , y=prot)) +
  geom_violin(trim=FALSE, fill='firebrick2', color = 'firebrick2')+
  geom_boxplot(width=0.1) + theme_minimal() + xlab("APOE Haplotype") + ylab(protein_name) + theme_classic() + 
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0)))

###################################################################################################

### Plot associations for interesting overlapping proteins between cog/imaging 

###################################################################################################

library(tidyverse)
library(ggplot2)
library(ggpubr)

data <-  read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/prot_file_150621.csv", check.names = F)

# Which assocs do we want to show 

# RBL2 - digit, gf, global grey matter vol
# HEXB 
# SMPD1
# ACY1

prot <- data

# RBL2 gf

  d1 <- prot[,"13565-2"] %>% as.data.frame()
  d2 <- prot[,"gf"] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("prot", "trait")
  d3 <- na.omit(d3)


  pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/imaging_plots/RBL2_gf.pdf"))
  ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("RBL2") + ylab("General Fluid Cognitive Ability")
  dev.off()


RBL2_gf <- ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("RBL2") + ylab("General Fluid Cognitive Ability")

# RBL2 GGM

  d1 <- prot[,"13565-2"] %>% as.data.frame()
  d2 <- prot[,"Global_GM_Volume"] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("prot", "trait")
  d3 <- na.omit(d3)


  pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/imaging_plots/RBL2_GGM.pdf"))
  ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("RBL2") + ylab("Global Grey Matter Volume")
  dev.off()

RBL2_GGM <-   ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("RBL2") + ylab("GGlobal Grey Matter Volume")


# TREM1 - 9266-1 - digit, non verbal, Global Grey Matter Volume

  d1 <- prot[,"9266-1"] %>% as.data.frame()
  d2 <- prot[,"Global_GM_Volume"] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("prot", "trait")
  d3 <- na.omit(d3)


  pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/imaging_plots/TREM1_GGM.pdf"))
  ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("TREM1") + ylab("Global Grey Matter Volume")
  dev.off()

TREM1_GGM <-   ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("TREM1") + ylab("Global Grey Matter Volume")

# TREM1 digit 

  d1 <- prot[,"9266-1"] %>% as.data.frame()
  d2 <- prot[,"digit_symbol"] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("prot", "trait")
  d3 <- na.omit(d3)


  pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/imaging_plots/TREM1_digit.pdf"))
  ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("TREM1") + ylab("Processing Speed")
  dev.off()

 TREM1_digit <- ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("TREM1") + ylab("Processing Speed")



# HEXB - brain age accel 

  d1 <- prot[,"15470-11"] %>% as.data.frame()
  d2 <- prot[,"brain_accel"] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("prot", "trait")
  d3 <- na.omit(d3)


  pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/imaging_plots/HEXB_BAA.pdf"))
  ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = 'firebrick2') +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("HEXB") + ylab("Brain Age Acceleration")
  dev.off()

HEXB_brain <-   ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = 'firebrick2') +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("HEXB") + ylab("Brain Age Acceleration")

# HEXB digit symbol 

  d1 <- prot[,"15470-11"] %>% as.data.frame()
  d2 <- prot[,"digit_symbol"] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("prot", "trait")
  d3 <- na.omit(d3)


  pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/imaging_plots/HEXB_digit.pdf"))
  ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("HEXB") + ylab("Processing Speed")
  dev.off()

HEXB_digit <-   ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("HEXB") + ylab("Processing Speed")


# NEU1 - g 15426-5

  d1 <- prot[,"15426-5"] %>% as.data.frame()
  d2 <- prot[,"g"] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("prot", "trait")
  d3 <- na.omit(d3)


  pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/imaging_plots/NEU1_g.pdf"))
  ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("NEU1") + ylab("g")
  dev.off()

NEU1_g <-   ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "lightskyblue3") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("NEU1") + ylab("g")


# NEU1 - gFA 15426-5

  d1 <- prot[,"15426-5"] %>% as.data.frame()
  d2 <- prot[,"gFA"] %>% as.data.frame()
  d3 <- cbind(d1,d2)
  names(d3) <- c("prot", "trait")
  d3 <- na.omit(d3)


  pdf(file = paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/plotting_results/imaging_plots/NEU1_gFA.pdf"))
  ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "firebrick2") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("NEU1") + ylab("General Fractional Anisotropy")
  dev.off()

NEU1_gFA <-   ggscatter(d3, x="prot" , y="trait", add = "reg.line", conf.int = TRUE, color = "firebrick2") +
  theme(legend.title = element_blank(),
    legend.position = "none", 
    axis.title.x = element_text(size = 27), 
    axis.text.x = element_text(size = 27),
    axis.text.y = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    plot.title = element_text(size = 27),
    axis.title.x.bottom = element_text(margin = margin(14, 0, 0, 0))) + xlab("NEU1") + ylab("General Fractional Anisotropy")


# Now patchwork them all together 
library(patchwork)

part1 <- heatmap 

part2 <- NEFL + BIRC2 + PAF + APOB

part3 <- RBL2_gf + RBL2_GGM + TREM1_digit + TREM1_GGM

part4 <- HEXB_digit + HEXB_brain + NEU1_g + NEU1_gFA

# pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/plot_for_NEURO_joint_V2.pdf", width = 30, height = 30)
# heatmap / (NEFL + BIRC2 + PAF + APOB + RBL2_gf + RBL2_GGM + TREM1_digit + TREM1_GGM + HEXB_digit + HEXB_brain + NEU1_g + NEU1_gFA) +
# plot_layout(heights = c(10,20))
# dev.off()

# Save out part two separately joinT
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/plot_for_NEURO_part_2.pdf", width = 30, height = 20)
NEFL + BIRC2 + PAF + APOB + RBL2_gf + RBL2_GGM + TREM1_digit + TREM1_GGM + HEXB_digit + HEXB_brain + NEU1_g + NEU1_gFA
dev.off()

# Save out part one separately joint
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/plot_for_NEURO_part_1.pdf", width = 17, height = 5)
heatmap
dev.off()

# Save out the rows for part two separately
# First APOE row 
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/plot_for_NEURO_APOE_row.pdf", width = 30, height = 7.5)
NEFL + BIRC2 + PAF + APOB + 
plot_layout(ncol = 4)
dev.off()


# Save out the rows for part two separately
# Next row with two example proteins on it for cog/imaging
pdf("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/plot_for_NEURO_COGIM_row.pdf", width = 30, height = 7.5)
RBL2_gf + RBL2_GGM + TREM1_digit + TREM1_GGM + 
plot_layout(ncol = 4)
dev.off()













