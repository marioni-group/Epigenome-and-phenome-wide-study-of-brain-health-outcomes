############################################################################################
############################################################################################
############################### Interpretation of study results ############################
############################################################################################
############################################################################################


############################################################################################################

### EWAS results summaries

############################################################################################################

library(tidyverse)
library(readxl)

cpgs = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_090621/no_eGFR_cistrans_table_unformatted_naming.csv") 
list1 <- unique(cpgs$SeqId)

# Do a quick count of the total cpgs for each protein marker and the total cpg signatures across proteins 
EWAS_phen <- cpgs

# Of the proteome, which proteins had the most amount of CpG signals?
	Gene <- EWAS_phen[c("SeqId", "gene")]
	count <- Gene %>% group_by(SeqId) %>% count(SeqId)
	count <- as.data.frame(count)
	count <- count[order(-count$n),]
	count <- left_join(count, Gene, by = "SeqId")
	count <- unique(count)
	count <- count[c(1,3,2)]
	names(count)[3] <- "N"

# Which proteins had pQTLs in Sun et al 
list <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/pQTLs_seqId_edited.xlsx")
list <- as.data.frame(list)
list <- list[-1,]

# 1981 pQTLs from Sun et al - for 1561 unique SeqIds 
# unique SNPs = 

count2 <- left_join(count, list, by = "SeqId")

# Work out which proteins had SNPs in our extraction (837 SNPs)
# read in SNPs extracted 
SNPs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr_joint_to_single_file/pQTLs_170521.csv", check.names = F)

# read in list of SNPs with SeqIds for proteins in the list for extraction
table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/pQTLs_in_extraction_list.csv")

# subset the table to only those SNPs we have extracted for protein regression available 
names <- colnames(SNPs)[-1]
names <- as.data.frame(names)
names(names)[1] <- "pQTL"

join <- left_join(names, table, by = "pQTL") # we have 1690 associations with biomarkers accross the 837 SNPs extracted 
names(join)[2] <- "SeqId"

# Add this info on what SNPs were available to the main table of EWAS protein counts 
count3 <- left_join(count2, join, by = "SeqId")

# Format table 
count <- count3[c(1,2,3,9,10,11,14,15,16)]
names(count) <- c("SeqId", "Gene of Protein", "Number of CpG Associations", "Sun et al Sentinel Variant", "Sun et al Chromosome", "Sun et al Position", "Sun et al Effect Allele", "Sun et al Other Allele", "pQTL extracted in STRADL")

# Save out suppl table 
write.csv(count, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/no_eFGR_PROTEIN_summary_with_pQTLs_linked.csv"), row.names = F)


###################

### EWAS CPG CAT ANNO

library(tidyverse)
library(readxl)

cpgs = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_090621/no_eGFR_cistrans_table_unformatted_naming.csv") 
list1 <- unique(cpgs$SeqId)

# Do a quick count of the total cpgs for each protein marker and the total cpg signatures across proteins 
EWAS_phen <- cpgs

	# Which were the most common cpgs in the proteome EWAS signature?
	probe <- EWAS_phen[c("Probe")]
	count <- probe %>% group_by(Probe) %>% count(Probe)
	count <- as.data.frame(count)
	count <- count[order(-count$n),]
	names(count)[2] <- "No_times_selected_across_proteins"

# Annotate to EWAS catalogue lookup 

# names(count)[1] <- "CpG_site"
# counts3 <- count(count, CpG_site)

# counts3 <- counts3[order(-counts3$n),]

# Read in catalogue file 
EWAS <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins/outputs/Glmnet_re_runs_for_all_analyses/EWAS_catalogue/EWAS_Catalog_03-07-2019.txt")

results <- count

# Set a place for the annotations 
results$Annotation <- "X"

# Get annotations for the cpgs of interest 

for (i in 1:1760){
	cpg <- results[i,1]
	anno <- EWAS
	anno_cpg <- anno[which(anno$CpG %in% cpg),]
	anno_cpg <- anno_cpg[which(anno_cpg$P < 3.6e-8),]
	trait <- anno_cpg$Trait %>% unique()
	str <- str_c(trait, collapse = ", ")
	results[i,3] <- str
}

# Now we want to add a column which lists which of the proteins each cpg associates with 
# pred needs to be a file that has the weights with predictor info 
pred <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_090621/no_eGFR_cistrans_table_unformatted_naming.csv")
names(pred)[10] <- "Short_name"
pred <- pred %>% as.data.frame()

unique(pred$Probe) # 1760 probes 

for (i in 1:2895){ # get dimensons of the pred results file 
	cpg <- as.character(results[i,1])
	data <- pred
	data_cpg <- data[which(data$Probe %in% cpg),]
	list <- data_cpg$Short_name %>% unique()
	str <- str_c(list, collapse = ", ")
	results[i,"EpiScore"] <- str

	list2 <- data_cpg$SeqId %>% unique()
	str2 <- str_c(list2, collapse = ", ")
	results[i,"SeqId"] <- str2
}


# Add in gene of CpG and trans/cis info 
names(pred)[5] <- "CpG_gene"
add <- pred[,c(3,5)]

count <- left_join(results, add, by = "Probe")

# Save off results file with annotations for full EWAS 
count <- count[c(1,2,4,6,3)]
count <- unique(count)
write.csv(count, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/no_eGFR_CPG_summary_annotated.csv"), row.names = F)

# Do a version for the filtered neuro proteins 


############################################################################################################

### Find which proteins overlap in the pheWAS between imaging/cognitive/APOE phenotypes at FDR & BON 

############################################################################################################

### BONFERRONI FIRST 

### APOE 
apoe <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated/APOE_JOINT_annotated_covs_bon.csv")
apoe$type <- "APOE"
list1 <- apoe[which(apoe$Status == "pass"),]
list1 <- list1[c("SeqId", "phenotype", "type")] # 11 proteins for APOE 

### IMAGING 
brain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/brain_accel_result_annotated_bon.csv")
brain$type <- "imaging"
brain <- brain[which(brain$Status == "pass"),]
list1 <- rbind(list1, brain[c("SeqId", "phenotype", "type")])

GGM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/Global_GM_Volume_result_annotated_bon.csv")
GGM$type <- "imaging"
GGM <- GGM[which(GGM$Status == "pass"),]
list1 <- rbind(list1, GGM[c("SeqId", "phenotype", "type")])

GFA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/gFA_result_annotated_bon.csv")
GFA$type <- "imaging"
GFA <- GFA[which(GFA$Status == "pass"),]
list1 <- rbind(list1, GFA[c("SeqId", "phenotype", "type")])

GMD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/gMD_result_annotated_bon.csv")
GMD$type <- "imaging"
GMD <- GMD[which(GMD$Status == "pass"),]
list1 <- rbind(list1, GMD[c("SeqId", "phenotype", "type")])

# CER <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/Cerebrum_WM_Volume_result_annotated_bon.csv")
# CER$type <- "imaging"
# CER <- CER[which(CER$Status == "pass"),]
# list1 <- rbind(list1, CER[c("SeqId", "phenotype", "type")])

# FAZ <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/Fazekas_Score_Total_result_annotated_bon.csv")
# FAZ$type <- "imaging"
# FAZ <- FAZ[which(FAZ$Status == "pass"),]
# list1 <- rbind(list1, FAZ[c("SeqId", "phenotype", "type")])

WM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_WMHV_110821/WMHV_annotated_FDR.csv")
WM$type <- "imaging"
WM <- WM[which(WM$Status == "pass"),]
list1 <- rbind(list1, WM[c("SeqId", "phenotype", "type")])

WBV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/WBV_With_Ventricles_result_annotated_bon.csv")
WBV$type <- "imaging"
WBV <- WBV[which(WBV$Status == "pass"),]
list1 <- rbind(list1, WBV[c("SeqId", "phenotype", "type")])


### COG SCORES
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/g_result_annotated_bon.csv")
g$type <- "cognitive"
g <- g[which(g$Status == "pass"),]
list1 <- rbind(list1, g[c("SeqId", "phenotype", "type")])

gf <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/gf_result_annotated_bon.csv")
gf$type <- "cognitive"
gf <- gf[which(gf$Status == "pass"),]
list1 <- rbind(list1, gf[c("SeqId", "phenotype", "type")])

digit_symbol <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/digit_symbol_result_annotated_bon.csv")
digit_symbol$type <- "cognitive"
digit_symbol <- digit_symbol[which(digit_symbol$Status == "pass"),]
list1 <- rbind(list1, digit_symbol[c("SeqId", "phenotype", "type")])

LM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/LM_result_annotated_bon.csv")
LM$type <- "cognitive"
LM <- LM[which(LM$Status == "pass"),]
list1 <- rbind(list1, LM[c("SeqId", "phenotype", "type")])

mr_correct <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/mr_correct_result_annotated_bon.csv")
mr_correct$type <- "cognitive"
mr_correct <- mr_correct[which(mr_correct$Status == "pass"),]
list1 <- rbind(list1, mr_correct[c("SeqId", "phenotype", "type")])

verbal_total <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/verbal_total_result_annotated_bon.csv")
verbal_total$type <- "cognitive"
verbal_total <- verbal_total[which(verbal_total$Status == "pass"),]
list1 <- rbind(list1, verbal_total[c("SeqId", "phenotype", "type")])

vocabulary <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/vocabulary_result_annotated_bon.csv")
vocabulary$type <- "cognitive"
vocabulary <- vocabulary[which(vocabulary$Status == "pass"),]
list1 <- rbind(list1, vocabulary[c("SeqId", "phenotype", "type")])

# Work out how many assocs for each phenotype and modality 
# apoe is already 11 as seen above 

dim(list1) # 119 total associations for phenos 

apoe <- list1 %>% filter(type == "APOE")
dim(apoe) # 11 APOE 

cog <- list1 %>% filter(type == "cognitive")
dim(cog) # 92 cognitive 

im <- list1 %>% filter(type == "imaging")
dim(im) # 16 imaging

# Look at overlap between phenos 

which(apoe$SeqId %in% cog$SeqId) # 0 
which(apoe$SeqId %in% im$SeqId) # 0 
which(cog$SeqId %in% im$SeqId) # 4

overlap <- cog[which(cog$SeqId %in% im$SeqId),]

#      SeqId    phenotype      type
# 13 7638-30            g cognitive
# 30 7638-30           gf cognitive
# 48 15426-5           gf cognitive
# 56 15426-5 digit_symbol cognitive

overlap2 <- im[which(im$SeqId %in% cog$SeqId),]

#     SeqId   phenotype    type
# 5 15426-5 brain_accel imaging
# 8 7638-30 brain_accel imaging

# So the overlap here at bonferroni is 2 proteins 


##############################################################################

### FDR SECOND 

### APOE 
apoe <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated/APOE_JOINT_annotated_covs_FDR.csv")
apoe$type <- "APOE"
list1 <- apoe[which(apoe$Status == "pass"),]
list1 <- list1[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")] # 11 proteins for APOE 

### IMAGING 
brain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/brain_accel_result_annotated_FDR.csv")
brain$type <- "imaging"
brain$phenotype <- "Brain Age Acceleration"
brain <- brain[which(brain$Status == "pass"),]
list1 <- rbind(list1, brain[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

GGM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/Global_GM_Volume_result_annotated_FDR.csv")
GGM$type <- "imaging"
GGM$phenotype <- "Global Grey Matter Volume"
GGM <- GGM[which(GGM$Status == "pass"),]
list1 <- rbind(list1, GGM[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

GFA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/gFA_result_annotated_FDR.csv")
GFA$type <- "imaging"
GFA$phenotype <- "General Fractional Anisotropy"
GFA <- GFA[which(GFA$Status == "pass"),]
list1 <- rbind(list1, GFA[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

GMD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/gMD_result_annotated_FDR.csv")
GMD$type <- "imaging"
GMD$phenotype <- "General Mean Diffusivity"
GMD <- GMD[which(GMD$Status == "pass"),]
list1 <- rbind(list1, GMD[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

# CER <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/Cerebrum_WM_Volume_result_annotated_FDR.csv")
# CER$type <- "imaging"
# CER$phenotype <- "Cerebrum White Matter Volume"
# CER <- CER[which(CER$Status == "pass"),]
# list1 <- rbind(list1, CER[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

# FAZ <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/Fazekas_Score_Total_result_annotated_FDR.csv")
# FAZ$type <- "imaging"
# FAZ$phenotype <- "Fazekas White Matter Hyperintensity Score"
# FAZ <- FAZ[which(FAZ$Status == "pass"),]
# list1 <- rbind(list1, FAZ[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

WM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_WMHV_110821/WMHV_annotated_FDR.csv")
WM$type <- "imaging"
WM <- WM[which(WM$Status == "pass"),]
list1 <- rbind(list1, WM[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

WBV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/WBV_With_Ventricles_result_annotated_FDR.csv")
WBV$type <- "imaging"
WBV$phenotype <- "Whole Brain Volume"
WBV <- WBV[which(WBV$Status == "pass"),]
list1 <- rbind(list1, WBV[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])


### COG SCORES
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/g_result_annotated_FDR.csv")
g$type <- "cognitive"
g$phenotype <- "General Cognitive Ability"
g <- g[which(g$Status == "pass"),]
list1 <- rbind(list1, g[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

gf <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/gf_result_annotated_FDR.csv")
gf$type <- "cognitive"
gf$phenotype <- "General Fluid Cognitive Ability"
gf <- gf[which(gf$Status == "pass"),]
list1 <- rbind(list1, gf[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

digit_symbol <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/digit_symbol_result_annotated_FDR.csv")
digit_symbol$type <- "cognitive"
digit_symbol$phenotype <- "Processing Speed"
digit_symbol <- digit_symbol[which(digit_symbol$Status == "pass"),]
list1 <- rbind(list1, digit_symbol[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

LM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/LM_result_annotated_FDR.csv")
LM$type <- "cognitive"
LM$phenotype <- "Logical Memory"
LM <- LM[which(LM$Status == "pass"),]
list1 <- rbind(list1, LM[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

mr_correct <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/mr_correct_result_annotated_FDR.csv")
mr_correct$type <- "cognitive"
mr_correct$phenotype <- "Non-Verbal Reasoning"
mr_correct <- mr_correct[which(mr_correct$Status == "pass"),]
list1 <- rbind(list1, mr_correct[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

verbal_total <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/verbal_total_result_annotated_FDR.csv")
verbal_total$type <- "cognitive"
verbal_total$phenotype <- "Verbal Reasoning"
verbal_total <- verbal_total[which(verbal_total$Status == "pass"),]
list1 <- rbind(list1, verbal_total[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

vocabulary <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated_cog/vocabulary_result_annotated_FDR.csv")
vocabulary$type <- "cognitive"
vocabulary$phenotype <- "Vocabulary"
vocabulary <- vocabulary[which(vocabulary$Status == "pass"),]
list1 <- rbind(list1, vocabulary[c("SeqId", "phenotype", "type", "Gene.Name.Name", "beta", "SE", "Pcalc", "UniProt.Full.Name")])

# Work out how many assocs for each phenotype and modality 
# apoe is already 11 as seen above 

dim(list1) # 644 total associations for phenos 

apoe <- list1 %>% filter(type == "APOE")
dim(apoe) # 11 APOE 

cog <- list1 %>% filter(type == "cognitive")
dim(cog) # 579 cognitive 

im <- list1 %>% filter(type == "imaging")
dim(im) # 54 imaging


# Look at overlap between phenos 

which(apoe$SeqId %in% cog$SeqId) # 2

overlap <- apoe[which(apoe$SeqId %in% cog$SeqId),]

#       SeqId phenotype type
# 4  17671-58      APOE APOE - ING4
# 11  2797-56      APOE APOE - APOB

overlap2 <- cog[which(cog$SeqId %in% apoe$SeqId),]

#        SeqId    phenotype      type
# 80   2797-56            g cognitive - APOB
# 203  2797-56           gf cognitive - APOB
# 333  2797-56 digit_symbol cognitive - APOB
# 559 17671-58   mr_correct cognitive - ING4


which(apoe$SeqId %in% im$SeqId) # 0 

which(cog$SeqId %in% im$SeqId) # 59

overlap3 <- cog[which(cog$SeqId %in% im$SeqId),]

# 59 proteins 

overlap4 <- im[which(im$SeqId %in% cog$SeqId),]

# 29 proteins 

# > unique(overlap4$SeqId)
#  [1] "15470-11"  "4763-31"   "15426-5"   "7638-30"   "17691-1"   "6379-62"
#  [7] "4811-33"   "8352-26"   "3343-1"    "11178-21"  "3616-3"    "5307-12"
# [13] "4876-32"   "11516-7"   "9235-3"    "3316-58"   "6478-2"    "15573-110"
# [19] "15468-14"  "9266-1"    "13565-2"   "15522-2"   "4374-45"   "14337-1"
# [25] "6480-1"    "5452-71"

# So 26 unique proteins overlap between imaging and cognitive variables in FDR 

# Try to make a table to show the overlapping proteins 

overlap3 <- overlap3[-3]
names(overlap3)[2] <- "Cognitive"

overlap4 <- overlap4[-3]
names(overlap4)[2] <- "Imaging"

join <- merge(overlap3, overlap4, by = "SeqId", all = TRUE)

# Look at whether direction of effect is same or different 
join$Direction <- ifelse(join$beta.y < 0 & join$beta.x < 0 | join$beta.y > 0 & join$beta.x > 0, "Matched", "Opposite")

# Look at which is present in our EWAS with signal 
cpgs = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_090621/no_eGFR_cistrans_table_unformatted_naming.csv") 
cpgs1 <- cpgs[1] %>% unique()
cpgs1$EWAS <- "present" # 153

join <- left_join(join, cpgs1, by = "SeqId")

# sort out naming and order for table save out 
join2 <- join[c(1,3,13,2,8,10,11,12,4,5,6,14,15)]
join2 <- join2 %>% group_by(Gene.Name.Name.x)
join2 <- join2[order(-join2$Pcalc.y),]
join2 <- as.data.frame(join2)
names(join2) <- c("SeqId", "Gene of Protein", "UniProt Full Name", "Cognitive Phenotype", "Imaging Phenotype", "Imaging Beta", "Imaging SE", "Imaging P", "Cognitive Beta", "Cognitive SE", "Cognitive P", "Direction of Effect", "EWAS Signal")
write.csv(join2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/overlap_protein_markers_644_26_unique_ordered_2.csv", row.names = F)


# Look at the full 644 and see how many are in our EWAS dataset 
length(which(cpgs1$SeqId %in% list1$SeqId)) # 26 proteins passing FDR had EWAS signal 
cpgs <- cpgs[which(cpgs$SeqId %in% list1$SeqId),] 
dim(cpgs) # 88 pQTMs across the 26 proteins 

trans <- cpgs[which(cpgs$Effect == "TRANS"),]
dim(trans) # 42

cis <- cpgs[which(cpgs$Effect == "CIS"),]
dim(cis) # 46

write.csv(cpgs, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/no_eGFR_cis_FDR_88_assocs_slice_filtered.csv", row.names = F)


# 16 trans and 11 cis proteins involved here 
# > unique(trans$gene)
#  [1] "SMPD1"    "CRTAM"    "SCUBE1"   "RBL2"     "CSF1R"    "GBP1"
#  [7] "HEXB"     "TREM2"    "MX1"      "ALDOB"    "PSAT1"    "IL18BP"
# [13] "TNFRSF1B" "ACY1"     "CST5"     "CD163"


> unique(trans$gene)
 [1] "SMPD1"  "CRTAM"  "SCUBE1" "RBL2"   "CSF1R"  "GBP1"   "HEXB"   "TREM2"
 [9] "MX1"    "ALDOB"  "PSAT1"  "IL18BP" "ACY1"   "CD163"

# So we lose TNFRSF1B from the trans associations 

# > unique(cis$gene)
#  [1] "CHI3L1"   "LILRA3"   "CRELD1"   "SIGLEC14" "IL18R1"   "SIGLEC5"
#  [7] "MX1"      "UGDH"     "FAHD1"    "OLFM2"    "PCSK7"


# Write out cis and trans associations that passed FDR and were neurologically-linked 

write.csv(trans, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/no_eGFR_trans_FDR_42_assocs_filtered.csv", row.names = F)

write.csv(cis, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/no_eGFR_cis_FDR_46_assocs_filtered.csv", row.names = F)

# Write out a formatted table with the pQTMs that are selected for neurological markers overall 
chrom <- cpgs
chrom <- chrom %>% group_by(chrom$Gene.Name.Name)
chrom <- as.data.frame(chrom)
# chrom <- chrom[order(chrom$p),]
chrom <- chrom[c(3,2,5,4,6:9,1,10,11,12,15:17,20)]
names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
  "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
# Add in the SNP info 
table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/PROTEIN_summary_with_pQTLs_linked.csv")
table <- table[c(1,4,9)]
join <- left_join(chrom, table, by = "SeqId")
write.csv(join, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/no_eGFR_neuro_slice_cistrans_table_pQTL_added.csv", row.names = F)

