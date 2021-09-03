####################################################################################
####################################################################################
############################## STRADL markers ######################################
####################################################################################
####################################################################################

# Prepping using Robs script as basis to keep parameters as consistent as possible
# He based it on Riccardo's initial processing of the STRADL data 

# ADDITION: this version regresses pQTLs from protein levels, to adjust for known
# genetic effects.

####################################################################################

# Merge target from W2/W3 and subset proteins to those in target file 

####################################################################################

## read in the somalogic data ##
setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/")

soma <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/GS+soma+QC+normalized.csv")

## read in Archie's master linker file ##
link <- read.csv("ST id linkage.csv")

## read in the wave 2 target file ##
w2_tar <- read.csv("ALL-wave2.csv") 
a <- grep("S", w2_tar$Sample_Name)
w2_tar <- w2_tar[a,1:4]
w2_tar$id <- gsub("S", "", w2_tar$Sample_Name)

## merge target file with link file using GS id ##
w2_tar_update <- merge(link, w2_tar, by="id")

## filter to those who passed somalogic qc ##
b <- which(w2_tar_update$st %in% soma$SampleId)
w2_tar_soma <- w2_tar_update[b,]

## filter to stradl id (same as proteomics id), gs id, and DNAm id ##
w2_target <- w2_tar_soma[,c("st","id","Sample_Sentrix_ID", "sex", "age")]
names(w2_target) <- c("Stradl_id","GS_id","DNAm_id", "sex", "age_stradl")

## read in wave 3 target file ##
w3_tar <- read.csv("ST.csv")
w3_tar$id <- gsub("ST", "", w3_tar$Sample_Name)

## merge with link file ##
w3_tar_update <- merge(link, w3_tar, by="id")

## filter to those who passed somalogic qc ##
b <- which(w3_tar_update$st %in% soma$SampleId)
w3_tar_soma <- w3_tar_update[b,]

## filter to stradl id (same as proteomics id), gs id, and DNAm id ##
w3_target <- w3_tar_soma[,c("st","id","Sample_Sentrix_ID", "sex", "age")]
names(w3_target) <- c("Stradl_id","GS_id","DNAm_id", "sex", "age_stradl")

## harmonise into a single file ##
w3_target$wave <- "w3"
w2_target$wave <- "w2"

stradl_DNAm_target <- rbind(w2_target, w3_target)

# #write.table(stradl_DNAm_target, file="U:/Datastore/IGMM/marioni-lab/STRADL/STRADL_DNAm_target_REM_17April2020.txt", quote=F, col.names=T, row.names=F, sep="\t")
# # i will write a copy of this to my prep folder so i know how it was generated for the 844 individuals:
# write.table(stradl_DNAm_target, file="/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/STRADL_DNAm_target_REM_17April2020_regen_140121.txt", quote=F, col.names=T, row.names=F, sep="\t")


# soma1 <- soma[which(soma$SampleId %in% stradl_DNAm_target$Stradl_id),]
# write.csv(soma1, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/soma_844_file.csv", row.names = F)

####################################################################################

## Read in w2 and w3 cell count info and add into target

####################################################################################

## WAVE 2 
stradl_w2 <- stradl_DNAm_target[stradl_DNAm_target$wave %in% "w2",]
w2_cell <- read.csv("wave2_cell.csv")
w2_cell$Sample_Name <- gsub(".*S", "", w2_cell$Sample_Name)
w2_cell$Sample_Name <- gsub(".*R", "", w2_cell$Sample_Name)
w2_cell <- w2_cell[-which(duplicated(w2_cell$Sample_Name)),]
stradl_w2 <- merge(stradl_w2, w2_cell, by.x="GS_id", by.y="Sample_Name", all.x=T)
for(i in 7:12){ 
stradl_w2[,i][stradl_w2[,i] %in% NA] <- mean(stradl_w2[,i],na.rm = T)
} 

write.csv(stradl_w2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_w2_cells.csv", row.names = F)


## WAVE 3 
stradl_w3 <- stradl_DNAm_target[stradl_DNAm_target$wave %in% "w3",]
w3_cell <- read.csv("w3_cell.csv")
w3_cell$Sample_Name <- gsub(".*ST", "", w3_cell$Sample_Name)
w3_cell <- w3_cell[-which(duplicated(w3_cell$Sample_Name)),]
stradl_w3 <- merge(stradl_w3, w3_cell, by.x = "GS_id", by.y="Sample_Name", all.x=T)

write.csv(stradl_w3, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_w3_cells.csv", row.names = F)


####################################################################################

## Read in w2 and w3 batch info and add into target

####################################################################################

## WAVE 2 

w2_batch <- read.csv("wave2_batch.csv")
stradl_w2 <- merge(stradl_w2, w2_batch, by.x="GS_id", by.y="Sample_Name", all.x=T)

## WAVE 3 

w3_batch <- read.csv("wave_3_original_samplesheet.csv")
w3_batch$Sample_Name <- gsub(".*ST", "", w3_batch$Sample_Name)
w3_batch = w3_batch[-which(duplicated(w3_batch$Sample_Name)),]
stradl_w3 <- merge(stradl_w3, w3_batch[,c("Sample_Name", "Batch_all")], by.x="GS_id",by.y="Sample_Name",all.x=T)
stradl_w3$Batch_all[stradl_w3$Batch_all == 1] <- 32
stradl_w3$Batch_all[stradl_w3$Batch_all == 3] <- 33
stradl_w3$Batch_all[stradl_w3$Batch_all == 4] <- 34
stradl_w3$Batch_all[stradl_w3$Batch_all == 5] <- 35
stradl_w3$Batch_all[stradl_w3$Batch_all == 6] <- 36
names(stradl_w3)[13] <- "Batch"
stradl <- rbind(stradl_w2, stradl_w3)

nrow(stradl) # this file has the WBCs and the batch for the 844 STRADL people

# Write out batch info 
batch <- stradl[c(1:3,13)]

write.csv(batch, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_batch.csv", row.names = F)

####################################################################################

## Add covariates into target file 

####################################################################################

soma_demo <- merge(link[,c("st","id", "age","sex")], soma, by.y = "SampleId",by.x="st") # merging in proteins with demographics 
demo <- read.csv("demographicsV2.csv")

soma_demo1 <- merge(soma_demo, demo[,c("st","st_age", "sample_dt")])
table(soma_demo1$age == soma_demo1$st_age)
soma_demo1$st_age <- NULL
library(stringi)
soma_demo1$study_site <- stri_sub(soma_demo1$st, 5,6)
soma_demo1$month_of_sample <- substring(soma_demo1$PlateRunDate, 4,5)
soma_demo1$month_of_draw <- substring(soma_demo1$sample_dt, 4,5)
soma_demo1$year_of_sample <- substring(soma_demo1$PlateRunDate, 7,10)
soma_demo1$year_of_draw <- substring(soma_demo1$sample_dt, 7,10)
soma_demo1$year_lag <- as.numeric(soma_demo1$year_of_sample) - as.numeric(soma_demo1$year_of_draw)
soma_demo1$month_lag <- as.numeric(soma_demo1$month_of_sample) - as.numeric(soma_demo1$month_of_draw)
soma_demo1$month_lag <- soma_demo1$month_lag/12
soma_demo1$lag_time <- soma_demo1$year_lag + soma_demo1$month_lag

eigen <- read.table("GS20K_ALL_MAF5_PCA.eigenvec")
head(eigen)
eigen <- eigen[which(eigen$V2 %in% soma_demo1$id),]
eigen$V1 <- NULL
names(eigen)[1] <- "id"
names(eigen)[2:21] <- paste0("PC", 1:20)
soma_demo1 <- merge(eigen, soma_demo1, by= "id")


####################################################################################

## Determine Association between Lag Time and Protein Variables 

####################################################################################

list1<- list() 
for(i in colnames(soma_demo1)[56:4639]){ 
  list1[[i]]<- summary(lm(soma_demo1[,i] ~ soma_demo1$lag_time))$coefficients[2,c(1,4)]
}

# linear regressions run above - each protein has been regressed to show its relationship with lag time variable 
# estimates and Pr(>|t|) values listed 

list2 <- do.call("rbind", list1)
list2 <- as.data.frame(list2)
list2$Protein <- row.names(list2)
names(list2)[1] <- "Beta"
names(list2)[2] <- "P.Value"
list3 <- list2[which(list2$P.Value < 0.05/4584),]
list3=list3[order(list3$P.Value),]
list3[list3$Beta <0 , "Protein"]


soma_demo1$lag_group <- 0
soma_demo1[soma_demo1$lag_time <= 2, "lag_group"] <- 1
soma_demo1[soma_demo1$lag_time > 2 & soma_demo1$lag_time <= 2.9, "lag_group"] <- 2
soma_demo1[soma_demo1$lag_time > 2.9 & soma_demo1$lag_time <= 3.5, "lag_group"] <- 3
soma_demo1[soma_demo1$lag_time > 3.5, "lag_group"] <- 4

# classify the associations into grouped brackets based on results 


####################################################################################

## Add in Creatinine and Estimate Glomerular Filtration Rate 

####################################################################################

creat <- read.csv("CRE corrected.csv")
names(creat)[5] <- "Creat"
soma_demo1 <- merge(creat[,c(1,5)], soma_demo1, by = "id", all.y = T)

females = soma_demo1[which(soma_demo1$sex %in% "F"),]
males = soma_demo1[which(soma_demo1$sex %in% "M"),] 

females_low <- females[which(females$Creat <= 62),] 
females_high <- females[which(females$Creat > 62),] 
females_na <- females[which(females$Creat %in% NA),] 

males_low <- males[which(males$Creat <= 80),] 
males_high <- males[which(males$Creat > 80),] 
males_na <- males[which(males$Creat %in% NA),] 

females_low$eGFR <- 141*((females_low$Creat/61.9)^-0.329)*(0.993^females_low$age)*1.018
females_high$eGFR <- 141*((females_high$Creat/61.9)^-1.209)*(0.993^females_high$age)*1.018
males_low$eGFR <- 141*((males_low$Creat/79.6)^-0.411)*(0.993^males_low$age)
males_high$eGFR <-141*((males_high$Creat/79.6)^-1.209)*(0.993^males_high$age)
females_na$eGFR <- NA
males_na$eGFR <- NA 

females_1 = rbind(females_low,females_high)
females_1 = rbind(females_1, females_na)
males_1 = rbind(males_low,males_high)
males_1 = rbind(males_1, males_na)
all = rbind(males_1, females_1)
soma_demo1 <- all

dim(soma_demo1) #  1065 4302

# write a copy out for the pheWAS eGFR correction 
# eGFR <- soma_demo1[c(1,4302)]
# write.csv(eGFR, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/eGFR_in_1065.csv", row.names = F)


set <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/Set_info_778.csv")
soma_demo2 = merge(soma_demo1, set, by = "id")
matcher = set$id
soma_demo2 = soma_demo2[match(matcher, soma_demo2$id), ]

id_ewas <- read.csv("IDs_778.csv")
id_gwas <- read.csv("IDs_1064.csv")

## Log Transform 
soma_demo2 <- soma_demo2[soma_demo2$id %in% id_ewas$x,]
ids = id_ewas$x 
soma_demo2 <- soma_demo2[match(ids,soma_demo2$id), ]
phenotypes <- soma_demo2
for(i in colnames(phenotypes)[57:4291]){ 
  phenotypes[,i]<- log(phenotypes[,i])
}

# Read in the pQTLs extracted for the sun et al list 
QTLs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/chr_extraction_imputed/run_all_chr_joint_to_single_file/pQTLs_170521.csv", check.names = F)
names(QTLs)[1] <- "id"

library(tidyverse)

# Join QTLs into the dataset 
phenotypes <- left_join(phenotypes, QTLs, by = "id")

# Read in the sun pQTL data 
library(data.table)
library(tidyverse)
library(readxl)
list <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/pQTL_list.xlsx")
list <- as.data.frame(list)

# format the seq ID column appropriately 
names(list)[2] <- "SOMAmer"
list$SeqId <- gsub("..$", "", list$SOMAmer)
library(stringr)
list$SeqId <- list$SeqId %>% str_replace_all("\\.", "-")
list$SeqId <- gsub('^.+?-(.*)', "\\1",list$SeqId)
write.csv(list, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/pQTLs_seqId.csv", row.names = F)

# Edit to make sure all names have converted SeqId correctly and read back in
list <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_pQTL_extraction_GS/pQTLs_seqId_edited.xlsx")
list <- as.data.frame(list)

# recreate order2 from pQTL extraction - which was the correct format for allele order 
order2 <- list[c(7,8,12,11,13)]
order2 <- order2[-1,]
order2 <- order2 %>% unite("new", 1:4, remove = FALSE)
order2 <- order2[c(6,1:5)] # 1980 protein - pQTL associations 
length(unique(order2$SeqId)) # 1561 unique proteins 

# Now format a table which can be called upon below for protein & sites 
table <- order2[c(1,2)]
names(table) <- c("Biomarker", "pQTL") # format of the pQTLs now matches the QTLs extracted and added to phenotype dataset here above for regressions

# Subset the table to pQTLs which are available in the dataset 
ov <- which(table$pQTL %in% colnames(phenotypes))
table <- table[ov,]

# write out tabe for joining in pQTL record suppl table mapping 
write.csv(table, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/pQTLs_in_extraction_list.csv", row.names = F)

# Now - due to the formula input not recognising _ below, we need to remove all _ from pQTL formats in table index and in the pQTL columns
table$pQTL <- gsub("_", "", table$pQTL)
table$pQTL <- sub("^", "A", table$pQTL)

# Format pQTLs without _'s so that they can be inputted into regression formula below
colnames(QTLs) <- gsub("_", "", colnames(QTLs))
colnames(QTLs) <- sub("^", "A", colnames(QTLs))

# Replace naming structure for proteins with seqids 
seq <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/260121_Matching_cohorts/STRADL_778_residualised_matched.csv", check.names = F)
colnames(phenotypes)[57:4291] <- colnames(seq[2:4236])


############################################################################################

## Regress Proteins onto Covariates - with eGFR included - AND ADDITIONAL SNP REGRESSION FOR THE PQTLS FOR EACH PROTEIN WITH KNOWN GENETIC ASSOCIATIONS
phenotypes_residualised <- phenotypes
phenotypes_residualised$eGFR[phenotypes_residualised$eGFR %in% NA] <- mean(phenotypes_residualised$eGFR, na.rm = T)

# Now in order to run the regressions we will need to rename the seqIds so they are without "-" separators, then rename after results run
table$Biomarker <- gsub("-", "", table$Biomarker)
table$Biomarker <- sub("^", "A", table$Biomarker)
names_saved <- colnames(phenotypes_residualised)[57:4291]
names_saved2 <- colnames(phenotypes_residualised)[4304:5140]
# sort the naming of proteins to match table 
colnames(phenotypes_residualised)[57:4291] <- gsub("-", "", colnames(phenotypes_residualised)[57:4291])
colnames(phenotypes_residualised)[57:4291] <- sub("^", "A", colnames(phenotypes_residualised)[57:4291])
#sort the naming of pQTLs to match table 
colnames(phenotypes_residualised)[4304:5140] <- gsub("_", "", colnames(phenotypes_residualised)[4304:5140])
colnames(phenotypes_residualised)[4304:5140] <- sub("^", "A", colnames(phenotypes_residualised)[4304:5140])

# Do a test protein to make sure loop for regressions is working 
duplicated(table$Biomarker)
# ov <- which(table$Biomarker == "9216-100")
# table  <- table[ov,]
# which(colnames(phenotypes_residualised) == "9216-100") # 4035 
i <- colnames(phenotypes_residualised)[4035]

for(i in colnames(phenotypes_residualised)[57:4291]){ 
	name <- as.character(names(phenotypes_residualised[i])) # index the protein name 
	sites <- table[which(table$Biomarker %in% name),] # get the sites (if any from the table)
	if (nrow(sites) == "0"){# most proteins wont have any matches in the table and can be residualised as per usual 
		phenotypes_residualised[,i]<- lm(phenotypes_residualised[,i] ~ age + factor(sex) + factor(study_site) + factor(lag_group) + eGFR + 
                                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, 
                                    na.action = na.exclude, data = phenotypes_residualised)$residuals
	}else {
		pQTL <- sites$pQTL # get the list of pQTLs that youll need to pull out to regress
		formula = paste0(names(phenotypes_residualised[i]), " ~ age + factor(sex) + factor(study_site) + factor(lag_group) + eGFR + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + ", paste0(pQTL, collapse=" + "))
		phenotypes_residualised[,i] <- as.numeric(phenotypes_residualised[,i])
		phenotypes_residualised[,i] <- scale(resid(lm(formula, data = phenotypes_residualised, na.action="na.exclude")))
	}
}

## Rank-Inverse Based Normaliation on Pre-corrected Protein Phenotypes
library(bestNormalize)
for(i in colnames(phenotypes_residualised)[57:4291]){ 
  phenotypes_residualised[,i]<- orderNorm(phenotypes_residualised[,i])$x.t
}

## Scale somalogic data 
phenotypes_residualised[,57:4291] <- apply(phenotypes_residualised[,57:4291], 2, scale)

## Replace naming of the proteins with the unformatted SeqIds as saved above 
names(phenotypes_residualised[57:4291]) <- names_saved

# Check 
colMeans(phe[3:4237])
apply(phe, 2, sd, na.rm = T)

### SAVE PHENO WITH EGFR + pQTLs regressed out 

# ## Save out somalogic file which has all covariates 
# write.csv(phenotypes_residualised, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_all_covariates_778_eGFR_inlcuded_pQTLs_regressed.csv", row.names = F)

## Save out somalogic file whih has no covariates 
phe <- phenotypes_residualised[,c(1,23,57:4291)]	
colnames(phe)[3:4237] <- names_saved
write.csv(phe, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778_eGFR_inlcuded_pQTLs_regressed.csv", row.names = F)



### Do a correlation check between SMPD1 and HEXB

prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778_eGFR_inlcuded_pQTLs_regressed.csv", check.names = F)

list <- c("10818-36", "15470-11")

prots <- prot[,which(colnames(prot) %in% list)]
names(prots) <-  c("SMPD1", "HEXB")

cor.test(prots$SMPD1, prots$HEXB)

# t = 34.817, df = 776, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.7518092 0.8068412
# sample estimates:
#       cor
# 0.7808354




###############################################################################################################

### ALTERNATE VERSION: REMOVE eGFR AND REGRESS AGAIN FOR RE-RUNS 

## Regress Proteins onto Covariates - without eGFR included - AND ADDITIONAL SNP REGRESSION FOR THE PQTLS FOR EACH PROTEIN WITH KNOWN GENETIC ASSOCIATIONS
phenotypes_residualised <- phenotypes
phenotypes_residualised$eGFR[phenotypes_residualised$eGFR %in% NA] <- mean(phenotypes_residualised$eGFR, na.rm = T)

# Now in order to run the regressions we will need to rename the seqIds so they are without "-" separators, then rename after results run
table$Biomarker <- gsub("-", "", table$Biomarker)
table$Biomarker <- sub("^", "A", table$Biomarker)
names_saved <- colnames(phenotypes_residualised)[57:4291]
names_saved2 <- colnames(phenotypes_residualised)[4304:5140]
# sort the naming of proteins to match table 
colnames(phenotypes_residualised)[57:4291] <- gsub("-", "", colnames(phenotypes_residualised)[57:4291])
colnames(phenotypes_residualised)[57:4291] <- sub("^", "A", colnames(phenotypes_residualised)[57:4291])
#sort the naming of pQTLs to match table 
colnames(phenotypes_residualised)[4304:5140] <- gsub("_", "", colnames(phenotypes_residualised)[4304:5140])
colnames(phenotypes_residualised)[4304:5140] <- sub("^", "A", colnames(phenotypes_residualised)[4304:5140])

# Do a test protein to make sure loop for regressions is working 
duplicated(table$Biomarker)
# ov <- which(table$Biomarker == "9216-100")
# table  <- table[ov,]
# which(colnames(phenotypes_residualised) == "9216-100") # 4035 
i <- colnames(phenotypes_residualised)[4035]

for(i in colnames(phenotypes_residualised)[57:4291]){ 
	name <- as.character(names(phenotypes_residualised[i])) # index the protein name 
	sites <- table[which(table$Biomarker %in% name),] # get the sites (if any from the table)
	if (nrow(sites) == "0"){# most proteins wont have any matches in the table and can be residualised as per usual 
		phenotypes_residualised[,i]<- lm(phenotypes_residualised[,i] ~ age + factor(sex) + factor(study_site) + factor(lag_group) + 
                                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, 
                                    na.action = na.exclude, data = phenotypes_residualised)$residuals
	}else {
		pQTL <- sites$pQTL # get the list of pQTLs that youll need to pull out to regress
		formula = paste0(names(phenotypes_residualised[i]), " ~ age + factor(sex) + factor(study_site) + factor(lag_group) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + ", paste0(pQTL, collapse=" + "))
		phenotypes_residualised[,i] <- as.numeric(phenotypes_residualised[,i])
		phenotypes_residualised[,i] <- scale(resid(lm(formula, data = phenotypes_residualised, na.action="na.exclude")))
	}
}

# subset to one protein for comparison between non pQTL version 
J <- phenotypes_residualised[,which(colnames(phenotypes_residualised) %in% "A141514")]

# > head(J)
#         [,1]
# 1  0.3950379
# 2 -0.9600435
# 3 -0.2341473
# 4  0.1991286
# 5 -0.8692787
# 6 -0.5309943


# # Test to make sure pQTL regression is working - comment out after 
# for(i in colnames(phenotypes_residualised)[57:4291]){ 
# 	name <- as.character(names(phenotypes_residualised[i])) # index the protein name 
# 	sites <- table[which(table$Biomarker %in% name),] # get the sites (if any from the table)
# 	if (nrow(sites) == "0"){# most proteins wont have any matches in the table and can be residualised as per usual 
# 		phenotypes_residualised[,i]<- lm(phenotypes_residualised[,i] ~ age + factor(sex) + factor(study_site) + factor(lag_group) + 
#                                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, 
#                                     na.action = na.exclude, data = phenotypes_residualised)$residuals
# 	}else {
# 		pQTL <- sites$pQTL # get the list of pQTLs that youll need to pull out to regress
# 		formula = paste0(names(phenotypes_residualised[i]), " ~ age + factor(sex) + factor(study_site) + factor(lag_group) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")
# 		phenotypes_residualised[,i] <- as.numeric(phenotypes_residualised[,i])
# 		phenotypes_residualised[,i] <- scale(resid(lm(formula, data = phenotypes_residualised, na.action="na.exclude")))
# 	}
# }

# # subset to one protein for comparison between non pQTL version 
# J <- phenotypes_residualised[,which(colnames(phenotypes_residualised) %in% "A141514")]

# # > head(J)
# #         [,1]
# # 1  0.1243286
# # 2 -0.8340210
# # 3 -0.1797023
# # 4 -0.1050765
# # 5 -0.4904108
# # 6 -0.7923543


## Rank-Inverse Based Normaliation on Pre-corrected Protein Phenotypes
library(bestNormalize)
for(i in colnames(phenotypes_residualised)[57:4291]){ 
  phenotypes_residualised[,i]<- orderNorm(phenotypes_residualised[,i])$x.t
}

## Scale somalogic data 
phenotypes_residualised[,57:4291] <- apply(phenotypes_residualised[,57:4291], 2, scale)

## Replace naming of the proteins with the unformatted SeqIds as saved above 
names(phenotypes_residualised[57:4291]) <- names_saved

# # Check 
# colMeans(phe[3:4237])
# apply(phe, 2, sd, na.rm = T)

### SAVE PHENO WITH EGFR + pQTLs regressed out 

# ## Save out somalogic file which has all covariates 
# write.csv(phenotypes_residualised, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_all_covariates_778_eGFR_inlcuded_pQTLs_regressed.csv", row.names = F)

## Save out somalogic file whih has no covariates 
phe <- phenotypes_residualised[,c(1,23,57:4291)]	
colnames(phe)[3:4237] <- names_saved
write.csv(phe, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778_eGFR_removed_pQTLs_regressed.csv", row.names = F)



# ### Do a correlation check between SMPD1 and HEXB

# prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778_eGFR_inlcuded_pQTLs_regressed.csv", check.names = F)

# list <- c("10818-36", "15470-11")

# prots <- prot[,which(colnames(prot) %in% list)]
# names(prots) <-  c("SMPD1", "HEXB")

# cor.test(prots$SMPD1, prots$HEXB)

# # t = 34.817, df = 776, p-value < 2.2e-16
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.7518092 0.8068412
# # sample estimates:
# #       cor
# # 0.7808354







