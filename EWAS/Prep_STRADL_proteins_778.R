####################################################################################
####################################################################################
############################## STRADL markers ######################################
####################################################################################
####################################################################################

# Prepping using Robs script as basis to keep parameters as consistent as possible
# He based it on Riccardo's initial processing of the STRADL data 

# We will need to update when the revised STRADL data comes through 

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

####################################################################################

## Context on where were at so far in the script 

####################################################################################

# This soma_demo1 file is my total group of people possible with protein data - 1061, with all the covariates apart from Set in it 
# For the methylation data, we will need to add set in when we regress proteins onto the covariates 

# there is a set file loaded here, which doesnt appear in the code previously taken from Rob
# im guessing its the set file that rob has saved in his folder with info for the 778 people 
# im going to load it in now 

set <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/Set_info_778.csv")

# Create a second soma_demo file called soma_demo2 which will be just the 778 with methylation
# and make sure set is merged into this file for those individuals below

soma_demo2 = merge(soma_demo1, set, by = "id")
matcher = set$id
soma_demo2 = soma_demo2[match(matcher, soma_demo2$id), ]

# Okay, so this soma_demo2 file has all the proteins and covariates needed for the EWAS with the 778 people only 

# soma_demo1 = all
# soma_demo2 = just EWAS


# the next part will transform it based on IDs in the EWAS and GWAS supplied, so ill load these in here so that it works:

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

## Regress Proteins onto Covariates 
phenotypes_residualised <- phenotypes
phenotypes_residualised$eGFR[phenotypes_residualised$eGFR %in% NA] <- mean(phenotypes_residualised$eGFR, na.rm = T)
for(i in colnames(phenotypes_residualised)[57:4291]){ 
  phenotypes_residualised[,i]<- lm(phenotypes_residualised[,i] ~ age + factor(sex) + factor(study_site) + factor(lag_group) + 
                                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, 
                                    na.action = na.exclude, data = phenotypes_residualised)$residuals
}


# ## Regress Proteins onto Covariates - with eGFR included 
# phenotypes_residualised <- phenotypes
# phenotypes_residualised$eGFR[phenotypes_residualised$eGFR %in% NA] <- mean(phenotypes_residualised$eGFR, na.rm = T)
# for(i in colnames(phenotypes_residualised)[57:4291]){ 
#   phenotypes_residualised[,i]<- lm(phenotypes_residualised[,i] ~ age + factor(sex) + factor(study_site) + factor(lag_group) + eGFR + 
#                                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, 
#                                     na.action = na.exclude, data = phenotypes_residualised)$residuals
# }

## Rank-Inverse Based Normaliation on Pre-corrected Protein Phenotypes
library(bestNormalize)
for(i in colnames(phenotypes_residualised)[57:4291]){ 
  phenotypes_residualised[,i]<- orderNorm(phenotypes_residualised[,i])$x.t
}

# rob <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/Phenotype_File_EWAS_778.csv")

# My file is now the same as Rob's file at this stage - good check for reproducibility 

## Scale somalogic data 
phenotypes_residualised[,57:4291] <- apply(phenotypes_residualised[,57:4291], 2, scale)

### SAVE PHNEO WITHOUR EGFR 

# ## Save out somalogic file which has all covariates 
write.csv(phenotypes_residualised, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_all_covariates_778_eGFR_removed.csv", row.names = F)

## Save out somalogic file whih has no covariates 
phe <- phenotypes_residualised[,c(1,23,57:4291)]	
write.csv(phe, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778_eGFR_removed.csv", row.names = F)

# ### SAVE PHENO WITH EGFR

# # ## Save out somalogic file which has all covariates 
# write.csv(phenotypes_residualised, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_all_covariates_778_eGFR_inlcuded.csv", row.names = F)

# ## Save out somalogic file whih has no covariates 
# phe <- phenotypes_residualised[,c(1,23,57:4291)]	
# write.csv(phe, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778_eGFR_inlcuded.csv", row.names = F)









# ## Log Transform 

# soma_demo2 <- soma_demo1[soma_demo1$id %in% id_ewas$x,] # gives 775 4302, so there are 775 which have their ID in the EWAS ID file 

# #soma_demo3 <- soma_demo2[soma_demo2$id %in% id_ewas$x,] # changed it to this? so that it keeps set in there for the regression below 

# ids = id_ewas$x # select IDs for EWAS in 778 

# soma_demo2 <- soma_demo2[match(ids,soma_demo2$id), ] 

# # There are 3 people with NA that are not in the EWAS IDs 

# phenotypes <- soma_demo2

# for(i in colnames(phenotypes)[57:4291]){ 
#   phenotypes[,i]<- log(phenotypes[,i])
# }

# dim(phenotypes) #  778 4303 - set is added still as 4303 column 



# ## Regress Proteins onto Covariates 

# phenotypes_residualised <- phenotypes
# phenotypes_residualised$eGFR[phenotypes_residualised$eGFR %in% NA] <- mean(phenotypes_residualised$eGFR, na.rm = T)

# for(i in colnames(phenotypes_residualised)[57:4291]){ 
#   phenotypes_residualised[,i]<- lm(phenotypes_residualised[,i] ~ age + factor(sex) + factor(study_site) + factor(Set) + factor(lag_group) + eGFR +
#                                      PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, 
#                                     na.action = na.exclude, data = phenotypes_residualised)$residuals
# }




# ## Rank-Inverse Based Normaliation on Pre-corrected Protein Phenotypes

# library(bestNormalize)

# for(i in colnames(phenotypes_residualised)[57:4291]){ 
#   phenotypes_residualised[,i]<- orderNorm(phenotypes_residualised[,i])$x.t
# }


# phenotypes_final <- phenotypes_residualised[,c(1,57:4291)]

# # id_ewas <- read.csv("IDs_778.csv")
# # id_gwas <- read.csv("IDs_1064.csv")



# # # phenotypes_ewas <- phenotypes_final[phenotypes_final$id %in% id_ewas$x,]
# # # ids = id_ewas$x 
# # # phenotypes_ewas <- phenotypes_ewas[match(ids,phenotypes_ewas$id), ]
# # # write.csv(phenotypes_final, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/Phenotype_File_EWAS.csv",row.names = F)



# # # phenotypes_gwas <- phenotypes_final[phenotypes_final$id %in% id_gwas$ID,]
# # # ids = id_gwas$ID 
# # # phenotypes_gwas <- phenotypes_gwas[match(ids,phenotypes_gwas$id), ]
# # # write.csv(phenotypes_gwas, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/Phenotype_File_GWAS_1064.csv",row.names= F)





# # ### This bit runs fine to get the covs for the EWAS (thanks rob)

# ## Prepare Covariates File 

# nrow(stradl)
# id <- read.csv("IDs_778.csv")
# covariates <- stradl[(stradl$GS_id %in% id$x),]
# ids = id$x
# covariates <- covariates[match(ids, covariates$GS_id),]
# covariates <- covariates[,c(1,7:13)]
# covariates[,2:8] <- apply(covariates[,2:8], 2, scale)
# covariates$GS_id <- NULL
# covariates <- covariates[1:778, -6]


# # This is the list of proteins Rob has been working with as a reference - make sure mine do not overlap 
# interesting_proteins <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Prep/List of interesting proteins.csv")

