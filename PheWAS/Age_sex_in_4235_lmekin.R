################################################################################################################################
################################################################################################################################
############################## Age and sex in the STRADL protein dataset #######################################################
################################################################################################################################
################################################################################################################################

# Aim - replicate previous age/sex/age*sex findings by Lehallier et al 
# By doing this, we will identify the age/sex proteomes for our own sample in addition to the replication assessment

setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Aging_sex_profiling/")

library(tidyverse)
library(readxl)

################################################################################################################################

### STRADL DATA 

################################################################################################################################

# Protein data (raw, without processing)
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/GS+soma+QC+normalized.csv")

# Annotation linker file for SeqIds
link <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/annotation.csv", check.names = F)

# Update naming so were working in SeqIds 
names <- colnames(link)
names(prot)[c(33:4267)] <- names

# Load target file 
target <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/STRADL_DNAm_target_REM_17April2020.txt")

# > dim(target)
# [1] 847   6

# Check for missing protein data 
proteins <- prot[c(13,33:4267)]

table(is.na(proteins))

# > table(is.na(proteins))
#   FALSE
# 4511340

################################################################################################################################

### TRANSFORM PROTEIN DATA AND JOIN FOR REGRESSIONS

################################################################################################################################

# For now, I will just log transform and scale as Rob did with protAA calc 

## Log Transform 
for(i in colnames(prot)[33:4267]){ 
  prot[,i]<- log(prot[,i])
}

## Rank-Inverse Based Normaliation
library(bestNormalize)
for(i in colnames(prot)[33:4267]){ 
  prot[,i]<- orderNorm(prot[,i])$x.t
}

################################################################################################################################

### LINEAR REGRESSIONS

################################################################################################################################

# Load demographics for all people in STRADL 
demo <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/demographicsV2.csv")
names(demo)[1] <- "SampleId"

# > dim(demo)
# [1] 1069   11

# Join phenotype info to protein dataset 
prot <- left_join(prot, demo, by = "SampleId")

# Add depression status into the dataset (read in file generated in depression covariate check, that is prepped with combined case/controls) 
dep <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/SCID_joint_to_GP_hospital_combined.csv")
dep <- dep[c(1,6)]
dep$combined <- as.character(dep$combined)
names(dep)[1] <- "SampleId"
prot <- left_join(prot, dep, by = "SampleId")
table(is.na(prot$combined))

prot$combined[is.na(prot$combined)] = 0

## Set marker names to loop through in models as proteins (x)
markers <- prot[c(33:4267)]
marker_names <- colnames(markers)

### LOAD IN COXME REQUIREMENTS 

library(survival)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

## Code that is already processed and read in as the ped file below
# ped = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree.csv")
# ped$father <- as.numeric(ped$father)
# ped$mother <- as.numeric(ped$mother)
# ped$father[ped$father==0] <- NA
# ped$mother[ped$mother==0] <- NA
# table(ped$sex)
# ped$sex <- as.numeric(ped$sex)
# ped$sex[ped$sex==2] <- 0
# ped$sex <- ped$sex+1

# Read in the prepped file to cluster 
ped <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree_formatted.csv")

# Create kinship matrix for GS
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid)) # Pedigree list with 26 total subjects in 5 families
kin_model <- kinship(kin) 


# Function to Extract Lmekin Results	
extract_coxme_table <- function (mod){
  #beta <- mod$coefficients #$fixed is not needed
  beta <- fixef(mod)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- beta/se
  p<- 1 - pchisq((beta/se)^2, 1)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

# Join in the GS id information into the protein dataset for models
ids <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/ST id linkage.csv")
names(ids)[1] <- "SampleId"
ids <- ids[c(1,2)]
library(tidyverse)
prot <- left_join(prot, ids, by = "SampleId")

#####################################################################################################################################

### BATCHES

length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 1:500){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b1.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b1.csv", row.names = F)


length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 501:1000){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b2.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b2.csv", row.names = F)


length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 1001:1500){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b3.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b3.csv", row.names = F)


length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 1501:2000){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# batch 4
# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b4.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b4.csv", row.names = F)



length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 2001:2500){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b5.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b5.csv", row.names = F)


length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 2501:3000){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b6.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b6.csv", row.names = F)



length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 3001:3500){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b7.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b7.csv", row.names = F)



length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 3501:4000){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b8.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b8.csv", row.names = F)


length <- 4235
results <- data.frame(SeqId = 1:length, Age_beta = 1:length, Age_SE = 1:length, Age_P = 1:length)
results2 <- data.frame(SeqId = 1:length, Sex_beta = 1:length, Sex_SE = 1:length, Sex_P = 1:length)

for(i in 4001:4235){
  prot_name <- as.character(colnames(markers[i]))
  
  mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

  print(i)
  print(prot_name)
  results[i,1] <- prot_name
  results[i,2] <- extract_coxme_table(mod)[2,1]
  results[i,3] <- extract_coxme_table(mod)[2,2]
  results[i,4] <- extract_coxme_table(mod)[2,4]

  results2[i,1] <- prot_name
  results2[i,2] <- extract_coxme_table(mod)[3,1]
  results2[i,3] <- extract_coxme_table(mod)[3,2]
  results2[i,4] <- extract_coxme_table(mod)[3,4]
}

# save a copy of results tables 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b9.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b9.csv", row.names = F)


#########################################################################################################################

# Sort the results for each model and combine the batches

# Age 

e1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b1.csv")
e2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b2.csv")
e3 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b3.csv")
e4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b4.csv")
e5 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b5.csv")
e6 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b6.csv")
e7 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b7.csv")
e8 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b8.csv")
e9 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/age_b9.csv")

e1 <- e1[1:500,]
e2 <- e2[501:1000,]
e3 <- e3[1001:1500,]
e4 <- e4[1501:2000,]
e5 <- e5[2001:2500,]
e6 <- e6[2501:3000,]
e7 <- e7[3001:3500,]
e8 <- e8[3501:4000,]
e9 <- e9[4001:4235,]

e <- rbind(e1, e2)
e <- rbind(e, e3)
e <- rbind(e, e4)
e <- rbind(e, e5)
e <- rbind(e, e6)
e <- rbind(e, e7)
e <- rbind(e, e8)
e <- rbind(e, e9)

results <- e

# Sex 

e1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b1.csv")
e2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b2.csv")
e3 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b3.csv")
e4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b4.csv")
e5 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b5.csv")
e6 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b6.csv")
e7 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b7.csv")
e8 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b8.csv")
e9 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/individual_results/sex_b9.csv")



e1 <- e1[1:500,]
e2 <- e2[501:1000,]
e3 <- e3[1001:1500,]
e4 <- e4[1501:2000,]
e5 <- e5[2001:2500,]
e6 <- e6[2501:3000,]
e7 <- e7[3001:3500,]
e8 <- e8[3501:4000,]
e9 <- e9[4001:4235,]

e <- rbind(e1, e2)
e <- rbind(e, e3)
e <- rbind(e, e4)
e <- rbind(e, e5)
e <- rbind(e, e6)
e <- rbind(e, e7)
e <- rbind(e, e8)
e <- rbind(e, e9)

results2 <- e

# Rank
results <- results[order(results$Age_P),]
results2 <- results2[order(results2$Sex_P),]

# Merge in protein info 
anno <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
anno <- as.data.frame(anno)
anno <- anno[c(1,18,13,4,2)]

results <- left_join(results, anno, by = "SeqId")

results2 <- left_join(results2, anno, by = "SeqId")

# Manually calculate p vals 
results$Age_Pcalc <- pchisq((results$Age_beta/results$Age_SE)^2, 1, lower.tail=F)
results2$Sex_Pcalc <- pchisq((results2$Sex_beta/results2$Sex_SE)^2, 1, lower.tail=F)

# Plot p vals extracted vs calculated 
plot(results$Age_P, results$Age_Pcalc) # these look good

# Do p value adjustment 
results$Age_P_adjust <- p.adjust(results$Age_Pcalc, method = "BH")
results2$Sex_P_adjust <- p.adjust(results2$Sex_Pcalc, method = "BH")

# Add identifier for those that pass/fail significance in each 
results$Age_Status <- ifelse(results$Age_P_adjust < 0.05, "pass", "fail")
results2$Sex_Status <- ifelse(results2$Sex_P_adjust < 0.05, "pass", "fail")

# Count those that pass
sig <- which(results$Age_P_adjust < 0.05)
sig <- results[sig,]
dim(sig) # 546 protein traits identified - with lmekin it is 583 - now 798 - now 801 

# Count those that pass
sig2 <- which(results2$Sex_P_adjust < 0.05)
sig2 <- results2[sig2,]
dim(sig2) # 607 protein traits identified - with lmekin it is 623 - now 814 - now 799 

# Order by new p vals calculated manually 
results <- results[order(results$Age_Pcalc),]
results2 <- results2[order(results2$Sex_Pcalc),]

# Save out results tables for further processing 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/age_results_file.csv", row.names = F)
write.csv(results2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/age_sex_dep_cov/sex_results_file.csv", row.names = F)

# Move these files over to datastore for the Fig 1 plotting and replication assessment 

###########################################################################

### REPLICATION FOR LEHALIER

# BLSA/GESTALT study - tanaka study: https://onlinelibrary.wiley.com/doi/full/10.1111/acel.12799
# INTERVAL/LonGenity study - lehalier study: https://www.nature.com/articles/s41591-019-0673-2#Sec8
# InCHIANTI study - https://elifesciences.org/articles/61073
# As lehalier is the largest to date, we will just use that one for the replication assessment 

### First - look at how many age assocs were also sex significant 

setwd("Y:/Danni/stradl_markers/pheWAS/03_lmekin_rank_inverse_280321/results_age_sex/")
age <- read.csv("pcalc_main_RI_age_results_file.csv")
sex <- read.csv("pcalc_main_RI_sex_results_file.csv")

library(tidyverse)
library(readxl)

# Merge age and sex results together
agesex <- left_join(age, sex, by = "SeqId")

# Work out how many age assocs were also sex assocs 
age_sig <- filter(agesex, Age_P_adjust < 0.05) # 583 age associations - now 798
agesexsig <- filter(age_sig, Sex_P_adjust < 0.05) # now 394

# of the 798, 394 were also associated with sex

### Now - COMPARE TO LEHALLIER RESULTS 
leh <- read_excel("Y:/Danni/stradl_markers/pheWAS/03_lmekin_rank_inverse_280321/results_age_sex/replication_assessment/INTERVAL_lehalier_results.xlsx")
leh <- as.data.frame(leh)

dim(leh) # 2925 - matches paper!

# # Sort out the naming situation for joining to our results - no longer needed thanks to linker used above 
library(stringr)
x = strsplit(leh$ID, ".", fixed = T)
head(x)
ids = lapply(x, function(x){paste0(x[2], "-", x[3])})
ids = unlist(ids)
leh$SeqId <- ids

# Now merge in the seqIds and make sure the conversion has been completed for all of them
anno <- read_excel("Y:/Protein_DNAm_Proxies/Manuscript_revision/Annotations_for_reference.xlsx")
anno <- as.data.frame(anno)
anno <- anno[c(1,18,13,4,2)]
anno$check_col <- anno$SeqId
anno <- anno[c(1,6,2:5)]
leh <- left_join(leh, anno, by = "SeqId")
write.csv(leh, "lehalier_renamed_joint_to_anno_using_extraction_method_chopped.csv", row.names = F) # check the conversions and manually edit any issues 

# # Read in the linker file to sort out the SeqIds 
# link <- read_excel("Y:/Danni/stradl_markers/pheWAS/03_lmekin_rank_inverse_280321/results_age_sex/replication_assessment/lehalier_linker_file.xlsx")
# link$SeqId <- as.character(link$SeqId)
# link <- as.data.frame(link)
# link$SeqId <- as.character(link$SeqId)
# link$SeqId <- gsub("\\.", "-", link$SeqId) # format the SeqIds to be consistent with ours 

# # Join the seqids into the lehalier human results file 
# leh <- left_join(leh, link, by = "ID")

# Read in the correct SeqIds file as edited to remove any issues with naming 
leh <- read_excel("lehalier_renamed_joint_to_anno_using_extraction_method_chopped_edited_excel_file.xlsx")
leh <- as.data.frame(leh)
length(unique(leh$SeqId)) # 2925 - all good

# Correct the NA values for 4 issue cells read in 
t <- na.omit(leh)
diff <-  which(!leh$SeqId %in% t$SeqId)
issue <- leh[c(2924:2925),]

# Add in the correct values from the table in place of the incorrect NA values 
issue[1,2] <- 1.97626258336499E-323
issue[1,3] <- 2.89028402817129E-320

issue[2,2] <- 9.88131291682493E-324
issue[2,3] <- 2.89028402817129E-320

# Join the right 4 cells back to main dataset 
leh <- rbind(t, issue) # 2925
leh <- leh[c(14,2:7)]
names(leh) <- c("SeqId", "leh_Age_P", "leh_Age_P_adj", "leh_Age_beta", "leh_Sex_P", "leh_Sex_P_adj", "leh_Sex_beta")
write.csv(leh, "lehalier_renamed_and_complete.csv", row.names = F)

joint <- agesex # Assign our results to a variable called 'joint' for age and sex 

# check to see how many proteins match across all tables 
length(which(leh$SeqId %in% joint$SeqId)) # 2438

### How many were the comparisons unavailable for?
length(which(!leh$SeqId %in% joint$SeqId)) # 487 from our file arent in leh file 

length(which(!joint$SeqId %in% leh$SeqId)) # 1800 from leh file arent in our file  

# Subset to common proteins across our results and their results files 
joint <- joint[which(joint$SeqId %in% leh$SeqId),]  # 2435  in leh file that are available as direct comparisons to ours 

j_sig <- joint %>% filter(Age_P_adjust < 0.05) # 476 possible matches for replication in ours 
dim(j_sig)

j_sig2 <- joint %>% filter(Sex_P_adjust < 0.05) # 450 possible matches for replication in ours 
dim(j_sig2)

# Do p values corrected version for replication
names(leh)[3] <- "adj_age"
names(leh)[6] <- "adj_sex"
leh$adj_sex <- as.numeric(leh$adj_sex)

# get lehanlier age sig assocs
age <- leh[which(leh$adj_age <= 0.05),] # 1,379 - matches their paper exactly  
# get lehallier sex sig assocs 
sex <- leh[which(leh$adj_sex <= 0.05),] # 1651 - matches their paper exactly 
# get the age and sex sig assocs 
agesex <- age %>% filter(adj_sex <= 0.05) # 895 - matches their paper exactly 

# See if we replicate the associations 
a_filter <- joint %>% filter(Age_P_adjust < 0.05) 

length(which(age$SeqId %in% a_filter$SeqId)) # 373 / 476 that match for age between sig results of both studies 

s_filter <- joint %>% filter(Sex_P_adjust < 0.05)

length(which(sex$SeqId %in% s_filter$SeqId)) # 381 / 450 that match for sex for sig results of both studies 

# Collate the data so we can check replications based on direction of effect for beta coeficient 
lehag <- leh[c(1,2,3,4)] # version with age results across the leh results 
lehage <- lehag[which(lehag$adj_age <= 0.05),] #1379
lehse <- leh[c(1,5,6,7)] # version with sex results across the leh results 
lehsex <- lehse[which(lehse$adj_sex <= 0.05),] #1651

# Now join our results to those which are present in lehalier study
j_sig <- j_sig[c(1:5,9,10)]
j_sig$SeqId <- as.character(j_sig$SeqId)
usage <- merge(j_sig, lehage, by = "SeqId", all = T) 
usage <- usage %>% filter(Age_P_adjust <= 0.05)
length(unique(usage$SeqId))
n_occur <- data.frame(table(usage$SeqId)) # its row , which is 2953-31 thats duplicated due to having 2 genes - but its same info so we can remove i think 
usage <- usage[-169,]
write.csv(usage, "age_replication_table.csv", row.names = F)

# Make sure all were available 
usage2 <- left_join(j_sig, lehag, by = "SeqId")
usage2 <- usage2[-37,]
write.csv(usage2, "age_replication_table_all.csv", row.names = F)

# Add a replication status for the age assocs 
usage2 <- usage2 %>% mutate(Replication = case_when(
  Age_P_adjust <= 0.05 & adj_age <= 0.05 ~ "1"))
check <- na.omit(usage2) # 372 - this is actually the true count for the replication match with the extra 2x removed 

# See if we can check the direction of effect has replicated too 
usage2 <- usage2 %>% mutate(Sign = case_when(
  Age_beta < 0 ~ "-",
  Age_beta > 0 ~ "+"))

usage2 <- usage2 %>% mutate(Sign2 = case_when(
  leh_Age_beta < 0 ~ "-",
  leh_Age_beta > 0 ~ "+"))

usage2$Replication_direction <- ifelse(usage2$Sign == usage2$Sign2, "1", NA)
check <- na.omit(usage2) # 340 of the 372 also replicated the direction of effect 

usage2 <- usage2 %>% mutate(Rep_both = case_when(
  Replication == "1" & Replication_direction == "1" ~ "1"))
write.csv(usage2, "age_replication_table_all_labelled_with_direction_of_effect.csv", row.names = F)

######################################
# Read back in table for plotting 
setwd("Y:/Danni/stradl_markers/pheWAS/03_lmekin_rank_inverse_280321/results_age_sex/")
library(tidyverse)
library(readxl)
library(ggplot2)
usage2 <- read.csv("age_replication_table_all_labelled_with_direction_of_effect.csv")

# Replace NA's with 0's 
usage2$Rep_both[is.na(usage2$Rep_both)] <- 0

# Plot effect sizes for the proteins where comparisons were available across the studies 
pdf("Plot_age_replication_BH_340_of_372_of_476_rep.pdf", width = 5, height = 5)
ggplot(usage2, aes(leh_Age_beta, Age_beta), scale="globalminmax") +
       geom_vline(xintercept = 0, linetype = 1) +
       geom_hline(yintercept = 0, linetype = 1) +
       geom_point() +
       theme_minimal() +
      xlab("Beta Lehallier et al") +
      ylab("Beta") +
       theme(plot.title = element_text(size = 14)) +
  # geom_abline(slope=1, intercept=0, linetype = 2) +
  geom_point(data=usage2[usage2$Rep_both =="1",], shape=21, fill="white", size=2)+
  geom_point(data=usage2[usage2$Rep_both !="1",], colour='red') 
dev.off()

# Convert to z-scores as effects across both studies 
# subtract mean and divide by sd

usage2$z_us <- (usage2$Age_beta - mean(usage2$Age_beta, na.rm = T)) / sd(usage2$Age_beta, na.rm = T)

usage2$z_leh <- (usage2$leh_Age_beta - mean(usage2$leh_Age_beta, na.rm = T)) / sd(usage2$leh_Age_beta, na.rm = T)

# Plot the z-transformed effect sizes where comparisons were available across the studies 
pdf("Plot_age_replication_BH_340_of_372_of_476_rep_z_scores.pdf", width = 5, height = 5)
ggplot(usage2, aes(z_leh, z_us), scale="globalminmax") +
       geom_vline(xintercept = 0, linetype = 1) +
       geom_hline(yintercept = 0, linetype = 1) +
       geom_point() +
       theme_minimal() +
      xlab("z-score effect Lehallier et al") +
      ylab("z-score effect") +
       theme(plot.title = element_text(size = 14)) +
   geom_abline(slope=1, intercept=0, linetype = 2) +
  geom_point(data=usage2[usage2$Rep_both =="1",], shape=21, fill="white", size=2)+
  geom_point(data=usage2[usage2$Rep_both !="1",], colour='red') + xlim(-10,10) + ylim(-10,10)
dev.off()

# save as variable for patchworking 
plot1 <- ggplot(usage2, aes(z_leh, z_us), scale="globalminmax") +
       geom_vline(xintercept = 0, linetype = 1) +
       geom_hline(yintercept = 0, linetype = 1) +
       geom_point() +
       theme_minimal() +
      xlab("Age z-score effect Lehallier et al") +
      ylab("Age z-score effect") +
       theme(plot.title = element_text(size = 14)) +
   geom_abline(slope=1, intercept=0, linetype = 2) +
  geom_point(data=usage2[usage2$Rep_both =="1",], shape=21, fill="white", size=1)+
  geom_point(data=usage2[usage2$Rep_both !="1",], colour='red') + xlim(-10,10) + ylim(-10,10)

# check on the odd outlier 
which(usage2$leh_Age_beta > 0.01) # row 6

#################################################
### REPEAT THE SAME FOR SEX REPLICATION
j_sig2 <- j_sig2[c(1,12:14,19,20)]
j_sig2$SeqId <- as.character(j_sig2$SeqId)
ussex <- merge(j_sig2, lehse, by = "SeqId", all = T) 
ussex2 <- ussex %>% filter(Sex_P_adjust <= 0.05) # 451
length(unique(ussex2$SeqId)) # 450
n_occur <- data.frame(table(ussex2$SeqId)) # its row , which is 2953-31 thats duplicated due to having 2 genes - but its same info so we can remove i think 
# usage <- usage[-169,]
# which(n_occur$Freq == "2")
# n <- n_occur[2295,]
write.csv(ussex2, "sex_replication_table.csv", row.names = F)

# Make sure all were available 
usage2 <- left_join(j_sig2, lehse, by = "SeqId")
usage2 <- usage2[-25,]
write.csv(usage2, "sex_replication_table_all.csv", row.names = F)

# Add a replication status for the age assocs 
usage2 <- usage2 %>% mutate(Replication = case_when(
  Sex_P_adjust <= 0.05 & adj_sex <= 0.05 ~ "1"))
check <- na.omit(usage2) # 380 - this is actually the true count for the replication match with the extra 2x removed 

# See if we can check the direction of effect has replicated too 
usage2 <- usage2 %>% mutate(Sign = case_when(
  Sex_beta < 0 ~ "-",
  Sex_beta > 0 ~ "+"))

usage2 <- usage2 %>% mutate(Sign2 = case_when(
  leh_Sex_beta < 0 ~ "-",
  leh_Sex_beta > 0 ~ "+"))

usage2$Replication_direction <- ifelse(usage2$Sign == usage2$Sign2, "1", NA)
check <- na.omit(usage2) # 372 of the 380 also replicated the direction of effect 

usage2 <- usage2 %>% mutate(Rep_both = case_when(
  Replication == "1" & Replication_direction == "1" ~ "1"))
write.csv(usage2, "sex_replication_table_all_labelled_with_direction_of_effect.csv", row.names = F)


######################################
# Read back in table for plotting 
setwd("Y:/Danni/stradl_markers/pheWAS/03_lmekin_rank_inverse_280321/results_age_sex/")
library(tidyverse)
library(readxl)
library(ggplot2)
usage2 <- read.csv("sex_replication_table_all_labelled_with_direction_of_effect.csv")

# Replace NA's with 0's 
usage2$Rep_both[is.na(usage2$Rep_both)] <- 0

# Plot the effect sizes for the proteins where comparisons were available across the studies 
pdf("Plot_sex_replication_BH_372_of_380_of_450_rep.pdf", width = 5, height = 5)
ggplot(usage2, aes(leh_Sex_beta, Sex_beta), scale="globalminmax") +
       geom_vline(xintercept = 0, linetype = 1) +
       geom_hline(yintercept = 0, linetype = 1) +
       geom_point() +
       theme_minimal() +
      xlab("Beta Lehallier et al") +
      ylab("Beta") +
       theme(plot.title = element_text(size = 14)) +
  # geom_abline(slope=1, intercept=0, linetype = 2) +
  geom_point(data=usage2[usage2$Rep_both =="1",], shape=21, fill="white", size=2)+
  geom_point(data=usage2[usage2$Rep_both !="1",], colour='red')
dev.off()

# Convert to z-scores as effects across both studies 
# subtract mean and divide by sd

usage2$z_us <- (usage2$Sex_beta - mean(usage2$Sex_beta, na.rm = T)) / sd(usage2$Sex_beta, na.rm = T)

usage2$z_leh <- (usage2$leh_Sex_beta - mean(usage2$leh_Sex_beta, na.rm = T)) / sd(usage2$leh_Sex_beta, na.rm = T)

# Plot the z-transformed effect sizes where comparisons were available across the studies 
pdf("Plot_sex_replication_BH_372_of_380_of_450_rep_z_scores.pdf", width = 5, height = 5)
ggplot(usage2, aes(z_leh, z_us), scale="globalminmax") +
       geom_vline(xintercept = 0, linetype = 1) +
       geom_hline(yintercept = 0, linetype = 1) +
       geom_point() +
       theme_minimal() +
      xlab("z-score effect Lehallier et al") +
      ylab("z-score effect") +
       theme(plot.title = element_text(size = 14)) +
   geom_abline(slope=1, intercept=0, linetype = 2) +
  geom_point(data=usage2[usage2$Rep_both =="1",], shape=21, fill="white", size=2)+
  geom_point(data=usage2[usage2$Rep_both !="1",], colour='red') + xlim(-9,9) + ylim(-9,9)
dev.off()

# save as variable for plotting via patchwork 
plot2 <- ggplot(usage2, aes(z_leh, z_us), scale="globalminmax") +
       geom_vline(xintercept = 0, linetype = 1) +
       geom_hline(yintercept = 0, linetype = 1) +
       geom_point() +
       theme_minimal() +
      xlab("Sex z-score effect Lehallier et al") +
      ylab("Sex z-score effect") +
       theme(plot.title = element_text(size = 14)) +
   geom_abline(slope=1, intercept=0, linetype = 2) +
  geom_point(data=usage2[usage2$Rep_both =="1",], shape=21, fill="white", size=1)+
  geom_point(data=usage2[usage2$Rep_both !="1",], colour='red') + xlim(-9,9) + ylim(-9,9)

# Patch together age and sex repl plots 
library(patchwork)
plot <- plot1 + plot2
pdf("Y:/Danni/stradl_markers/Plots/SUPPL_replication_age_sex/plot_joint_age_sex.pdf", width = 12, height = 6)
plot
dev.off()

#########################################################################################################

### NOW COLLATE THE RESULTS TABLES FOR REPORTING IN SUPPLEMENTARY TABLE

# Input the age and sex results for all 4235 proteins
# Mark out which of these replicated with lehalier for association and direction of effect 
# Save a table with both age and sex proteomes and whether they replicated 

library(tidyverse)
library(readxl)

setwd("Y:/Danni/stradl_markers/pheWAS/03_lmekin_rank_inverse_280321/results_age_sex/")
age <- read.csv("pcalc_main_RI_age_results_file.csv")
sex <- read.csv("pcalc_main_RI_sex_results_file.csv")

# select just the info we need for the table 
age <- age[c(1,5,6,2,3,9,10)]
sex <- sex[c(1,5,6,2,3,9,10)]

# add col for those with FDR P < 0.05 in our study 
age$Age_FDR_pass <- ifelse(age$Age_P_adjust <= 0.05, "yes", "no")
check <- age %>% filter(age$Age_FDR_pass == "yes")
dim(check) # 798 

sex$Sex_FDR_pass <- ifelse(sex$Sex_P_adjust <= 0.05, "yes", "no")
check <- sex %>% filter(sex$Sex_FDR_pass == "yes")
dim(check) # 798 

## Read in the replication tables with direction of effect indications  
rep1 <- read.csv("age_replication_table_all_labelled_with_direction_of_effect.csv")
rep2 <- read.csv("sex_replication_table_all_labelled_with_direction_of_effect.csv")

check <- rep1 %>% filter(Rep_both == "1")
dim(check) # 340 

check <- rep2 %>% filter(Rep_both == "1")
dim(check) # 372

# Create a column signifying whether status = replicated vs unreplicated
rep1$Rep_both[is.na(rep1$Rep_both)] <- 0 
rep1$Age_replication <- ifelse(rep1$Rep_both == "1", "Replicated", "Unreplicated")

rep2$Rep_both[is.na(rep2$Rep_both)] <- 0 
rep2$Sex_replication <- ifelse(rep2$Rep_both == "1", "Replicated", "Unreplicated")

# Select just the seqId and replication assessment for joining 
rep1 <- rep1[c(1,16)]
rep2 <- rep2[c(1,15)]

# Join to main results file for age and sex proteomes 
age_join <- left_join(age, rep1, by = "SeqId")
sex_join <- left_join(sex, rep2, by = "SeqId")

# Set the NA values to "Comparison unavailable"
age_join$Age_replication[is.na(age_join$Age_replication)] <- "Comparison unavailable"
sex_join$Sex_replication[is.na(sex_join$Sex_replication)] <- "Comparison unavailable"

# Join these tables together 
join <- left_join(age_join, sex_join, by = "SeqId")
join <- join[c(1:9,12:17)]

# Rename and save 
names(join) <- c("SeqId", "Entrez Gene Name", "Uniprot Full Name",
  "Age Beta", "Age SE", "Age P", "Age FDR P", "Age FDR P < 0.05", "Age Replication",
  "Sex Beta", "Sex SE", "Sex P", "Sex FDR P", "Sex FDR P < 0.05", "Sex Replication")

write.csv(join, "suppl_table_age_sex.csv", row.names = F)


