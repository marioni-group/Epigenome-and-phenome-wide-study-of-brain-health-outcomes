################################################################################################################################
################################################################################################################################
############################## Imaging in the STRADL protein dataset ###########################################################
################################################################################################################################
################################################################################################################################

# Run proteome-wide studies of imaging phenotypes

# 070621 - models set off 

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
# table(is.na(proteins))

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

# Join phenotype info to protein dataset 
prot <- left_join(prot, demo, by = "SampleId")

# Join in the GS id linkage so that further phenotypes can be joined in 
IDs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/ST id linkage.csv")
names(IDs)[2] <- "GS_id"
names(IDs)[1] <- "SampleId"
IDs <- IDs[c(1,2)]

prot <- left_join(prot, IDs, by = "SampleId")

# Join in APOE and cognitive phenotypes 
APOE <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/APOE_snps/apoe_GS.csv")
names(APOE)[1] <- "GS_id"

prot <- left_join(prot, APOE, by = "GS_id")

# Cognitive scores - read in and join to dataset 
cog <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/st_cognitive.csv")
names(cog)[1] <- "SampleId"

prot <- left_join(prot, cog, by = "SampleId")

# Now join in the processed cognitive data from the code daniel shared with me - composite gf and g scores 
comp <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/cog1_240321.csv")
names(comp)[1] <- "SampleId"
comp <- comp[c(1,5,6,7)]

prot <- left_join(prot, comp, by = "SampleId")

# Add depression status into the dataset (read in file generated in depression covariate check, that is prepped with combined case/controls) 
dep <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/SCID_joint_to_GP_hospital_combined.csv")
dep <- dep[c(1,6)]
dep$combined <- as.character(dep$combined)
names(dep)[1] <- "SampleId"
prot <- left_join(prot, dep, by = "SampleId")
table(is.na(prot$combined))
prot$combined[is.na(prot$combined)] = 0


# Add in eGFR as covariate 
eGFR <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/eGFR_in_1065.csv")
names(eGFR)[1] <- "GS_id"
prot <- left_join(prot, eGFR, by = "GS_id")
prot$eGFR[prot$eGFR %in% NA] <- mean(prot$eGFR, na.rm = T)

# Add BMI as a covariate 
demo_table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/protAA/demo_table.csv")
BMI <- demo_table[c(14,54)]
prot = left_join(prot, BMI, by = "GS_id")
library(imputeTS)
prot$BMI = na_mean(prot$BMI) 

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
prot <- left_join(prot, ids, by = "SampleId")

# Check how many have each of the apoe phenotypes 
test <- prot
test$apoe %>% unique() #  "e3e4" "e3e3" "e2e3" NA     "e4e4" "e2e4" "e2e2"
table(is.na(test$apoe)) # 15 missing NA values 
outcomes <- test$apoe %>% as.data.frame() 
names(outcomes)[1] <- "X"
count(outcomes, X)

prot <- prot %>% mutate(APOE = case_when(
  apoe == "e4e4" ~ 2,
  apoe == "e3e4" ~ 2,
  apoe == "e3e3" ~ 1,
  apoe == "e2e2" ~ 0,
  apoe == "e2e3" ~ 0))

test <- prot
outcomes <- test$APOE %>% as.data.frame() 
names(outcomes)[1] <- "X"
count(outcomes, X)

# # # Write out to use in heatmaps 
# write.csv(prot, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/file_for_heatmaps_cog.csv", row.names = F)

#################################################################################################################

### ADD IN IMAGING DATA 

# Read in variables and get a sense of numbers 
brain_age <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_Brain_Measures_BrainAge.csv")
age <- brain_age[c(1,5)]
names(age)[1] <- "ID"

vol <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_Brain_Measures_Volumes.csv")
# of 1068 individuals with information on volumes, there are 36 NA values in the dataset 

# WMHV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_Brain_Measures_WMHV_Aberdeen.csv")
# # There are 474 individuals with WMHV

faz <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_Brain_Measures_Fazekas.csv")
# 1063 have faz scores for WMH, with one missing NA value 

WM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Imaging_data/STRADL_Brain_Measures_gFA-gMD.csv")


### Join
vol <- full_join(vol, faz, by = "ID", all = TRUE)
vol <- full_join(vol, WM, by = "ID", all = TRUE)
vol <- full_join(vol, age, by = "ID", all = TRUE)

### Join the imaging data into the protein dataset for regressions 
# We want imaging data for all of the individuals with protein data, so will leftjoin 
names(vol)[1] <- "SampleId"
prot <- left_join(prot, vol, by = "SampleId")

# Work out roughly who has both imaging and protein data available 

table(is.na(prot$WBV_With_Ventricles)) # 944 
table(is.na(prot$Fazekas_Score_Total)) # 944

# write for case study plots for one individual 
write.csv(prot, "/Cluster_Filespace/Marioni_Group/Danni/Case_study/STRADL_010621.csv", row.names = F)

# write for protAA comparisons
im <- prot[c(7,4299:4314)]
write.csv(im, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/protAA/imaging.csv", row.names = F)

# Calculate brain age acceleration score 
prot$brain_accel <- resid(lm(Brain_age ~ st_age, na.action = na.exclude, data = prot))

#+ prot$eGFR + prot$BMI + as.factor(prot$ever_smoke)

# Save out prot file for plots of brain imaging & cog variables 
write.csv(prot, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/prot_file_150621.csv", row.names = F)


# Check trends for gFA and gMD with age 

cor.test(prot$gFA, prot$st_age)
cor.test(prot$gMD, prot$st_age)




#####################

# Run by phenotype 

length <- 4235
phenotype <- c("brain_accel")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("brain_accel")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + prot$eGFR + prot$BMI + as.factor(prot$ever_smoke) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result_with_covs.csv"))
}



length <- 4235
phenotype <- c("gFA")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result.csv"))
}





length <- 4235
phenotype <- c("gFA")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + prot$eGFR + prot$BMI + as.factor(prot$ever_smoke) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result_with_covs.csv"))
}





length <- 4235
phenotype <- c("gMD")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result.csv"))
}




length <- 4235
phenotype <- c("gMD")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + prot$eGFR + prot$BMI + as.factor(prot$ever_smoke) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result_with_covs.csv"))
}



length <- 4235
phenotype <- c("Fazekas_Score_Total")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result.csv"))
}




length <- 4235
phenotype <- c("Fazekas_Score_Total")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site) + prot$eGFR + prot$BMI + as.factor(prot$ever_smoke) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result_with_covs.csv"))
}




length <- 4235
phenotype <- c("Cerebrum_WM_Volume")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result.csv"))
}




length <- 4235
phenotype <- c("Cerebrum_WM_Volume")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + prot$eGFR + prot$BMI + as.factor(prot$ever_smoke) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result_with_covs.csv"))
}




length <- 4235
phenotype <- c("WBV_With_Ventricles")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result.csv"))
}




length <- 4235
phenotype <- c("WBV_With_Ventricles")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + prot$eGFR + prot$BMI + as.factor(prot$ever_smoke) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result_with_covs.csv"))
}




length <- 4235
phenotype <- c("WBV_No_Ventricles")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result.csv"))
}




length <- 4235
phenotype <- c("WBV_No_Ventricles")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + prot$eGFR + prot$BMI + as.factor(prot$ever_smoke) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result_with_covs.csv"))
}




length <- 4235
phenotype <- c("Global_GM_Volume")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result.csv"))
}




length <- 4235
phenotype <- c("Global_GM_Volume")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + as.factor(prot$Site)*prot$Estimated_ICV + prot$eGFR + prot$BMI + as.factor(prot$ever_smoke) + as.factor(prot$Edited) + as.factor(prot$Batch) + (1|prot$id), na.action = na.exclude, data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- pheno
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]
    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/", pheno, "_result_with_covs.csv"))
}




##############################################################################################################

### PROCESS THE RESULTS FOR IMAGING MODELS (with and without covs)

##############################################################################################################

library(tidyverse)
library(readxl)

# Lets try one as an example and then we can work it up on loop 
location <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results/"

location_output <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/"
loop = list.files(paste0(location), ".")

pheno_count <- data.frame(Phenotype = 1:18, Associations = 1:18, Protein_Gene = 1:18)
pheno_count2 <- data.frame(Phenotype = 1:18, Associations = 1:18, Protein_Gene = 1:18)
pheno_count3 <- data.frame(Phenotype = 1:18, Associations = 1:18, Protein_Gene = 1:18)
pheno_count4 <- data.frame(Phenotype = 1:18, Associations = 1:18, Protein_Gene = 1:18)

list_proteins <- list()
list_genes <- list()
list_proteins2 <- list()
list_genes2 <- list()

for(i in 1:length(loop)){

  ### DO FDR version
  j <- loop[i]
  title = gsub(".csv", "", j)
  title <- as.character(title)

  results <- read.csv(paste0(location, j))
  results <- results[-1]

  # manually calculate P values
  results$Pcalc <- pchisq((results$beta/results$SE)^2, 1, lower.tail=F)

  results <- results[order(results$Pcalc),] # order by P value 
  results$P_adjust <- p.adjust(results$Pcalc, method = "BH") # FDR adjustment 
  results$Status <- ifelse(results$P_adjust < 0.05, "pass", "fail") # Add identifier for those that pass/fail significance in each 
  results <- results[-6]

  anno <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
  anno <- as.data.frame(anno)
  anno <- anno[c(1,18,13,4,2)] # annotations 

  results$SeqId <- as.character(results$SeqId)
  anno$SeqId <- as.character(anno$SeqId)

  joint <- left_join(results, anno, by = "SeqId") # join in annotations
  joint <- joint[c(1,9,2:8,10:12)] # order results 

  write.csv(joint, paste0(location_output, title, "_annotated_FDR.csv"), row.names = F)

  # Output the summaries to seqID and gene name tables 
  sig <- which(results$P_adjust < 0.05) # count those that pass 
  count <- length(sig)[1]
  pheno_count[i,1] <- title
  pheno_count[i,2] <- count
  pheno_count3[i,1] <- title
  pheno_count3[i,2] <- count

  sign <- joint[sig,] # filter to those with just significant P adjusted values 
  names(sign)[2] <- "Gene_name"
  genes <- sign$Gene_name
  proteins <- sign$SeqId
  str <- str_c(proteins, collapse = " ")
  str2 <- str_c(genes, collapse = " ")
  
  pheno_count[i,3] <- str
  pheno_count3[i,3] <- str2

  list_proteins[i] <- list(proteins)
  list_genes[i] <- list(genes)
  print(dim(results))

  ############################################################################

  ### DO Bonferroni version
  # 1.180638e-05
  # 0.0000180638
  j <- loop[i]
  title = gsub(".csv", "", j)
  title <- as.character(title)

  results <- read.csv(paste0(location, j))
  results <- results[-1]

  # manually calculate P values
  results$Pcalc <- pchisq((results$beta/results$SE)^2, 1, lower.tail=F)

  results <- results[order(results$Pcalc),] # order by P value 
  results$Status <- ifelse(results$Pcalc < 0.0000180638, "pass", "fail") 
  results <- results[-6]

  anno <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
  anno <- as.data.frame(anno)
  anno <- anno[c(1,18,13,4,2)] # annotations 

  results$SeqId <- as.character(results$SeqId)
  anno$SeqId <- as.character(anno$SeqId)

  joint <- left_join(results, anno, by = "SeqId") # join in annotations
  joint <- joint[c(1,9,2:8,10:11)] # order results 

  write.csv(joint, paste0(location_output, title, "_annotated_bon.csv"), row.names = F)

  # Output the summaries to seqID and gene name tables 
  sig <- which(results$Pcalc < 0.0000180638) # count those that pass 
  count <- length(sig)[1]
  pheno_count2[i,1] <- title
  pheno_count2[i,2] <- count
  pheno_count4[i,1] <- title
  pheno_count4[i,2] <- count

  sign <- joint[sig,] # filter to those with just significant P adjusted values 
  names(sign)[9] <- "Gene_name"
  genes <- sign$Gene_name
  proteins <- sign$SeqId
  str <- str_c(proteins, collapse = ", ")
  str2 <- str_c(genes, collapse = " ")
  
  pheno_count2[i,3] <- str
  pheno_count4[i,3] <- str2

  list_proteins2[i] <- list(proteins)
  list_genes2[i] <- list(genes)

  print(dim(results))
}

# Write out summaries of significant proteins/genes for each marker in tables 
write.csv(pheno_count, paste0(location_output, "IM_FDR_proteins.csv"), row.names = F)

write.csv(pheno_count2, paste0(location_output, "IM_BON_proteins.csv"), row.names = F)

write.csv(pheno_count3, paste0(location_output, "IM_FDR_genes.csv"), row.names = F)

write.csv(pheno_count4, paste0(location_output, "IM_BON_genes.csv"), row.names = F)


######################################################################################

### COLLATE TO A SUMMARY TABLE FOR SUPPL RESULTS 

# Brain age accel 
# Global GM vol
# gFA
# gMD
# Cerebrum total vol 
# Fazekas score total 
# WBV with ventricles 

# We want 1-4235 proteins on the left, then each table with phenotype joined in on the right in layers in this order 

brain <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/brain_accel_result_annotated_FDR.csv")
brain <- brain[c(1,2,10,3:8)]
brain$phenotype <- "Brain Age Acceleration"
names(brain) <- c("SeqId", "Gene of Protein", "UniProt Full Name", "Phenotype 1", "n", "Beta", "SE", "P", "FDR P")

GGM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/Global_GM_Volume_result_annotated_FDR.csv")
GGM <- GGM[c(1,3:8)]
GGM$phenotype <- "Global Grey Matter Volume"
names(GGM) <- c("SeqId", "Phenotype 2", "n", "Beta", "SE", "P", "FDR P")
GGM <- GGM[match(brain$SeqId, GGM$SeqId),]
GGM <- GGM[-1]

GFA <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/gFA_result_annotated_FDR.csv")
GFA <- GFA[c(1,3:8)]
GFA$phenotype <- "General Fractional Anisotropy"
names(GFA) <- c("SeqId", "Phenotype 3", "n", "Beta", "SE", "P", "FDR P")
GFA <- GFA[match(brain$SeqId, GFA$SeqId),]
GFA <- GFA[-1]

GMD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/gMD_result_annotated_FDR.csv")
GMD <- GMD[c(1,3:8)]
GMD$phenotype <- "General Mean Diffusivity"
names(GMD) <- c("SeqId", "Phenotype 4", "n", "Beta", "SE", "P", "FDR P")
GMD <- GMD[match(brain$SeqId, GMD$SeqId),]
GMD <- GMD[-1]

# CER <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/Cerebrum_WM_Volume_result_annotated_FDR.csv")
# CER <- CER[c(1,3:8)]
# CER$phenotype <- "Cerebrum White Matter Volume"
# names(CER) <- c("SeqId", "Phenotype 5", "n", "Beta", "SE", "P", "FDR P")
# CER <- CER[match(brain$SeqId, CER$SeqId),]
# CER <- CER[-1]

# FAZ <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/Fazekas_Score_Total_result_annotated_FDR.csv")
# FAZ <- FAZ[c(1,3:8)]
# FAZ$phenotype <- "Fazkeas White Matter Score"
# names(FAZ) <- c("SeqId", "Phenotype 5", "n", "Beta", "SE", "P", "FDR P")
# FAZ <- FAZ[match(brain$SeqId, FAZ$SeqId),]
# FAZ <- FAZ[-1]

WM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_WMHV_110821/WMHV_annotated_FDR.csv")
WM <- WM[c(1,3:8)]
WM$phenotype <- "White Matter Hyperintensity Volume"
names(WM) <- c("SeqId", "Phenotype 5", "n", "Beta", "SE", "P", "FDR P")
WM <- WM[match(brain$SeqId, WM$SeqId),]
WM <- WM[-1]

WBV <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/WBV_With_Ventricles_result_annotated_FDR.csv")
WBV <- WBV[c(1,3:8)]
WBV$phenotype <- "Whole Brain Volume"
names(WBV) <- c("SeqId", "Phenotype 6", "n", "Beta", "SE", "P", "FDR P")
WBV <- WBV[match(brain$SeqId, WBV$SeqId),]
WBV <- WBV[-1]

# Join together MATCHED files 

join <- cbind(brain, GGM)
join <- cbind(join, GFA)
join <- cbind(join, GMD)
# join <- cbind(join, CER)
join <- cbind(join, WM)
join <- cbind(join, WBV)

# Write out suppl table file 
write.csv(join, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_070621_all/imaging_results_annotated/no_eGFR_IM_joint_suppl_table.csv", row.names = F)
write.csv(join, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/no_eGFR_IM_joint_suppl_table.csv", row.names = F)


















