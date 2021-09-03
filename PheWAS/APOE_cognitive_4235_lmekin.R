################################################################################################################################
################################################################################################################################
############################## Cog and APOE in the STRADL protein dataset ######################################################
################################################################################################################################
################################################################################################################################

# Run proteome-wide studies of APOE and cognitive - with scaling 

# setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Aging_sex_profiling/")

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

# Scale age data prior to running model 
# prot[,4269:4270] <- apply(prot[,4269:4270], 2, scale)

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
prot <- left_join(prot, ids, by = "SampleId")


# ### ADDED: a new cognitive decline measure to run against 

# # read in GS baseline cog prep 
# base <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/STRADL/cog_score_prep/cog1_081220.csv")
# colnames(base) <- paste("Base", colnames(base), sep = "_")
# names(base)[1] <- "GS_id"

# # Join to the STRADL dataset 
# prot <- left_join(prot, base, by = "GS_id")

# # Get difference in g score 
# prot$gdiff <- prot$g - prot$Base_g
# prot$gdiff2 <-  prot$Base_g - prot$g

# prot$diffp <- (prot$gdiff / prot$Base_g) * 100
# prot$diff2p <- (prot$gdiff2 / prot$Base_g) * 100

## Lets try to divvy up the cognitive scores 
# phenotype <- c("digit_symbol")
# phenotype <- c("verbal_total")
# phenotype <- c( "gf")
# phenotype <- c("g")
# phenotype <- c("LM")
# phenotype <- c("mr_correct")
# phenotype <- c("vocabulary")

length <- 4235
phenotype <- c("digit_symbol")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

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
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/cog/individual/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("verbal_total")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

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
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/cog/individual/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("gf")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

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
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/cog/individual/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("g")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

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
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/cog/individual/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("LM")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

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
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/cog/individual/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("mr_correct")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

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
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/cog/individual/", pheno, "_result.csv"))
}


length <- 4235
phenotype <- c("vocabulary")

for (j in 1:length(phenotype)){
  pheno <- phenotype[j]
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))
  
      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

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
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/cog/individual/", pheno, "_result.csv"))
}





##############################################################################################################

### RUN APOE MODELS 

##############################################################################################################

# Check how many have each of the apoe phenotypes 

test <- prot
test$apoe %>% unique() #  "e3e4" "e3e3" "e2e3" NA     "e4e4" "e2e4" "e2e2"
table(is.na(test$apoe)) # 15 missing NA values 

outcomes <- test$apoe %>% as.data.frame() 
names(outcomes)[1] <- "X"
count(outcomes, X)

#      X   n
# 1 e2e2   5
# 2 e2e3 121
# 3 e2e4  22
# 4 e3e3 633
# 5 e3e4 234
# 6 e4e4  35
# 7 <NA>  15

# Split aope into bins based on the presence of e4 and e2 risk alleles 

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

#    X   n
# 1  0 126
# 2  1 633
# 3  2 269
# 4 NA  37


# Write out to use in heatmaps 
write.csv(prot, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/file_for_heatmaps_cog.csv", row.names = F)

# Now we can run APOE status through the loops as numeric continuum 
length <- 4235
set <- c("set")
phenotype <- c("APOE")

for (j in 1:length(set)){
  pheno <- "APOE"
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1:1000){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))

      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- paste0(pheno)
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]

    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/APOE/individual/e2_b1.csv"), row.names = F)
}


# Now we can run APOE status through the loops as numeric continuum 
length <- 4235
set <- c("set")
phenotype <- c("APOE")

for (j in 1:length(set)){
  pheno <- "APOE"
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 1001:2000){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))

      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- paste0(pheno)
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]

    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/APOE/individual/e2_b2.csv"), row.names = F)
}


# Now we can run APOE status through the loops as numeric continuum 
length <- 4235
set <- c("set")
phenotype <- c("APOE")

for (j in 1:length(set)){
  pheno <- "APOE"
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 2001:3000){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))

      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- paste0(pheno)
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]

    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/APOE/individual/e2_b3.csv"), row.names = F)
}


# Now we can run APOE status through the loops as numeric continuum 
length <- 4235
set <- c("set")
phenotype <- c("APOE")

for (j in 1:length(set)){
  pheno <- "APOE"
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 3001:4000){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))

      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- paste0(pheno)
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]

    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/APOE/individual/e2_b4.csv"), row.names = F)
}


# Now we can run APOE status through the loops as numeric continuum 
length <- 4235
set <- c("set")
phenotype <- c("APOE")

for (j in 1:length(set)){
  pheno <- "APOE"
  results <- data.frame(SeqId = 1:length, phenotype = 1:length, n = 1:length, beta = 1:length, SE = 1:length, P = 1:length)
  for(i in 4001:4235){
    tryCatch({ 
      prot_name <- as.character(colnames(markers[i]))

      mod <- lmekin(scale(prot[,prot_name]) ~ scale(prot$st_age) + as.factor(prot$sex) + scale(prot[,pheno]) + as.factor(prot$combined) + (1|prot$id), data = prot, varlist = kin_model*2)

      print(i)
      print(prot_name)
      print(pheno)

      results[i,1] <- prot_name
      results[i,2] <- paste0(pheno)
      results[i,3] <- mod$n
      results[i,4] <- extract_coxme_table(mod)[4,1]
      results[i,5] <- extract_coxme_table(mod)[4,2]
      results[i,6] <- extract_coxme_table(mod)[4,4]

    }, error = function(e) cat("skipped"))
  }
  write.csv(results, paste0("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/APOE/individual/e2_b5.csv"), row.names = F)
}


##############################################################################################################

### PROCESS APOE MODELS

##############################################################################################################

# Will need to join up the batches of APOE runs for e2 and e4 

e21 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/APOE/individual/e2_b1.csv")
e22 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/APOE/individual/e2_b2.csv")
e23 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/APOE/individual/e2_b3.csv")
e24 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/APOE/individual/e2_b4.csv")
e25 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/APOE/individual/e2_b5.csv")


# Cut the info you need and bind together (should have put NA's in but too late now)
e21 <- e21[1:1000,]
e22 <- e22[1001:2000,]
e23 <- e23[2001:3000,]
e24 <- e24[3001:4000,]
e25 <- e25[4001:4235,]

e2 <- rbind(e21, e22)
e2 <- rbind(e2, e23)
e2 <- rbind(e2, e24)
e2 <- rbind(e2, e25)


# Save out combined files into main repowhere processing code can run the annotations and pass/fail FDR 
write.csv(e2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled/APOE_JOINT.csv")

# Run the loop to generate annotated results files above with only these 2 APOE files in the folder and youll have them added to results
# APOE is processed in the covs APOE script 

# Do extraction using loop in the _covs run script for all phenos 

# Results are in this location
# /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_scaled_310521/results_scaled_annotated/

