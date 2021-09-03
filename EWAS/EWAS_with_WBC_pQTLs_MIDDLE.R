############################################################################################
############################################################################################
############################### EWA MOA 778 with genetic ORM ###############################
############################################################################################
############################################################################################

# Run 04/06/21

# EWAS with MOA on 778 individuals, using a genetic ORM to adjust for relatedness
# Accounting for covariates in a fully adjusted model:
# age
# Sex
# Wave
# Batch
# Depression status 
# 5x WBC imputed 

# Known pQTLs have also been regressed from the protein levels - residuals are used in 
# the models that therefore have known genetic effects (from Sun et al) adjusted for.

# This model will produce the version of EWAS that will be compared to lifestyle covariate version.

############################################################################################################

### Methylation file - 778

############################################################################################################

# Taken as per the WBC adjusted basic model without covariates of smoking or BMI dated 280521
# Only the protein levels have been altered (i.e. regressing known SNP effects from Sun et al).

## Read in methylation file used for BayesR and transpose so CpGs are columns and IDs are rows
library(data.table) 
w2 = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave2-STRADL-mvals.rds")
w3 = readRDS("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave3-STRADL-mvals.rds")
w2 = w2[which(row.names(w2) %in% row.names(w3)),]
w3 = w3[which(row.names(w3) %in% row.names(w2)),]

## Combine w2 and w3 
meth = cbind(w2,w3)

## Add in GS ids instead of DNAm_ids 
stradl_ids <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/STRADL_DNAm_target_REM_17April2020.txt",header =T)
id.778 <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/IDs_778.csv")

## Subset to the people in the ID file 
stradl_ids <- stradl_ids[which(stradl_ids$GS_id %in% id.778$x),]
ids = stradl_ids$DNAm_id
meth <- meth[,match(ids, colnames(meth))]
table(stradl_ids$DNAm_id == colnames(meth))
colnames(meth) <- stradl_ids$GS_id

## Prepare phenotypes file for covariates that we will use to pre-correct data 
stradl_ids <- stradl_ids[,c("GS_id", "age_stradl", "sex", "wave", "Stradl_id")]

# Add batch 
w2_batch <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave2_batch.csv") 
w2_batch$Batch[w2_batch$Batch == 1]<- "w1_1"
w2_batch$Batch[w2_batch$Batch == 2]<- "w1_2"
w2_batch$Batch[w2_batch$Batch == 3]<- "w1_3"
w2_batch$Batch[w2_batch$Batch == 4]<- "w1_4"
w2_batch$Batch[w2_batch$Batch == 5]<- "w1_5"
w2_batch$Batch[w2_batch$Batch == 6]<- "w1_6"
w3_batch <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/wave_3_original_samplesheet.csv") 
w2_batch <- w2_batch[,c("Sample_Name","Batch")]
names(w2_batch) <- c("GS_id", "Batch")
w3_batch <- w3_batch[,c("Sample_Name","Batch_all")]
names(w3_batch) <- c("GS_id", "Batch")
w3_batch$GS_id <- gsub("ST", "", w3_batch$GS_id)
batch = rbind(w2_batch, w3_batch)
batch = batch[-which(duplicated(batch$GS_id)),]

phenotypes= merge(stradl_ids, batch, by = "GS_id")

# Add WBC info into phenotypes for regressions 
w2cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_w2_cells.csv")
w3cells <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/stradl_w3_cells.csv")
merged <- rbind(w2cells, w3cells)
names(merged)[1] <- "GS_id"
WBC <- merged[c(1,3,7:12)]

phenotypes = merge(phenotypes, WBC, by = "GS_id")

# Add in BMI info into phenotypes for regressions 
demo_table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/protAA/demo_table.csv")
BMI <- demo_table[c(14,54)]

phenotypes = merge(phenotypes, BMI, by = "GS_id")
library(imputeTS)
phenotypes$BMI = na_mean(phenotypes$BMI) # 3 individuals with no BMI assigned mean imputed 

# Join in epismoker info for those in STRADL to linker file, then to cov file 
epi <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/stradl_epismoker.rds")
names(epi)[1] <- "DNAm_id"

phenotypes = merge(phenotypes, epi, by = "DNAm_id")

# Add depression status (based on SCID and incident case diagnoses pre STRADL sampling), as cov 
dep <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/Depression_case_controls_200521/prep_files/SCID_joint_to_GP_hospital_combined.csv")
dep <- dep[c(1,6)]
dep$combined <- as.character(dep$combined)
names(dep)[1] <- "Stradl_id"
library(tidyverse)
which(phenotypes$Stradl_id %in% dep$Stradl_id) # 2 people with no data, which we will need to assign 0 as they have no depression status of case recorded

phenotypes = left_join(phenotypes, dep, by = "Stradl_id")
phenotypes$combined[is.na(phenotypes$combined)] <- 0

# Get the set of finalised phneotypes for regressions 
phenotypes = phenotypes[c(2:5,7:16)]

# Match ID order between phenotype and meth files 
ids= colnames(meth)
phenotypes = phenotypes[match(ids, as.character(phenotypes$GS_id)),]
table(colnames(meth) == phenotypes$GS_id)

# Write out covariates to use in heatmaps 
# write.csv(phenotypes, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/corrplots/phenotypes_file_778_depression_added_240521.csv", row.names = F)

# Transpose methylation data
meth = t(meth)

# Resdualise meth matrix 
for(i in 1:(ncol(meth))){
  print(i)
  meth[,i] <- resid(lm(meth[,i] ~ phenotypes$age_stradl ++ as.factor(phenotypes$sex) + as.factor(phenotypes$Batch) + as.factor(phenotypes$wave) + as.factor(phenotypes$combined) + phenotypes$Bcell + phenotypes$CD4T + phenotypes$CD8T + phenotypes$Gran + phenotypes$NK, na.action = na.exclude))
} 

meth<-as.data.frame(meth)

## Combined IDs of individuals 
id <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/IDs_778.csv")
fam <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/GWAS_Somalogic_data.fam")

# Subset to 778
fam1 = fam[fam$V2 %in% id$x,]
## Make sure order of IDs in Family ID and regular ID file match 
matcher = id$x
fam1 = fam1[match(matcher,fam1$V2),]
table(fam1$V2 == id$x)

# match meth file to the order in the 778_id file 
ids = id$x
meth <- meth[match(ids, row.names(meth)),]
table(fam1$V2 == row.names(meth))

## Combine IDs for FID and IID and rearrange columns so it goes FID, IID, cpg1, cpg2... 
tmp <- meth
tmp$IID <- fam1$V2 
tmp$FID <- fam1$V1 
tmp <- tmp[,c(ncol(tmp), ncol(tmp)-1, 1:(ncol(tmp)-2))]
tmp$FID <- as.character(tmp$FID)
tmp$IID <- as.character(tmp$IID)

# test3 <- tmp[c(1:3)]
# identical(test3$IID, as.character(id$x)) # true 

# Add GS_ to IDs 
tmp$FID<-sub("^","GS_", tmp$FID)
tmp$IID<-sub("^","GS_", tmp$IID)

# This file has everything regressed 
fwrite(tmp, file="/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521.txt", row.names=F, quote = F, sep=' ')

# Save order of IDs from the methylation file 
test3 <- tmp[c(1:3)]
write.csv(test3, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/record_of_IDS_for_methylation_778_meth_order_full_280521_v2.csv")

# # meth file basic order 
# o1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/record_of_IDS_for_methylation_778_meth_order_basic_280521_v2.csv")
# # meth file full order 
# o2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/record_of_IDS_for_methylation_778_meth_order_full_280521_v2.csv")
# # Check they match up in terms of IDs 
# identical(o1$X, o2$X)
# identical(o1$FID, o2$FID)
# identical(o1$IID, o2$IID)

############################################################################################################

### Phenotype file - 778

############################################################################################################

# Adapted from the fully adjusted model phenotype script, but with the pQTLs regressed phenotype file

# Read in order of meth files we are matching to 
order <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/record_of_IDS_for_methylation_778_meth_order_full_240521_v2.csv")
# order_old <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/record_of_IDS_for_methylation_778_meth_order.csv")

pheno2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778_eGFR_removed_pQTLs_regressed.csv", check.names = F)
pheno2 <- pheno2[-2]


# Add FAM and ID columns 
id <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/IDs_778.csv") # 778 people in EWAS
fam <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/GWAS_Somalogic_data.fam") # 1064 total LBC family IDs included, with other info 

## Make sure order of IDs in Family ID and regular ID file match 
# So the second column of the fam file is the participant ID, which we will make sure matches here 
fam1 = fam[fam$V2 %in% id$x,] # subset the fam file to the participant IDs which are in the ID file 
matcher = id$x # create order based on the id file 
fam1 = fam1[match(matcher,fam1$V2),] # match order of participant IDs between the 2 files (first column is fam ID)

# Check to make sure the id in the proteins file is the same as the id in the fam1 file 
identical(pheno2$id, fam1$V2) # TRUE

# Name properly 
osca_dat <- fam1[c(1,2)]
names(osca_dat) <- c("FID", "IID")

pheno2$IID <- osca_dat$IID
pheno2$FID <- osca_dat$FID

# Order so proteins are first 
pheno2 <- pheno2[c(2:4236,1,4237,4238)]

# Add GS_ to the IDs 
pheno2$FID<-sub("^","GS_", pheno2$FID)
pheno2$IID<-sub("^","GS_", pheno2$IID)

identical(pheno2$IID, order$IID)

# Make sure the order of the phenotype file is matched to the main meth file order 
pheno2 <- pheno2[match(order$IID, pheno2$IID),]

identical(pheno2$IID, order$IID)

# Write out phenotype files for a heritability run in the 4235 with eGFR included as covariate (2 locations - one with WBC and one without)

location <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/"

## Write out protein files so FID, IID, phenotype 
for(i in 1:350){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch1/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 351:700){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch2/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 701:1050){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch3/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 1051:1400){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch4/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 1401:1750){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch5/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 1751:2100){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch6/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 2101:2450){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch7/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 2451:2800){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch8/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 2801:3150){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch9/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 3151:3500){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch10/",name,".phen"), col.names = F, row.names=F, quote = F, sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 3501:3850){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch11/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}

## Write out protein files so FID, IID, phenotype 
for(i in 3851:4235){ 
name <- names(pheno2)[i]
write.table(pheno2[,c("FID","IID",name)], file = paste0(location, "batch12/",name,".phen"), col.names = F, row.names=F, quote = F,sep=' ')
}



# ############################################################################################################

# ### MAKE BINARY FILES AND ORM - 778 

# ############################################################################################################

# All preps follow the WBC adjusted basic model script.


### Version one 

# Make Binary Methylation File
osca_Linux --efile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521.txt --methylation-m --make-bod --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521

# Regenerate the opi file with correct annotations 
anno = readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
opi1 = read.table("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521.opi", header=F, stringsAsFactors=F)
opi <- anno[opi1$V2,c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
opi$chr <- gsub("chr", "", opi$chr)
opi$chr <- gsub("X", "23", opi$chr)
opi$chr <- gsub("Y", "24", opi$chr)
opi$chr <- as.numeric(opi$chr)
opi$UCSC_RefGene_Name  <- as.factor(opi$UCSC_RefGene_Name )
opi$strand <- as.factor(opi$strand)
opi[which(opi$UCSC_RefGene_Name ==""), "UCSC_RefGene_Name"] <- NA
write.table(opi, file="/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521.opi",  # make sure to keep the same filename as before
                col.names=F, 
                row.names=F, 
                quote=F, sep='\t')

head /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521.opi


# Make ORM from methylation data 
osca_Linux --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --make-orm --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521


###########################################################################################

# BATCH EWAS 

cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch1/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch2/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch3/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch4/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch5/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch6/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch7/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch8/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch9/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch10/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch11/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



cd /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch12/
for i in *.phen 

do

  A=$( echo $i | cut -d"/" -f3)
  B=$( echo $A | cut -d. -f1)

osca_Linux --moa --orm /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/relatedness_testing_140421/daniel_orm/ormtest778 --befile /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080421/test_files_comparison_090421/methylation_778_meth_order_full_280521 --pheno $i --fast-linear --out ${B}_MOA_778_with_ORM --methylation-m


done 



####################################################################################################

### PROCESSING THE EWAS RESULTS FOR BATCHES ABOVE - MOA 778 with genetic ORM - WBC - AND pQTLs - but no smoking/BMI

####################################################################################################

# Split this up for speed into 12 screens 

# Go through each protein and subset to the significant sites 

library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch1/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"

batch <- c("b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11", "b12")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch2/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"


batch <- c("b2")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch3/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"

batch <- c("b3")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}



library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch4/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"


batch <- c("b4")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}



library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch5/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"


batch <- c("b5")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch6/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"


batch <- c("b6")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}



library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch7/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"


batch <- c("b7")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch8/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"


batch <- c("b8")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}



library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch9/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"


batch <- c("b9")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}



library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch10/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"


batch <- c("b10")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}



library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch11/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"


batch <- c("b11")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}



library(data.table)

list <- c("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/batch12/")

sig_hits <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"


batch <- c("b12")

for(i in 1:length(list)){

  path <- list[i]
  batch <- batch[i]
  setwd(paste0(path))

  L <- list.files(".", ".mlma")

  # results tables for each protein (all vs sig)
  for(j in 1:length(L)){
    file <- fread(L[j], header = T)
    file <- as.data.frame(file)
    name <- as.character(L[j])
    name <- gsub("_.*", "", name)
    file$SeqId <- name
    write.csv(file, paste0(path, name, "_results_table.csv"), row.names = F)
    sub <- file[-which(file$p > 0.000000036),]
    write.csv(sub, paste0(path, name, "_results_sig.csv"), row.names = F)
    write.csv(sub, paste0(sig_hits, batch, "_", name, "_results_sig.csv"), row.names = F)
  }
}


# Collate top hits for proteins 

path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"

setwd(path)

L <- list.files(".", ".csv")
L # 4230 converged 

files <- lapply(L, read.csv)
names <- as.character(L)
batch <- gsub("_.*", "", names)
marker <- gsub("_results.*", "", names)
marker <- gsub(".*_", "", marker)
names(files) <- marker
osca <- do.call(rbind, files)
osca <- osca[c(9,1:8)]

# > dim(osca)
# [1] 8834       - number of pQTL adjusted pQTMs 

length(unique(osca$SeqId))

# [1] 529 proteins with signals 

# length(files) # 4230

# annotations file add in
library(tidyverse) 
library(readxl)
anno <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/heritability_050321/outputs_combined/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
anno <- as.data.frame(anno)
anno <- anno[c(1,18,4,13)]
comb <- left_join(osca, anno, by = "SeqId")

# # Rank to see top proteins 
# comb2 <- comb[order(comb$p),]

# Save 
write.csv(comb, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_240521_results/no_eGFR_MIDDLE_EWAS_thr_3_6.csv", row.names = F)

# Filter to CpG-protein correction value of adjusted p 0.05/4235/772619 = 1.528098e-11

# > 0.05 / 4235 / 772619
# [1] 1.528098e-11

# 0.00000000001528098

comb2_filt2 <- comb[which(comb$p < 0.00000000001528098),]

# Save 
write.csv(comb2_filt2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_240521_results/no_eGFR_MIDDLE_EWAS_thr_bonfer.csv", row.names = F)

# I will use the latter as the correction threshold for now 

dim(comb2_filt2) # 2974  - at cpg/protein adjusted significance 

# Order and save an ordered version
 comb2 <- comb2_filt2[order(comb2_filt2$p),]

 write.csv(comb2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_240521_results/no_eGFR_MIDDLE_EWAS_thr_bonfer_ordered.csv", row.names = F)

##############################################################

### WORK OUT WHICH MIDDLE PROTEINS DIDNT CONVERGE 

# Read in protein file which has GS id and stradl id, then all 4,235 proteins that have been preprocessed
pheno1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/proteins_4235_778_eGFR_included.csv", check.names = F)
pheno2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/Phenotype_file_778_eGFR_removed.csv", check.names = F)
names(pheno2)[1] <- "ID"
pheno2 <- pheno2[-2]
colnames(pheno2) <- colnames(pheno1)


path <- "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_110821_middle/batches/sig_hits/"

setwd(path)

L <- list.files(".", ".csv")
L # 4226 converged 

files <- lapply(L, read.csv)
names <- as.character(L)
batch <- gsub("_.*", "", names)
marker <- gsub("_results.*", "", names)
marker <- gsub(".*_", "", marker)
names(files) <- marker
osca <- do.call(rbind, files)
osca <- osca[c(9,1:8)]

# load pheno2 above and do this to work out difference 
list <- names(files)
list2 <- pheno2[-which(colnames(pheno2) %in% list)]

# Plots - of proteins look okay 
proteins <- colnames(list2)[1:8]
proteins

### PROTEINS THAT DIDNT CONVERGE IN MIDDLE MODEL:
# [1] "ID"       "15509-2"  "15584-9"  "19141-22" "4407-10"  "6402-8"   NA
# [8] NA







