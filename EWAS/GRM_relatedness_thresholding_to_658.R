############################################################################################

## Exclude those who are related based on 5% upper limits of GRM in STRADL 

############################################################################################

# Daniel sent this location for the GRM to use:
# /Cluster_Filespace/Marioni_Group/Daniel/Somalogic_GWAS/grm_1000_full

# Prune the GRM for relatedness by a cutoff of 0.05
gcta64 --grm /Cluster_Filespace/Marioni_Group/Daniel/Somalogic_GWAS/grm_1000_full --grm-cutoff 0.05 --make-grm --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/trimming_GRM_180321/grm_pruned_0.05_related

# 181 inidviduals removed from the GRM and 884 kept from the 1000 we had 

# Extract the GRM subject id of all the singletons by a cutoff of 0.05
gcta64 --grm /Cluster_Filespace/Marioni_Group/Daniel/Somalogic_GWAS/grm_1000_full --grm-singleton 0.05  --out /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/trimming_GRM_180321/grm_singletons

# Save IDs for the 884 unrelated in STRADL 
test <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/trimming_GRM_180321/grm_singletons.singleton.txt", header = F) # this is the IDs of all my 884 single people 



### NOW GET A LIST OF THE 658 UNRELATED INDIVIDUALS TO SUBSET BY 

# read in protein example dataset with GS id listed for the 778 individuals 
pheno2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/prep/proteins_the_whole_lot_total.csv", check.names = F)

# Add FAM and ID columns 
id <- read.csv("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/IDs_778.csv") # 778 people in EWAS
fam <- read.table("/Cluster_Filespace/Marioni_Group/Rob/Somalogic/GWAS_Somalogic_data.fam") # 1064 total LBC family IDs included, with other info 

## Make sure order of IDs in Family ID and regular ID file match 
# So the second column of the fam file is the participant ID, which we will make sure matches here 
fam1 = fam[fam$V2 %in% id$x,] # subset the fam file to the participant IDs which are in the ID file 
matcher = id$x # create order based on the id file 
fam1 = fam1[match(matcher,fam1$V2),] # match order of participant IDs between the 2 files (first column is fam ID)

# Check to make sure the id in the proteins file is the same as the id in the fam1 file 
# Check to make sure the id file order matches the pheno2 order 
identical(pheno2$ID, id$x) # TRUE 
identical(pheno2$ID, fam1$V2) # TRUE

# Name properly 
osca_dat <- fam1[c(1,2)]
names(osca_dat) <- c("FID", "IID")

pheno2$IID <- osca_dat$IID
pheno2$FID <- osca_dat$FID

identical(pheno2$ID, pheno2$IID) # TRUE 

# Order so proteins are first 
pheno2 <- pheno2[c(2:4236,1,4237,4238)]

# Check subset to make sure matches correctly 
test <- pheno2[c(4235:4238)]

# Subset phenotype file to only those that are unrelated (658) 
# read in the unrelated IIDs generated above 
test <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/trimming_GRM_180321/grm_singletons.singleton.txt", header = F) # this is the IDs of all my 884 single people 
ov <- which(pheno2$IID %in% test$V2)
pheno2 <- pheno2[ov,]

# > dim(pheno2)
# [1]  658 4238 

# Check to make sure its matched okay 
test2 <- pheno2[c(4235:4238)]

# write out the IDs for the 658 people so the protein and meth files can be filtered in prep steps 
ID <- pheno2[c(4237:4238)]
write.csv(ID, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/658_unrelated_IDs.csv", row.names = F)

