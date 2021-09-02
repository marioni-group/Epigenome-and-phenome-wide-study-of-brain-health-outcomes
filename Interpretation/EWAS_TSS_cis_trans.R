############################################################################################
############################################################################################
################# Lookup of cis trans effects for EWAS #####################################
############################################################################################
############################################################################################

# Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_circus/")

# Libraries
library(tidyverse)
library(readxl)

# Load in annotational information for the STRADL proteins
library(readxl)
anno <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/heritability_050321/outputs_combined/SOMAscan_Assay_v4_Annotations_version3.3.2_DG.xlsx")
anno <- as.data.frame(anno)
anno <- anno[c(1,18,4,13)]
# comb <- left_join(osca, anno, by = "SeqId")

# Load in the significant EWAS results with proteins listed to get unique proteins
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080821_full/results/no_eGFR_FULL_EWAS_thr_bonfer.csv")
unique <- unique(cpgs$SeqId)

length(unique)
# [1] 153

unique <- cpgs[c(1,10,11)] %>% unique() # xtracted into dataset with gene names 
tab <- unique
names(tab)[2] <- "gene" 

# get gene info set up from external source (Use Ensembl to Extract Relevant Features)
library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl, page = "structure")
attributes[grep("transcript", attributes$description, ignore.case = TRUE), ]

# get list of proteins gene names 
list <- tab$gene

## Get Transcription Start Sites for soma proteins 
tss <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", values = c(list),
             mart = ensembl)

# see how many this worked for
length(unique(tss$external_gene_name)) # 145

# See if there are unique TSS extracted for these proteins 
dim(tss) # 1420
length(unique(tss$transcription_start_site)) # 1162

# which are the gene names that didnt extract
ov <- which(!list %in% tss$external_gene_name)
list2 <- list[ov]


# [1] "CASC4"   "HDGFRP3" "NT5C3L"  "CECR1"

# [1] "HDGFRP3" "NT5C3L"  "CECR1"   "CASC4"

# HDGFL3

# NT5C3B

# ADA2 

# GOLM2 

# Try these aliases for the genes as per lookup on gene cards site 

list3 <- c("HDGFL3","NT5C3B","ADA2","GOLM2")

## Get Transcription Start Sites for remaining soma proteins 
tss2 <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", values = c(list3),
             mart = ensembl)

# see how many this worked for
length(unique(tss2$external_gene_name)) # 4

# Now we have a complete extracted set of tss information - join them together 

tss <- rbind(tss, tss2)

length(unique(tss$external_gene_name)) # 149

## set up code to create biomart dataframe with gene, start and end columns
biomart_df <- tss

names(biomart_df)[8] <- "gene"
names(biomart_df)[3] <- "start"
names(biomart_df)[4] <- "end"
names(biomart_df)[2] <- "chromosome_ensembl"

# Again, lets do a check to see if each protein has its own TSS extracted that is unique 

subset <- biomart_df[,c("transcription_start_site", "gene")]
subset <- unique(subset$gene)

## separate data into genes on positive strand and genes on negative strand

biomart_df_pos <- biomart_df[biomart_df$strand == 1,]
biomart_df_neg <- biomart_df[biomart_df$strand == -1,]

# sanity check to see if overlap is present/correct number of genes are present

a = unique(biomart_df_pos$gene)
b = unique(biomart_df_neg$gene)
which(a %in% b)
length(a)
length(b)

## create output data frames for genes on positive strand and those on negative strand to be r-bound later
   
   
    ##out_df1 == genes on positive strand
   
out_df1 <- matrix(nrow=length(unique(biomart_df_pos$gene)), ncol=2)
colnames(out_df1) <- c("start", "end")
rownames(out_df1) <- unique(biomart_df_pos$gene)


    ##out_df2 == genes on negative strand
   
out_df2 <- matrix(nrow=length(unique(biomart_df_neg$gene)), ncol=2)
colnames(out_df2) <- c("start", "end")
rownames(out_df2) <- unique(biomart_df_neg$gene)

## loop through gene names to fill in min 5'end (start) and max 3'end (end) as on positive strand you want furthest away 5'(tss) and biggest distance to 3' end
## fill out out_df1 based on the above

for(gene in unique(biomart_df_pos$gene)) {
tmp <- biomart_df_pos[which(biomart_df_pos$gene==gene), ]
out_df1[gene,"start"] <-  min(tmp$start)
out_df1[gene,"end"] <-  max(tmp$end)
}

## loop through gene names to fill in min 3'end (end) and max 5'end (start) as on negative strand you want furthest away 3'(tss (relative to positive strand)) and biggest distance to 5' end
## fill out out_df2 based on the above

for(gene in unique(biomart_df_neg$gene)) {
tmp <- biomart_df_neg[which(biomart_df_neg$gene==gene), ]
out_df2[gene,"start"] <-  max(tmp$end)
out_df2[gene,"end"] <-  min(tmp$start)
}

## row bind the positive and negative strand dataframes to give one out_df
out_df <- rbind(out_df1, out_df2)


## Tidy the files before saving for new row with protein names (out_df), then one for length of transcripts in biomart_df for future reference

    ## set as data frames
    out_df <- as.data.frame(out_df)
    biomart_df <- as.data.frame(biomart_df)

out_df$Gene <- row.names(out_df)
biomart_df$transcript_length <- abs(biomart_df$start - biomart_df$end)

# Check unique start sites here 
unique(out_df$start) # 149
plot(out_df$start)

## Get Chromosome Number for Protein 
biomart_df$chromosome_ensembl2 <-  gsub("CHR_HSCHR", "", biomart_df$chromosome_ensembl)
biomart_df$chromosome_ensembl2 <- gsub("_.*", "", biomart_df$chromosome_ensembl2)
biomart_df <- biomart_df[,c(2,8,10)]
names(biomart_df) <- c("Chromosome_Biomarker", "Gene", "Chromosome_edit")
biomart_df = biomart_df[-which(duplicated(biomart_df$Gene)),]

# Join biomart_df file to the out_df file 
out_df <- merge(out_df, biomart_df, by = "Gene")

# # Check unique start sites here 
# unique(out_df$start) # 154
# plot(out_df$start)

names(out_df)[2] <- "Biomarker_TSS_Start"
names(out_df)[3] <- "Biomarker_TSS_End"
names(out_df)[1] <- "gene"

# Eyeball to make sure chr conversion has occured - fix any remaining formatting issues 
out_df$Chromosome_edit <- gsub("KIR", "", out_df$Chromosome_edit)
out_df$Chromosome_edit <- gsub("LRC", "", out_df$Chromosome_edit)
out_df$Chromosome_Biomarker <- NULL

# Replace the 4 trouble genes with the somamer name for joining as they are now extracted 
out_df$gene <- gsub("HDGFL3", "HDGFRP3", out_df$gene)
out_df$gene <- gsub("NT5C3B", "NT5C3L", out_df$gene)
out_df$gene <- gsub("ADA2", "CECR1", out_df$gene)
out_df$gene <- gsub("GOLM2", "CASC4", out_df$gene)

# [1] "HDGFRP3" "NT5C3L"  "CECR1"   "CASC4"
# list3 <- c("HDGFL3","NT5C3B","ADA2","GOLM2")

# Merge back into main dataset with sig proteins 
tab2 <- left_join(tab, out_df, by = "gene") 

table(is.na(tab2))
# FALSE
#   918

write.csv(tab2, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_090621/no_eGFR_TSS_table.csv", row.names = F)


#############################################################################################

# Now look at mapping CIS and TRANS associations for the cpgs that associated with the proteins 

# Read in the cpg positional info 
cpgs <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/EWAS_080821_full/results/no_eGFR_FULL_EWAS_thr_bonfer.csv")

# Read in the protein gene annotations as generated above 
prot <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_090621/no_eGFR_TSS_table.csv")

# Merge gene positional info into the table with weights based on the gene names 
table <- left_join(cpgs, prot, by = c("SeqId"))
dim(table) # 2895 - this looks good 

# See whether we still have unique TSS sites 
length(unique(table$Biomarker_TSS_Start)) # 149 protein start sites 
length(unique(table$bp)) # 1760 cpg sites
length(unique(table$Probe)) # 1760 cpgs

# Align column names to cis trans code 

chrom <- table
names(chrom)[15] <- "gene_start"
names(chrom)[16] <- "gene_end"
names(chrom)[4] <- "Cpg_pos"

chrom$gene_start <- as.numeric(chrom$gene_start)
chrom$gene_end <- as.numeric(chrom$gene_end)

# Look at which matches are above and below in the range 
chrom$diff <- chrom$gene_start - chrom$Cpg_pos 
chrom$diff2 <- chrom$gene_start + chrom$Cpg_pos

# Change negative differences to absolute values 
chrom$diff <- abs(chrom$diff)
chrom$diff2 <- abs(chrom$diff2)

# match chromosome naming columns
names(chrom)[2] <- "chr_cpg"
names(chrom)[17] <- "Chr"

# Do everything in one go to assign CIS or TRANS 
chrom <- chrom %>% mutate(Effect = case_when(
  chr_cpg == Chr & (chrom$diff < 10000000 | chrom$diff2 < 10000000) ~ "CIS",
  chr_cpg == Chr & (!chrom$diff < 10000000 | chrom$diff < 10000000) ~ "TRANS",
  !chr_cpg == Chr ~ "TRANS"
))

# Make sure you get same results as decomposed version below (yes, it matches up)
test <- chrom %>% filter(Effect == "CIS") 
dim(test)
test <- chrom %>% filter(Effect == "TRANS") 
dim(test)


# Write out table for use in next script with basic naming 
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_090621/no_eGFR_cistrans_table_unformatted_naming.csv", row.names = F)


# Write the table for the suppl file
chrom <- chrom[c(3,2,5,4,6:9,1,10,11,12,15:18,20)]
names(chrom) <- c("CpG", "Chromsome of CpG", "Gene of CpG", "CpG Position", "Orientation",
  "Beta", "SE", "P", "SeqId", "Gene of Protein", "UniProt", "UniProt Full Name", "Gene Start", "Gene End", "Chromosome of Gene", "Association Type")
chrom <- chrom[order(chrom$P),]
write.csv(chrom, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_090621/no_eGFR_cistrans_table.csv", row.names = F)




