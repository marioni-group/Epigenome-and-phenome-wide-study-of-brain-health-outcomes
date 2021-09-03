############################################################################################
############################### Plotting EWAS results ######################################
############################################################################################
############################################################################################

############################################################################################################

### DO PLOT FOR JOINT NEURO MARKERS, AT FDR CORRECTION - TRANS ONLY AND CIS ONLY 

############################################################################################################


all <- read.csv("Y:/Danni/Cluster_Files/OSCA/cistrans_table_unformatted_naming.csv")
# cis <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/no_eGFR_cis_FDR_46_assocs_filtered.csv")
# trans <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Interpretation/interpretation_090621/no_eGFR_trans_FDR_42_assocs_filtered.csv")

cpgs <- all

names(cpgs)[4] <- "Gene_of_Hit"
cpgs[which(cpgs$Gene_of_Hit %in% ""),"Gene_of_Hit"] <- "Unannotated"

## Make first bed file - positions of CpGss
bed1 <- matrix(nrow = nrow(cpgs), ncol = 4)
bed1 <- as.data.frame(bed1)
names(bed1) <- c("chr", "start", "end", "value1")

## Update file naming 
names(cpgs)[1] <- "Chromosome_of_Hit"
names(cpgs)[3] <- "Position_of_Hit"

## Set up loop to populate bed1 file 
list = cpgs$Probe
for(i in 1:length(list)){ 
## get chromosome 
chr = cpgs[i, "Chromosome_of_Hit"]
chr = paste0("chr", chr)
## get start/end 
start.end <- cpgs[i, "Position_of_Hit"]
## value = beta
value <- cpgs[i, "b"]
## populate the table
bed1[i, "chr"] <- chr
bed1[i, "start"] <- start.end
bed1[i, "end"] <- start.end+5e5
bed1[i, "value1"] <- value
}

## Create bed2 file - target of interest (genes of interest)
bed2 <- matrix(nrow = nrow(cpgs), ncol = 3)
bed2 <- as.data.frame(bed2)
names(bed2) <- c("chr", "start", "end")

# Update naming for somamers 
names(cpgs)[12] <- "Chromosome_of_Somamer"
names(cpgs)[10] <- "Position_of_Somamer"
names(cpgs)[9] <- "Gene_of_Somamer"

# Get the info for somamers 
for(i in 1:length(list)){ 
## get chromosome 
chr = cpgs[i, "Chromosome_of_Somamer"]
chr = paste0("chr", chr)
## get start/end 
start.end <- cpgs[i, "Position_of_Somamer"]
## populate the table
bed2[i, "chr"] <- chr
bed2[i, "start"] <- start.end
bed2[i, "end"] <- start.end+5e5
}

## Create annotation file 
names1 <- bed1
names2 <- bed2
## give value1 column in bed1 (now names1) gene names of hit + denote its a SNP 
names1$value1 <- cpgs$Gene_of_Hit
names1$value1 <- paste(names1$value1, "(CpG)")
## give value1 column in bed2 (now names2) gene names of Somamer + denote its a Somamer
names2$value1 <- 0
names2$value1 <- cpgs$Gene_of_Somamer
names2$value1 <- paste(names2$value1, "")
## create bed file - annotation 
bed = rbind(names1, names2)

# # Write out file and edit so that each CpG is only annotated to a single key gene
# write.csv(bed, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_circus/circos_240521/circus_g_dataset_cog_joint.csv", row.names = F)
# write.csv(bed, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/EWAS_4000/cis_trans_circus/circos_240521/circus_g_dataset_cog_joint_trans_only.csv", row.names = F)
write.csv(bed, "circ_bed.csv", row.names = F)


# Read it in for plotting for the report 
library(readxl)


bed <- read_excel("circ_bed_GDF.xlsx")
bed <- as.data.frame(bed)



pdf("GDF_circ.pdf", width = 15, height = 15)
## initalize plot 
circos.initializeWithIdeogram()
circos.genomicTrackPlotRegion(bed, ylim = c(0, 1), track.height = 0.1, bg.border = NA)
i_track = get.cell.meta.data("track.index") # remember this empty track, we'll come back
circos.genomicTrackPlotRegion(bed, ylim = c(0, 1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicText(region, value, y = 1, labels.column = 1,
                                                   facing = "clockwise", adj = c(1, 0.5),
                                                   posTransform = posTransform.text, cex = 0.6, padding =2.0)
                              }, track.height = 0.1, bg.border = NA)

tr_track = get.cell.meta.data("track.index") # position transformation track
# because `circos.genomicPosTransformLines` is implemented by
# `circos.trackPlotRegion`, it accepts `track.index` argument.
circos.genomicPosTransformLines(bed,
                                posTransform = function(region, value)
                                  posTransform.text(region, y = 1, labels = value[[1]],
                                                    cex = 1, padding = 1, track.index = tr_track),
                                direction = "inside", track.index = i_track
)

## link the trans CpGss to their target Somamer/gene 
circos.genomicLink(bed1,bed2, col = rand_color(nrow(bed1)), border= NA)
dev.off()
