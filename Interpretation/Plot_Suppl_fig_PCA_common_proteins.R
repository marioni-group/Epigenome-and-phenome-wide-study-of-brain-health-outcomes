

### PLOT FOR COMMON PROTEINS PCA 

library(tidyverse)
library(ggplot2)
library(readxl)
library(psych)
library(ggcorrplot)
library(cowplot)
library(tidyverse)

# Read in the protein file 
prot <- read.csv("Y:/Danni/Cluster_Files/prot_file_110821_WMHV.csv", check.names = F)

# Read in "common" file from the updated beta plots script for the plot with 25 proteins 
# List the top 25 proteins (26 somamers)
list <- unique(common$SeqId)

# subset prot file to include just these somamers 
prot2 <- prot[,which(colnames(prot) %in% list)]

# Assign variable to joint 
joint <- prot2

# 4876-32 - F9 
# 5307-12 - F9 

joint <- joint[,-which(colnames(joint) %in% "5307-12")]

# Scale data and get PC scores
scores_pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))$scores
variance_pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))$Vaccounted
pca <- principal(scale(joint), rotate="none", nfactors=ncol(joint))

scores <- as.data.frame(scores_pca)
var1 <- as.data.frame(variance_pca)

var <- var1[1,] # get variance 
cum <- var1[5,] # get cumulative variance

var <- gather(var) # gather
var$col <- ifelse(var$value >= 1, "darkgrey", "orange")
var$num <- as.integer(1:ncol(joint))
var$mes <- as.numeric(var$value)

cum <- gather(cum) # gather
cum$num <- as.integer(1:ncol(joint))
cum$mes <- as.numeric(cum$value)

cor = cor(joint)

## Read in seq-id conversion file 
anno <- read_excel("Y:/Protein_DNAm_Proxies/Annotations_for_reference.xlsx")
anno <- as.data.frame(anno)
anno <- anno[,c(1,4,18)]

## subset seq-ids
anno1 = anno[which(anno$SeqId %in% colnames(cor)),] 
## match up their order 
ids = colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] 
anno1 = anno1[match(ids, anno1$SeqId),]

## check they match
table(anno1$SeqId == colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] )

## replace seq-ids with gene names 
names(anno1)[3] <- "Name"
colnames(cor)[which(colnames(cor) %in% anno1$SeqId)] <- as.character(anno1$Name)
row.names(cor)[which(row.names(cor) %in% anno1$SeqId)] <- as.character(anno1$Name)

# a <- if(ncol(joint) > 12) { 
a <- ggplot(var, aes(num, mes, col = col)) +
geom_point(size = 4) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
scale_color_manual(values=c("#E69F00", "#999999")) 
# } else { 
  
# ggplot(var, aes(num, mes, col = col)) +
# geom_point(size = 4) + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
# xlab("Principal component") + ylab("Eigenvalue") + geom_hline(yintercept=1, color = "darkgrey", size=1) + 
# scale_color_manual(values=c("#E69F00", "#999999")) + scale_x_continuous(n.breaks = ncol(joint))
# # }


# b <- if(ncol(joint) > 12) { 
b <- ggplot(cum, aes(num, mes)) +
geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
xlab("Principal component") + ylab("Cumulative proportion") + ylim(0,1) 
# } else if(ncol(joint) > 4) { ggplot(cum, aes(num, mes)) +
#   geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
#   xlab("Principal component") + ylab("Cumulative proportion") + ylim(0,1) +  scale_x_continuous(n.breaks = ncol(joint))} else { 
#     ggplot(cum, aes(num, mes)) +
#       geom_bar(stat = "identity", fill = "steelblue") + theme(legend.position = "none", axis.text=element_text(size=13), axis.title=element_text(size=13)) +
#       xlab("Principal component") + ylab("Cumulative proportion") + ylim(0,1) 
    
#     }

# # c <- if(ncol(joint) > 12){  
# c <- ggcorrplot(cor, 
#            hc.order = TRUE)

c <- ggcorrplot(cor, 
           hc.order = TRUE,
           type = "lower")


c <- c + theme(
  panel.background = element_rect(fill = "white", colour = "white",
                                size = 2, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white")
  )

#            lab_size = 1.5,
#            colors = c("blue", "white", "red")) +  theme(legend.title = element_blank(), text = element_text(size = 11), axis.text.x=element_text(size=9), axis.text.y=element_text(size=9))


# +labs(title = list_diseases[i]) + 
 
# } else { 
#   ggcorrplot(cor, 
#              hc.order = TRUE, 
#              lab = TRUE,
#              type = "lower",
#              lab_size = 2.5,
#              colors = c("blue", "white", "red")) +labs(title = list_diseases[i]) + 
#     theme(legend.title = element_blank(), text = element_text(size = 11), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
  
#   } 

p1 = plot_grid(c,a,b, nrow = 1, labels = c("a", "b", "c"), rel_widths = c(0.85,0.5,0.5))
# list_plots[[i]] <- p1

# print(i)
# } 

pdf("Y:/Danni/stradl_markers/Plots/SUPPL_PCA/FIGURE_lower_25_proteins.pdf", height =6 , width = 17)
p1
dev.off()

# pdf("U:/Protein_DNAm_Proxies/Work_and_code_post_KORA/Supp_Figs/Fig3_Diseases1-3.pdf", height = 11, width = 10.5)
# plot_grid(list_plots[[1]],list_plots[[2]],list_plots[[3]], nrow = 3)
# dev.off()





