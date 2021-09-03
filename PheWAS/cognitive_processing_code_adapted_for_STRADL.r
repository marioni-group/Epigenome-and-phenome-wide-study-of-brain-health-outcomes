#######################
### read in IQ data ###
#######################

cog = read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/cognitive.csv", header=T)
cog = read.csv("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/st_cognitive.csv", header=T)

cog$LM <- cog$logical_mem_1 + cog$logical_mem_2

##############################################
### recode 0 to NA for each cognitive test ###
##############################################

cog$LM[cog$LM==0] <- NA
cog$digit_symbol[cog$digit_symbol==0] <- NA
cog$vocabulary[cog$vocabulary==0] <- NA
cog$verbal_total[cog$verbal_total==0] <- NA

############################################################################
### recode points +- 3.5 SDs from the mean to NA for each cognitive test ### 
############################################################################

low = mean(cog$LM, na.rm=T) - 3.5*sd(cog$LM, na.rm=T)
high = mean(cog$LM, na.rm=T) + 3.5*sd(cog$LM, na.rm=T)
table(cog$LM < low | cog$LM > high)
cog$LM[cog$LM < low | cog$LM > high] <- NA

low = mean(cog$verbal_total, na.rm=T) - 3.5*sd(cog$verbal_total, na.rm=T)
high = mean(cog$verbal_total, na.rm=T) + 3.5*sd(cog$verbal_total, na.rm=T)
table(cog$verbal_total < low | cog$verbal_total > high)
cog$verbal_total[cog$verbal_total < low | cog$verbal_total > high] <- NA

low = mean(cog$digit_symbol, na.rm=T) - 3.5*sd(cog$digit_symbol, na.rm=T)
high = mean(cog$digit_symbol, na.rm=T) + 3.5*sd(cog$digit_symbol, na.rm=T)
table(cog$digit_symbol < low | cog$digit_symbol > high)
cog$digit_symbol[cog$digit_symbol < low | cog$digit_symbol > high] <- NA

low = mean(cog$vocabulary, na.rm=T) - 3.5*sd(cog$vocabulary, na.rm=T)
high = mean(cog$vocabulary, na.rm=T) + 3.5*sd(cog$vocabulary, na.rm=T)
table(cog$vocabulary < low | cog$vocabulary > high)
cog$vocabulary[cog$vocabulary < low | cog$vocabulary > high] <- NA


# # Do versions which arent scales 
# library(psych)

# cog$g1 = principal(cog[,c("digit_symbol","verbal_total","vocabulary","LM")], factors=1, rotate="none", na.action="na.exclude")$score
# cog$gf1 = principal(cog[,c("digit_symbol","verbal_total","LM")], factors=1, rotate="none", na.action="na.exclude")$score

# names(cog)[14:15] <- c("gf", "g")


#################################
### Create gf and g using PCA ###
#################################

library(psych)

cog$g1 = scale(principal(cog[,c("digit_symbol","verbal_total","vocabulary","LM")], factors=1, rotate="none", na.action="na.exclude")$score)
cog$gf1 = scale(principal(cog[,c("digit_symbol","verbal_total","LM")], factors=1, rotate="none", na.action="na.exclude")$score)

names(cog)[14:15] <- c("gf", "g")

cog1 <- cog[,c("st","digit_symbol","verbal_total","vocabulary","LM","gf","g")]

write.csv(cog1, "/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/Data_for_inputs/cog1_240321.csv", row.names = F)


### Make versions of g and gf that arent scaled 
# /Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/pheWAS/pheWAS_240521/APOE_cog_decline_310521/



              PC1   h2   u2 com
digit_symbol 0.58 0.33 0.67   1
verbal_total 0.72 0.51 0.49   1
vocabulary   0.67 0.45 0.55   1
LM           0.69 0.47 0.53   1

                PC1
SS loadings    1.76
Proportion Var 0.44


              PC1   h2   u2 com
digit_symbol 0.74 0.55 0.45   1
verbal_total 0.68 0.46 0.54   1
LM           0.72 0.51 0.49   1

                PC1
SS loadings    1.52
Proportion Var 0.51


