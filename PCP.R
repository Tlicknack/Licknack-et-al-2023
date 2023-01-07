  # Paramecium PCP: PCP.R

library(tidyverse)
library(pRoloc)
library(reshape2)

### Important objects:
load("/ParameciumPCP/RData/ptetPCP_superANDunsuper.RData")
###

#-------------------------------------------------------------------------------------------------------------------------
##### Raw Data Processing

setwd("/ParameciumPCP/")

ptetPCP = read.csv(file = "/ParameciumPCP/tables/ptetPCP.csv", header=T) 
dim(ptetPCP)  # 12579 proteins X 248 columns

### Remove cRaP/Kleb Proteins
ptetPCP = rbind(ptetPCP[grep(pattern = "PTET*", x = ptetPCP$Accession),],   
                ptetPCP[grep(pattern = "*lcl", x = ptetPCP$Accession),])   
dim(ptetPCP)  # 11887 proteins X 248 columns (34 proteins removed)

# Write fasta
library(seqinr)
ptProts = read.fasta("/Paramecium_Proteomes-Predicted/ptetraurelia_mac_51_annotation_v2.0.protein.fa", seqtype = "AA", as.string = T)
# write.fasta(sequences = getSequence(ptProts[ptetPCP[grep(pattern = "PTET*", x = ptetPCP$Accession),'Accession']], as.string = T), 
#             names = ptetPCP[grep(pattern = "PTET*", x = ptetPCP$Accession),'Accession'], 
#             file.out = "/ParameciumPCP/sequences/ptetPCP-ALLPROTEINSUNFILTERED.fasta")

### Change Mitochondrial ORF names  (ORFinder output is weird)
ptetPCP[grep(pattern = "*lcl", x = ptetPCP$Accession),'Accession'] = ptetPCP[grep(pattern = "*lcl", x = ptetPCP$Accession),'Description']
ptetPCP[grep(pattern = "ORF*", x = ptetPCP$Accession),'Accession'] = unlist(lapply(ptetPCP[grep(pattern = "ORF*", x = ptetPCP$Accession),'Accession'], strsplit, " "))[grep(pattern = "ORF*", x = unlist(lapply(ptetPCP[grep(pattern = "ORF*", x = ptetPCP$Accession),'Accession'], strsplit, " ")))]
ptetPCP[grep(pattern = "ORF*", x = ptetPCP$Accession),'Accession'] = unlist(lapply(strsplit(ptetPCP[grep(pattern = "ORF*", x = ptetPCP$Accession),'Accession'], ":"), "[", 1))

### Remove Low Confidence
ptetPCP = ptetPCP[-which(ptetPCP$Protein.FDR.Confidence..Combined == "Low"),]
dim(ptetPCP)  # 11595 proteins X 248 columns (333 proteins removed)

### Remove Proteins with < 10 PSMs (arbitrary, but 3 is too low!)
ptetPCP = ptetPCP[-which(ptetPCP$X..PSMs < 10),]
dim(ptetPCP)  # 9545 proteins X 248 columns (2317 proteins removed)

### Remove Proteins without Unique Peptides 
ptetPCP = ptetPCP[-which(ptetPCP$X..Unique.Peptides == 0),]
dim(ptetPCP)  # 9544 proteins X 248 columns (1 protein removed)

nrow(ptetPCP[grep("PTET*", ptetPCP$Accession),])  # 9513 p. tetraurelia proteins
nrow(ptetPCP[grep("ORF*", ptetPCP$Accession),])   # 31 mito-ORFS

### Raw abundance file
# Notes: 1 5K is missing (corrupt), 1 12K was rerun due to machine error (F146 --> F148)
ptetPCP_raw = ptetPCP %>% dplyr::select(colnames(ptetPCP)[c(2,1,4:9)], colnames(ptetPCP)[34:140])
dim(ptetPCP_raw)  # 9544 proteins by 115 columns
# write.csv(ptetPCP_raw, "/ParameciumPCP/tables/ptetPCP_raw.csv", row.names = F, quote = F)

### Split Replicate Experiments
ptetPCP_raw_MAC = ptetPCP_raw[, grep(pattern = "*Mac", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_MAC = as.data.frame(rowMeans(ptetPCP_raw_MAC[,c(1:3)]))
ptetPCP_raw_rep2_MAC = as.data.frame(rowMeans(ptetPCP_raw_MAC[,c(4:6)]))
ptetPCP_raw_rep3_MAC = as.data.frame(rowMeans(ptetPCP_raw_MAC[,c(7:9)]))

ptetPCP_raw_300g = ptetPCP_raw[, grep(pattern = "*300g", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_300g = as.data.frame(rowMeans(ptetPCP_raw_300g[,c(1:3)]))
ptetPCP_raw_rep2_300g = as.data.frame(rowMeans(ptetPCP_raw_300g[,c(4:6)]))
ptetPCP_raw_rep3_300g = as.data.frame(rowMeans(ptetPCP_raw_300g[,c(7:9)]))

ptetPCP_raw_1K = ptetPCP_raw[, grep(pattern = "*1K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_1K = as.data.frame(rowMeans(ptetPCP_raw_1K[,c(1:3)]))
ptetPCP_raw_rep2_1K = as.data.frame(rowMeans(ptetPCP_raw_1K[,c(4:6)]))
ptetPCP_raw_rep3_1K = as.data.frame(rowMeans(ptetPCP_raw_1K[,c(7:9)]))

ptetPCP_raw_3K = ptetPCP_raw[, grep(pattern = "*3K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_3K = as.data.frame(rowMeans(ptetPCP_raw_3K[,c(1:3)]))
ptetPCP_raw_rep2_3K = as.data.frame(rowMeans(ptetPCP_raw_3K[,c(4:6)]))
ptetPCP_raw_rep3_3K = as.data.frame(rowMeans(ptetPCP_raw_3K[,c(7:9)]))

ptetPCP_raw_5K = ptetPCP_raw[, grep(pattern = "*..Sample..5K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_5K = as.data.frame(rowMeans(ptetPCP_raw_5K[,c(1:3)]))
ptetPCP_raw_rep2_5K = as.data.frame(rowMeans(ptetPCP_raw_5K[,c(4:6)]))
ptetPCP_raw_rep3_5K = as.data.frame(rowMeans(ptetPCP_raw_5K[,c(7:8)]))  # Missing one file

ptetPCP_raw_9K = ptetPCP_raw[, grep(pattern = "*..Sample..9K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_9K = as.data.frame(rowMeans(ptetPCP_raw_9K[,c(1:3)]))
ptetPCP_raw_rep2_9K = as.data.frame(rowMeans(ptetPCP_raw_9K[,c(4:6)]))
ptetPCP_raw_rep3_9K = as.data.frame(rowMeans(ptetPCP_raw_9K[,c(7:9)]))

ptetPCP_raw_12K = ptetPCP_raw[, grep(pattern = "*12K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_12K = as.data.frame(rowMeans(ptetPCP_raw_12K[,c(1:3)]))
ptetPCP_raw_rep2_12K = as.data.frame(rowMeans(ptetPCP_raw_12K[,c(4:6)]))
ptetPCP_raw_rep3_12K = as.data.frame(rowMeans(ptetPCP_raw_12K[,c(7:9)]))

ptetPCP_raw_15K = ptetPCP_raw[, grep(pattern = "*15K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_15K = as.data.frame(rowMeans(ptetPCP_raw_15K[,c(1:3)]))
ptetPCP_raw_rep2_15K = as.data.frame(rowMeans(ptetPCP_raw_15K[,c(4:6)]))
ptetPCP_raw_rep3_15K = as.data.frame(rowMeans(ptetPCP_raw_15K[,c(7:9)]))

ptetPCP_raw_30K = ptetPCP_raw[, grep(pattern = "*30K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_30K = as.data.frame(rowMeans(ptetPCP_raw_30K[,c(1:3)]))
ptetPCP_raw_rep2_30K = as.data.frame(rowMeans(ptetPCP_raw_30K[,c(4:6)]))
ptetPCP_raw_rep3_30K = as.data.frame(rowMeans(ptetPCP_raw_30K[,c(7:9)]))

ptetPCP_raw_79K = ptetPCP_raw[, grep(pattern = "*79K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_79K = as.data.frame(rowMeans(ptetPCP_raw_79K[,c(1:3)]))
ptetPCP_raw_rep2_79K = as.data.frame(rowMeans(ptetPCP_raw_79K[,c(4:6)]))
ptetPCP_raw_rep3_79K = as.data.frame(rowMeans(ptetPCP_raw_79K[,c(7:9)]))

ptetPCP_raw_120K = ptetPCP_raw[, grep(pattern = "*120K", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_120K = as.data.frame(rowMeans(ptetPCP_raw_120K[,c(1:3)]))
ptetPCP_raw_rep2_120K = as.data.frame(rowMeans(ptetPCP_raw_120K[,c(4:6)]))
ptetPCP_raw_rep3_120K = as.data.frame(rowMeans(ptetPCP_raw_120K[,c(7:9)]))

ptetPCP_raw_SUP = ptetPCP_raw[, grep(pattern = "*Sup", colnames(ptetPCP_raw))]
ptetPCP_raw_rep1_SUP = as.data.frame(rowMeans(ptetPCP_raw_SUP[,c(1:3)]))
ptetPCP_raw_rep2_SUP = as.data.frame(rowMeans(ptetPCP_raw_SUP[,c(4:6)]))
ptetPCP_raw_rep3_SUP = as.data.frame(rowMeans(ptetPCP_raw_SUP[,c(7:9)]))

ptetPCP_raw_rep1 = cbind(ptetPCP_raw[,c(1:8)], ptetPCP_raw_rep1_MAC, ptetPCP_raw_rep1_300g, ptetPCP_raw_rep1_1K, ptetPCP_raw_rep1_3K, ptetPCP_raw_rep1_5K, ptetPCP_raw_rep1_9K, ptetPCP_raw_rep1_12K, ptetPCP_raw_rep1_15K, ptetPCP_raw_rep1_30K, ptetPCP_raw_rep1_79K, ptetPCP_raw_rep1_120K, ptetPCP_raw_rep1_SUP)
ptetPCP_raw_rep2 = cbind(ptetPCP_raw[,c(1:8)], ptetPCP_raw_rep2_MAC, ptetPCP_raw_rep2_300g, ptetPCP_raw_rep2_1K, ptetPCP_raw_rep2_3K, ptetPCP_raw_rep2_5K, ptetPCP_raw_rep2_9K, ptetPCP_raw_rep2_12K, ptetPCP_raw_rep2_15K, ptetPCP_raw_rep2_30K, ptetPCP_raw_rep2_79K, ptetPCP_raw_rep2_120K, ptetPCP_raw_rep2_SUP)
ptetPCP_raw_rep3 = cbind(ptetPCP_raw[,c(1:8)], ptetPCP_raw_rep3_MAC, ptetPCP_raw_rep3_300g, ptetPCP_raw_rep3_1K, ptetPCP_raw_rep3_3K, ptetPCP_raw_rep3_5K, ptetPCP_raw_rep3_9K, ptetPCP_raw_rep3_12K, ptetPCP_raw_rep3_15K, ptetPCP_raw_rep3_30K, ptetPCP_raw_rep3_79K, ptetPCP_raw_rep3_120K, ptetPCP_raw_rep3_SUP)
colnames(ptetPCP_raw_rep1) = c("Accession", "Confidence-FDR", "Coverage", "nPeptides", "nPSM", "nUniquePeptides","nAA", "nRazorPeptides", 
                               "MAC-1", "300g-1", "1K-1", "3K-1", "5K-1", "9K-1", "12K-1", "15K-1", "30K-1", "79K-1", "120K-1", "SUP-1")
colnames(ptetPCP_raw_rep2) = c("Accession", "Confidence-FDR", "Coverage", "nPeptides", "nPSM", "nUniquePeptides","nAA", "nRazorPeptides", 
                               "MAC-2", "300g-2", "1K-2", "3K-2", "5K-2", "9K-2", "12K-2", "15K-2", "30K-2", "79K-2", "120K-2", "SUP-2")
colnames(ptetPCP_raw_rep3) = c("Accession", "Confidence-FDR", "Coverage", "nPeptides", "nPSM", "nUniquePeptides","nAA", "nRazorPeptides", 
                               "MAC-3", "300g-3", "1K-3", "3K-3", "5K-3", "9K-3", "12K-3", "15K-3", "30K-3", "79K-3", "120K-3", "SUP-3")

### Remove + Output MAC or SUP only proteins
source("/ParameciumPCP/scripts/fun_naFirstLast.R")
ptetPCP_raw_rep1_MAConly = fun_naFirstLast(ptetPCP_raw_rep1, "first")  # 27 proteins removed
ptetPCP_raw_rep2_MAConly = fun_naFirstLast(ptetPCP_raw_rep2, "first")  # 27 proteins removed
ptetPCP_raw_rep3_MAConly = fun_naFirstLast(ptetPCP_raw_rep3, "first")  # 35 proteins removed
ptetPCP_raw_rep1 = ptetPCP_raw_rep1[-which(ptetPCP_raw_rep1$Accession %in% ptetPCP_raw_rep1_MAConly$Accession),]
ptetPCP_raw_rep2 = ptetPCP_raw_rep2[-which(ptetPCP_raw_rep2$Accession %in% ptetPCP_raw_rep2_MAConly$Accession),]
ptetPCP_raw_rep3 = ptetPCP_raw_rep3[-which(ptetPCP_raw_rep3$Accession %in% ptetPCP_raw_rep3_MAConly$Accession),]

ptetPCP_raw_MAConly = intersect(intersect(ptetPCP_raw_rep1_MAConly$Accession, 
                                          ptetPCP_raw_rep2_MAConly$Accession), 
                                          ptetPCP_raw_rep3_MAConly$Accession)  # 8 shared
# write.csv(ptetPCP_raw_MAConly, file = "/ParameciumPCP/tables/ptetPCP_MAConly.csv", quote = F, row.names = F)

ptetPCP_raw_rep1_SUPonly = fun_naFirstLast(ptetPCP_raw_rep1, "last")  # 177 proteins removed
ptetPCP_raw_rep2_SUPonly = fun_naFirstLast(ptetPCP_raw_rep2, "last")  # 121 proteins removed
ptetPCP_raw_rep3_SUPonly = fun_naFirstLast(ptetPCP_raw_rep3, "last")  # 146 proteins removed
ptetPCP_raw_rep1 = ptetPCP_raw_rep1[-which(ptetPCP_raw_rep1$Accession %in% ptetPCP_raw_rep1_SUPonly$Accession),]
ptetPCP_raw_rep2 = ptetPCP_raw_rep2[-which(ptetPCP_raw_rep2$Accession %in% ptetPCP_raw_rep2_SUPonly$Accession),]
ptetPCP_raw_rep3 = ptetPCP_raw_rep3[-which(ptetPCP_raw_rep3$Accession %in% ptetPCP_raw_rep3_SUPonly$Accession),]

ptetPCP_raw_SUPonly = intersect(intersect(ptetPCP_raw_rep1_SUPonly$Accession, 
                                          ptetPCP_raw_rep2_SUPonly$Accession), 
                                          ptetPCP_raw_rep3_SUPonly$Accession)
# write.csv(ptetPCP_raw_SUPonly, file = "/ParameciumPCP/tables/ptetPCP_SUPonly.csv", quote = F, row.names = F)

### Write to CSV
write.csv(ptetPCP_raw_rep1, file = "/ParameciumPCP/tables/ptetPCP_raw_rep1.csv", quote = F, row.names = F)
write.csv(ptetPCP_raw_rep2, file = "/ParameciumPCP/tables/ptetPCP_raw_rep2.csv", quote = F, row.names = F)
write.csv(ptetPCP_raw_rep3, file = "/ParameciumPCP/tables/ptetPCP_raw_rep3.csv", quote = F, row.names = F)

### Convert to MSN + Remove Unimportant R Objects
ptetPCP_raw_rep1_msn = readMSnSet2(file = "/ParameciumPCP/tables/ptetPCP_raw_rep1.csv",
                                   ecol = c("MAC.1", "X300g.1", "X1K.1", "X3K.1", "X5K.1", "X9K.1", 
                                            "X12K.1", "X15K.1", "X30K.1", "X79K.1", "X120K.1", "SUP.1"), 
                                   fnames = "Accession")
ptetPCP_raw_rep2_msn = readMSnSet2(file = "/ParameciumPCP/tables/ptetPCP_raw_rep2.csv",
                                   ecol = c("MAC.2", "X300g.2", "X1K.2", "X3K.2", "X5K.2", "X9K.2", 
                                            "X12K.2", "X15K.2", "X30K.2", "X79K.2", "X120K.2", "SUP.2"), 
                                   fnames = "Accession")
ptetPCP_raw_rep3_msn = readMSnSet2(file = "/ParameciumPCP/tables/ptetPCP_raw_rep3.csv",
                                   ecol = c("MAC.3", "X300g.3", "X1K.3", "X3K.3", "X5K.3", "X9K.3", 
                                            "X12K.3", "X15K.3", "X30K.3", "X79K.3", "X120K.3", "SUP.3"), 
                                   fnames = "Accession")
# rm(list=setdiff(ls(), c("ptetPCP_raw_rep1_msn", "ptetPCP_raw_rep2_msn", "ptetPCP_raw_rep3_msn")))

##### Data Sets
### Remove NAs and Normalize
ptetPCP_raw_rep1_norm = normalize(ptetPCP_raw_rep1_msn, method = "sum")
ptetPCP_raw_rep2_norm = normalize(ptetPCP_raw_rep2_msn, method = "sum")
ptetPCP_raw_rep3_norm = normalize(ptetPCP_raw_rep3_msn, method = "sum")

### Impute NBAVG and Normalize
ptetPCP_raw_rep1_nbavg_norm = normalize(impute(ptetPCP_raw_rep1_msn, method = "nbavg"), method = "sum")
ptetPCP_raw_rep2_nbavg_norm = normalize(impute(ptetPCP_raw_rep2_msn, method = "nbavg"), method = "sum")
ptetPCP_raw_rep3_nbavg_norm = normalize(impute(ptetPCP_raw_rep3_msn, method = "nbavg"), method = "sum")

### Impute Zeros
ptetPCP_raw_rep1_nbavg_norm_zero = impute(ptetPCP_raw_rep1_nbavg_norm, method = "zero")
ptetPCP_raw_rep2_nbavg_norm_zero = impute(ptetPCP_raw_rep2_nbavg_norm, method = "zero")
ptetPCP_raw_rep3_bavg_norm_zero = impute(ptetPCP_raw_rep3_nbavg_norm, method = "zero")

### Concatenate above datasets
# normalized
tmp1a =  merge(exprs(ptetPCP_raw_rep1_norm), exprs(ptetPCP_raw_rep2_norm), by=0)
rownames(tmp1a) = tmp1a$Row.names
tmp1a$Row.names = NULL 
tmp1b = merge(tmp1a, exprs(ptetPCP_raw_rep3_norm), by=0)
rownames(tmp1b) = tmp1b$Row.names
tmp1c = fData(ptetPCP_raw_rep1_norm)[which((fData(ptetPCP_raw_rep1_norm)$'Accession') %in% (rownames(tmp1b))),]
tmp1d = fData(ptetPCP_raw_rep2_norm)[which((fData(ptetPCP_raw_rep2_norm)$'Accession') %in% (rownames(tmp1b))),]
tmp1e = fData(ptetPCP_raw_rep3_norm)[which((fData(ptetPCP_raw_rep3_norm)$'Accession') %in% (rownames(tmp1b))),]
tmp1f = intersect(intersect(tmp1c, tmp1d), tmp1e)
ptetPCP_raw_all_norm = merge(tmp1f, tmp1b, by=0)
ptetPCP_raw_all_norm$Row.names = NULL

# nbavg imputed and normalized
tmp2a =  merge(exprs(ptetPCP_raw_rep1_nbavg_norm), exprs(ptetPCP_raw_rep2_nbavg_norm), by=0)
rownames(tmp2a) = tmp2a$Row.names
tmp2a$Row.names = NULL 
tmp2b = merge(tmp2a, exprs(ptetPCP_raw_rep3_nbavg_norm), by=0)
rownames(tmp2b) = tmp2b$Row.names
tmp2c = fData(ptetPCP_raw_rep1_nbavg_norm)[which((fData(ptetPCP_raw_rep1_nbavg_norm)$'Accession') %in% (rownames(tmp2b))),]
tmp2d = fData(ptetPCP_raw_rep2_nbavg_norm)[which((fData(ptetPCP_raw_rep2_nbavg_norm)$'Accession') %in% (rownames(tmp2b))),]
tmp2e = fData(ptetPCP_raw_rep3_nbavg_norm)[which((fData(ptetPCP_raw_rep3_nbavg_norm)$'Accession') %in% (rownames(tmp2b))),]
tmp2f = intersect(intersect(tmp2c, tmp2d), tmp2e)
ptetPCP_raw_all_nbavg_norm = merge(tmp2f, tmp2b, by=0)
ptetPCP_raw_all_nbavg_norm$Row.names = NULL

# nbavg imputed and normalized and zero imputed
tmp3a =  merge(exprs(ptetPCP_raw_rep1_nbavg_norm_zero), exprs(ptetPCP_raw_rep2_nbavg_norm_zero), by=0)
rownames(tmp3a) = tmp3a$Row.names
tmp3a$Row.names = NULL 
tmp3b = merge(tmp3a, exprs(ptetPCP_raw_rep3_bavg_norm_zero), by=0)
rownames(tmp3b) = tmp3b$Row.names
tmp3c = fData(ptetPCP_raw_rep1_nbavg_norm_zero)[which((fData(ptetPCP_raw_rep1_nbavg_norm_zero)$'Accession') %in% (rownames(tmp3b))),]
tmp3d = fData(ptetPCP_raw_rep2_nbavg_norm_zero)[which((fData(ptetPCP_raw_rep2_nbavg_norm_zero)$'Accession') %in% (rownames(tmp3b))),]
tmp3e = fData(ptetPCP_raw_rep3_bavg_norm_zero)[which((fData(ptetPCP_raw_rep3_bavg_norm_zero)$'Accession') %in% (rownames(tmp3b))),]
tmp3f = intersect(intersect(tmp3c, tmp3d), tmp3e)
ptetPCP_raw_all_bavg_norm_zero = merge(tmp3f, tmp3b, by=0)
ptetPCP_raw_all_bavg_norm_zero$Row.names = NULL

### Remove Aberrations (Flat Lines caused by Imputation Issue; should have all values == 0.08333...)
ptetPCP_raw_all_bavg_norm_zero = ptetPCP_raw_all_bavg_norm_zero %>% filter(round(MAC.1, 4) != round(0.0833, 4)) %>%
  filter(round(X15K.1, 4) != round(0.0833, 4)) %>%  filter(round(SUP.1, 4) != round(0.0833, 4)) %>%
  filter(round(MAC.2, 4) != round(0.0833, 4)) %>% filter(round(X15K.2, 4) != round(0.0833, 4)) %>%
  filter(round(SUP.2, 4) != round(0.0833, 4)) %>% filter(round(MAC.3, 4) != round(0.0833, 4)) %>%
  filter(round(X15K.3, 4) != round(0.0833, 4)) %>% filter(round(SUP.3, 4) != round(0.0833, 4))

ptetPCP_raw_all_nbavg_norm = ptetPCP_raw_all_nbavg_norm %>% filter(round(MAC.1, 4) != round(0.0833, 4)) %>%
  filter(round(X15K.1, 4) != round(0.0833, 4)) %>%  filter(round(SUP.1, 4) != round(0.0833, 4)) %>%
  filter(round(MAC.2, 4) != round(0.0833, 4)) %>% filter(round(X15K.2, 4) != round(0.0833, 4)) %>%
  filter(round(SUP.2, 4) != round(0.0833, 4)) %>% filter(round(MAC.3, 4) != round(0.0833, 4)) %>%
  filter(round(X15K.3, 4) != round(0.0833, 4)) %>% filter(round(SUP.3, 4) != round(0.0833, 4))

##### 
###### Data Sets:
# Organellar Classification:    ptetPCP_raw_all_bavg_norm_zero, ptetPCP_raw_rep1_nbavg_norm_zero, ptetPCP_raw_rep2_nbavg_norm_zero, ptetPCP_raw_rep3_nbavg_norm_zero

### Write to new CSV
write.csv(ptetPCP_raw_all_bavg_norm_zero, file = "/ParameciumPCP/tables/ptetPCP_nbavg_norm_zero.csv", quote = F, row.names = F)
write.csv(ptetPCP_raw_all_nbavg_norm, file     = "/ParameciumPCP/tables/ptetPCP_nbavg_norm.csv", quote = F, row.names = F)
write.csv(ptetPCP_raw_all_norm, file           = "/ParameciumPCP/tables/ptetPCP_norm.csv", quote = F, row.names = F)
write.csv(data.frame(Accession = ptetPCP_raw_all_norm$Accession, NAs = complete.cases(ptetPCP_raw_all_norm)), 
          "/ParameciumPCP/tables/ptetPCP-MissingValues.csv", quote = F, row.names = F)  # if protein has missing values

### Read in new MSNs
ptetPCP_raw_all_norm_NA = filterNA(readMSnSet2("/ParameciumPCP/tables/ptetPCP_norm.csv",
                                               ecol = c("MAC.1", "X300g.1", "X1K.1", "X3K.1", "X5K.1", "X9K.1", 
                                                        "X12K.1", "X15K.1", "X30K.1", "X79K.1", "X120K.1", "SUP.1", 
                                                        "MAC.2", "X300g.2", "X1K.2", "X3K.2", "X5K.2", "X9K.2", 
                                                        "X12K.2", "X15K.2", "X30K.2", "X79K.2", "X120K.2", "SUP.2",
                                                        "MAC.3", "X300g.3", "X1K.3", "X3K.3", "X5K.3", "X9K.3", 
                                                        "X12K.3", "X15K.3", "X30K.3", "X79K.3", "X120K.3", "SUP.3"),                                           
                                               fnames = "Accession"))
plotNA(ptetPCP_raw_all_norm_NA)  # read in without filterNA for this. save

ptetPCP_raw_all_nbavg_norm_NA = filterNA(readMSnSet2("/ParameciumPCP/tables/ptetPCP_nbavg_norm.csv",
                                                     ecol = c("MAC.1", "X300g.1", "X1K.1", "X3K.1", "X5K.1", "X9K.1", 
                                                              "X12K.1", "X15K.1", "X30K.1", "X79K.1", "X120K.1", "SUP.1", 
                                                              "MAC.2", "X300g.2", "X1K.2", "X3K.2", "X5K.2", "X9K.2", 
                                                              "X12K.2", "X15K.2", "X30K.2", "X79K.2", "X120K.2", "SUP.2",
                                                              "MAC.3", "X300g.3", "X1K.3", "X3K.3", "X5K.3", "X9K.3", 
                                                              "X12K.3", "X15K.3", "X30K.3", "X79K.3", "X120K.3", "SUP.3"),                                           
                                                     fnames = "Accession"))

ptetPCP_raw_all_nbavg_norm_zero_NA = filterNA(readMSnSet2("/ParameciumPCP/tables/ptetPCP_nbavg_norm_zero.csv",
                                                          ecol = c("MAC.1", "X300g.1", "X1K.1", "X3K.1", "X5K.1", "X9K.1", 
                                                                   "X12K.1", "X15K.1", "X30K.1", "X79K.1", "X120K.1", "SUP.1", 
                                                                   "MAC.2", "X300g.2", "X1K.2", "X3K.2", "X5K.2", "X9K.2", 
                                                                   "X12K.2", "X15K.2", "X30K.2", "X79K.2", "X120K.2", "SUP.2",
                                                                   "MAC.3", "X300g.3", "X1K.3", "X3K.3", "X5K.3", "X9K.3", 
                                                                   "X12K.3", "X15K.3", "X30K.3", "X79K.3", "X120K.3", "SUP.3"),                                           
                                                          fnames = "Accession"))

### Summary
nrow(fData(ptetPCP_raw_all_norm_NA))             # 2345 proteins
nrow(fData(ptetPCP_raw_all_nbavg_norm_NA))       # 4353 proteins
nrow(fData(ptetPCP_raw_all_nbavg_norm_zero_NA))  # 9026 proteins

### Clean up
rm(list=setdiff(ls(), 
                c("ptetPCP_raw_rep1_msn", "ptetPCP_raw_rep2_msn", "ptetPCP_raw_rep3_msn", 
                  "ptetPCP_raw_rep1_nbavg_norm", "ptetPCP_raw_rep2_nbavg_norm", "ptetPCP_raw_rep3_nbavg_norm",
                  "ptetPCP_raw_rep1_nbavg_norm_zero", "ptetPCP_raw_rep2_nbavg_norm_zero", "ptetPCP_raw_rep3_nbavg_norm_zero",
                  "ptetPCP_raw_all_nbavg_norm_zero_NA", 
                  "ptetPCP_raw_all_norm_NA", 
                  "ptetPCP_raw_all_nbavg_norm_NA")))

#save.image("/ParameciumPCP/RData/ptetPCP_MSN.RData")
load("/ParameciumPCP/RData/ptetPCP_MSN.RData")

#-------------------------------------------------------------------------------------------------------------------------
##### Visualize Data
# Correlation matrix
library("corrplot")
corrplot(cor(exprs(filterNA(ptetPCP_raw_all_norm_NA))), type = "lower", order = "hclust", 
         tl.col = "black", tl.srt = 45)
corrplot(cor(exprs(filterNA(ptetPCP_raw_all_nbavg_norm_NA))), type = "lower", order = "hclust", 
         tl.col = "black", tl.srt = 45)
corrplot(cor(exprs(filterNA(ptetPCP_raw_all_nbavg_norm_zero_NA))), type = "lower", order = "hclust", 
         tl.col = "black", tl.srt = 45)

# tSNE
set.seed(42)
plot2D(ptetPCP_raw_all_norm_NA, fcol = NULL, method = "t-SNE", col = "black")
set.seed(42)
plot2D(ptetPCP_raw_all_nbavg_norm_NA, fcol = NULL, method = "t-SNE", col = "black")
set.seed(42)
plot2D(ptetPCP_raw_all_nbavg_norm_zero_NA, fcol = NULL, method = "t-SNE", col = "black")

# PCs
plot2D(ptetPCP_raw_all_norm_NA, fcol = NULL, col = "black")
plot2D(ptetPCP_raw_all_norm_NA, fcol = NULL, col = "black", dims = c(1,3))
plot2D(ptetPCP_raw_all_norm_NA, fcol = NULL, col = "black", dims = c(1,4))
plot2D(ptetPCP_raw_all_norm_NA, method = "scree")
plot2D(ptetPCP_raw_all_nbavg_norm_zero_NA, method = "scree")

# PCA
plot2D(ptetPCP_raw_all_norm_NA, fcol = NULL, col = "black")
plot2D(ptetPCP_raw_all_nbavg_norm_NA, fcol = NULL, col = "black")
plot2D(ptetPCP_raw_all_nbavg_norm_zero_NA, fcol = NULL, col = "black")

# venn Diagram of descriptions
library(VennDiagram)
t1 = fData(ptetPCP_raw_rep1_msn)[complete.cases(exprs(ptetPCP_raw_rep1_msn)),]
t2 = fData(ptetPCP_raw_rep2_msn)[complete.cases(exprs(ptetPCP_raw_rep2_msn)),]
t3 = fData(ptetPCP_raw_rep3_msn)[complete.cases(exprs(ptetPCP_raw_rep3_msn)),]
venn.diagram(
  x = list(t1$Accession, t2$Accession, t3$Accession),
  category.names = c("Experiment 1" , "Experiment 2 " , "Experiment 3"), print.mode = c("raw", "percent"),
  filename = '/ParameciumPCP/plots/vennDiagram_noImputation.png',
  output=TRUE, main = "Proteins Identfied in Each Experiment (Missing Values Removed)"
)
t1 = fData(ptetPCP_raw_rep1_nbavg_norm)[complete.cases(exprs(ptetPCP_raw_rep1_nbavg_norm)),]
t2 = fData(ptetPCP_raw_rep2_nbavg_norm)[complete.cases(exprs(ptetPCP_raw_rep2_nbavg_norm)),]
t3 = fData(ptetPCP_raw_rep3_nbavg_norm)[complete.cases(exprs(ptetPCP_raw_rep3_nbavg_norm)),]
venn.diagram(
  x = list(t1$Accession, t2$Accession, t3$Accession),
  category.names = c("Experiment 1" , "Experiment 2 " , "Experiment 3"), print.mode = c("raw", "percent"),
  filename = '/ParameciumPCP/plots/vennDiagram_nbavg.png',
  output=TRUE, main = "Proteins Identfied in Each Experiment (NBAVG Imputation)"
)

#-------------------------------------------------------------------------------------------------------------------------
###### Compare PSMs and Missing Values
psmSummary = cbind(as.data.frame(unclass(summary(fData(ptetPCP_raw_all_nbavg_norm_zero_NA)$'nPSM'))),
                   as.data.frame(unclass(summary(fData(ptetPCP_raw_all_nbavg_norm_NA)$'nPSM'))),
                   as.data.frame(unclass(summary(fData(ptetPCP_raw_all_norm_NA)$'nPSM'))))
colnames(psmSummary) = c("Hybrid Imputation", "NBAVG Imputation", "No Imputation")
write.csv(psmSummary, file = "/ParameciumPCP/tables/psmSummary.csv", quote = F)

ptetPCP = read.csv(file = "/ParameciumPCP/tables/ptetPCP_raw.csv", header=T) 
ptetPCP$count_na = rowSums(is.na(ptetPCP))
ggplot(ptetPCP, aes(x=X..PSMs, y=count_na)) + geom_point() + geom_smooth() +  
  scale_x_log10() + ylim(c(0,115)) + 
  xlab("Nmmber of Peptide Spectral Matches per Protein") + 
  ylab("Number of Missing Values (out of 108)")


#-------------------------------------------------------------------------------------------------------------------------
##### Attach Marker Proteins
markedNEWNEW = addMarkers(ptetPCP_raw_all_nbavg_norm_zero_NA, 
                          markers = "/ParameciumPCP/tables/markers/markersGen7Final.csv", 
                          mcol = "markers", fcol = "Accession", verbose = T)

set.seed(42)
plot2D(markedNEWNEW, fcol = "markers", method = "t-SNE")
#addLegend(markedNEWNEW, fcol = "markers", where = "topleft", bty = "n", cex = .6)

# Marker Protein Resolution
mrkHClust(markedNEWNEW, fcol = "markers")
qs = QSep(markedNEWNEW)
plot(qs)
levelPlot(qs)

#save.image("/ParameciumPCP/RData/ptetPCP_markers.RData")
load("/ParameciumPCP/RData/ptetPCP_markers.RData")

# Compare with T. gondii
library(pRolocdata)
data("Barylyuk2020ToxoLopit")
plot(QSep(Barylyuk2020ToxoLopit))

#-------------------------------------------------------------------------------------------------------------------------
##### plot marker distribution profiles
comps = unique(fData(markedNEWNEW)$'markers')
source("/ParameciumPCP/scripts/source_groupColors.R")

for(compp in comps){
  svg(file=paste("/ParameciumPCP/plots/markers/", compp, ".svg", sep=""))
  plotDist(markedNEWNEW[fData(markedNEWNEW)[fData(markedNEWNEW)$'markers' == compp,'Accession']], 
           pcol = source_groupColors[compp], ylim = c(0,1))
  graphics.off()
}

#-------------------------------------------------------------------------------------------------------------------------
##### Perform Supervised Classification
w = table(getMarkers(markedNEWNEW, verbose = TRUE))  # class weight
w <- 1/w[names(w) != "unknown"]
svmMarked  = svmOptimisation(filterNA(markedNEWNEW), fcol = "markers", 
                             times = 100, xval = 5, classWeights = w, verbose = FALSE)                        
svmresMarked = svmClassification(filterNA(markedNEWNEW), fcol = "markers", assessRes = svmMarked)  
processingData(svmresMarked)  
preds1 = getPredictions(svmresMarked, fcol = "svm", mcol = "markers")   
medSVM = median(fData(svmresMarked)$svm.scores)  #median svm score: .53
preds2 = getPredictions(svmresMarked, fcol = "svm", t = medSVM, mcol = "markers")  # 1st quartile 0.37
firstQuartSVM = summary(fData(svmresMarked)$svm.scores)[2]
preds3 = getPredictions(svmresMarked, fcol = "svm", t = firstQuartSVM, mcol = "markers") 
thirdQuartSVM = summary(fData(svmresMarked)$svm.scores)[5]
preds4 = getPredictions(svmresMarked, fcol = "svm", t = thirdQuartSVM, mcol = "markers") 

ts1 = tapply(fData(preds1)$'svm.scores', fData(preds1)$'svm', median) 
write.csv(round(data.frame(ts1), 2), "/ParameciumPCP/tables/svmScores.csv", row.names = T, quote = F)

# F1
getParams (svmMarked)
levelPlot(svmMarked)

# Summary Table
tmp1 = merge(data.frame(table(fData(preds2)$'markers')), data.frame(table(fData(preds2)$'svm')), by = "Var1")
preds2_compare = merge(tmp1, data.frame(table(fData(preds2)$'svm.pred')), by = 'Var1')
colnames(preds2_compare) = c("Organellar Compartment", "# Markers", "# Predicted", "# Predicted (post-filtering)")
preds2_compare = preds2_compare %>% 
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))

write.csv(preds2_compare, "/ParameciumPCP/tables/SVMpredictions_Summary.csv", row.names = F, quote = F)

#save.image("/ParameciumPCP/RData/ptetPCP_svm.RData")
load("/ParameciumPCP/RData/ptetPCP_svm.RData")

# tSNE
set.seed(42)
plot2D(preds2, fcol = "svm", method = "t-SNE")
#addLegend(preds2, fcol = "svm", where = "topleft", bty = "n", cex = .6)
set.seed(42)
plot2D(preds2, fcol = "svm.pred", method = "t-SNE")
#addLegend(preds2, fcol = "svm.pred", where = "topleft", bty = "n", cex = .6)
set.seed(42)
plot2D(preds3, fcol = "svm.pred", method = "t-SNE")
#addLegend(preds3, fcol = "svm.pred", where = "topright", bty = "n", cex = .6)

plot3D(preds2, fcol = "svm.pred")
rgl.snapshot('/ParameciumPCP/plots/PCA-3D-svmpred.png', fmt = 'png')

# #-------------------------------------------------------------------------------------------------------------------------
##### TAGM Classification
# # TAGM-MAP
# load("/ParameciumPCP/RData/ptetPCP_tagm-map.RData")
# set.seed(42)
# mappars <- tagmMapTrain(markedNEWNEW)
# tagm_map <- tagmPredict(markedNEWNEW, mappars)
# fData(tagm_map)
# mappars
# table(fData(tagm_map)$tagm.map.allocation)
# plotEllipse(tagm_map, mappars)

# # TAGM-MCMC
# library(coda)
# load("/ParameciumPCP/RData/tagm-mcmc.RData")
# # p = tagmMcmcTrain(markedNEWNEW, numIter = 20000, burnin = 10000, thin = 20, numChains = 6)
# # p = tagmMcmcTrain(object = markedNEWNEW, fcol = "markers", method = "MCMC", numIter = 3000L, burnin = 300L, thin = 15L, mu0 = NULL, lambda0 = 0.01, nu0 = NULL, S0 = NULL, beta0 = NULL, u = 2, v = 10, numChains = 5L, , BPPARAM = BiocParallel::bpparam())
# # p = tagmMcmcTrain(object = markedNEWNEW, fcol = "markers", method = "MCMC", numIter = 2000L, burnin = 200L, thin = 10L, mu0 = NULL, lambda0 = 0.01, nu0 = NULL, S0 = NULL, beta0 = NULL, u = 2, v = 10, numChains = 4L, BPPARAM = BiocParallel::SerialParam())
# # p = tagmMcmcTrain(object = markedNEWNEW, fcol = "markers", method = "MCMC", numIter = 1000L, burnin = 100L, thin = 5L, mu0 = NULL, lambda0 = 0.01, nu0 = NULL, S0 = NULL, beta0 = NULL, u = 2, v = 10, numChains = 4L, BPPARAM = BiocParallel::bpparam())
# # save.image("tagm-mcmc.RData")
# p
# out <- mcmc_get_outliers(p)
# 
# lapply(out, summary)  # Chains: 1) 1679, 2) 1575, 3) 1565, 4) 1636, 5) 1677
# gelman.diag(out, transform = FALSE)  # 3.93 == BAD
# for (i in seq_len(length(p)))
#   plot(out[[i]], main = paste("Chain", i), auto.layout = FALSE, col = i)  # Chain 2 is awful. 3 isn't great. maybe 1 and 4 are ok?
# 
# gelman.diag(out[c(1,4,5)], transform = FALSE)  # much better but still bad
# 
# meanAlloc <- mcmc_get_meanComponent(p)
# for (i in seq_len(length(p)))
#   plot(meanAlloc[[i]], main = paste("Chain", i), auto.layout = FALSE, col = i)
# 
# meanoutProb <- mcmc_get_meanoutliersProb(p)
# gelman.diag(meanoutProb[c(1,4,5)])
# 
# geweke_test(out)  # chain's 2 and 3 haven't converged according to this test
# 
# out_converged <- p[-c(1,4,5)]
# 
# e14Tagm_converged_thinned <- mcmc_thin_chains(out_converged, freq = 4)
# e14Tagm_converged_pooled <- mcmc_pool_chains(out_converged)
# p = tagmMcmcProcess(e14Tagm_converged_pooled)
# mcmcmc <- tagmPredict(object = markedNEWNEW,
#                         params = p,
#                         probJoint = TRUE)
# data.frame(table(fData(mcmcmc)$'tagm.mcmc.allocation'))
# 
# set.seed(42)
# plot2D(mcmcmc, fcol = "tagm.mcmc.allocation",
#        cex = fData(mcmcmc)$tagm.mcmc.probability,
#        main = "TAGM MCMC allocations", method = "t-SNE")
# 
# # Summary
# data.frame(table(as.character(fData(mcmcmc)$'svm'), fData(mcmcmc)$'tagm.mcmc.allocation'))
# 
# # Visualize
# prot = "PTET.51.1.P0510184"
# plot(p, prot)
# plotDist(mcmcmc[prot], pcol = "red")
# fData(mcmcmc[prot])

#-------------------------------------------------------------------------------------------------------------------------
##### Perform Unsupervised Clustering (Kmeans)
set.seed(42)
preds2_km = MLearn( ~ ., filterNA(preds2), kmeansI, centers=17)

preds2_km = as.data.frame(unlist(preds2_km@RObject$cluster))

colnames(preds2_km) = "Cluster"
write.csv(preds2_km, "/ParameciumPCP/tables/Kmeans/ptetPCP_km17.csv", quote = F, row.names = T) 

preds2 = addMarkers(preds2, markers = "/ParameciumPCP/tables/Kmeans/ptetPCP_km17.csv", 
                          mcol = "Cluster", fcol = "Accession", verbose = T)
#save.image("/ParameciumPCP/RData/ptetPCP_superANDunsuper.RData")
load("/ParameciumPCP/RData/ptetPCP_superANDunsuper.RData")

# Optimal K?
# set.seed(42)
# k.max <- 100
# wss <- sapply(1:k.max, function(k){kmeans(exprs(preds2), k, nstart=50,iter.max = 15 )$tot.withinss})
# wss
# plot(1:k.max, wss, type="b", pch = 19, frame = FALSE)  # meh

# tSNE
library("randomcoloR")
set.seed(42)
knnColors = distinctColorPalette(17)

set.seed(42)
plot2D(preds2, fcol = "Cluster", method = "t-SNE", col = knnColors)
#addLegend(preds2, fcol = "Cluster", where = "topright", bty = "n", cex = .7)

#####
##### Kmeans vs SVM
knnsvm = fData(preds2)
knnsvm$Cluster = as.character(knnsvm$Cluster)
knnsvm$Cluster = factor(knnsvm$Cluster, levels=c("12", 
                                                 "2", 
                                                 "13", 
                                                 "1", 
                                                 "11", 
                                                 "15", 
                                                 "9", 
                                                 "7", 
                                                 "14", 
                                                 "5",
                                                 "16", 
                                                 "3", 
                                                 "17", 
                                                 "8", 
                                                 "6", 
                                                 "4", 
                                                 "10"))

    # manually tried to make a diagonal line
ggplot(knnsvm) + geom_count(aes(x=Cluster, y=svm.pred, size = ..n..)) + 
  guides(color="legend") + theme_bw() + theme(legend.position="none") + 
  xlab("De Novo Cluster") + ylab("Predicted Organellar Compartment")

#-------------------------------------------------------------------------------------------------------------------------
##### Write FASTA
library(seqinr)
preds2_mod = fData(preds2)[grep(pattern = "PTET.51.*", x = fData(preds2)$'Accession'),]
ptProts = read.fasta("/ParameciumPCP/sequences/ptet-proteome-UtoC.fa", seqtype = "AA", as.string = T)


write.fasta(sequences = getSequence(ptProts[preds2_mod$Accession], as.string = T), ### all protein sequences
            names = getName(ptProts[preds2_mod$Accession]), 
            file.out = "/ParameciumPCP/sequences/ptetPCP_allPredictions.fasta")

comps = unique(preds2_mod$svm.pred)
for(compp in comps){ # one fasta per compartment- median filtered (Used for MEME)
  write.fasta(sequences = getSequence(ptProts[preds2_mod[which(preds2_mod$svm.pred == compp),'Accession']], as.string = T), 
              names = getName(ptProts[preds2_mod[which(preds2_mod$svm.pred == compp),'Accession']]), 
              file.out = paste("/ParameciumPCP/sequences/", compp, "_predictions.fasta", sep=""))}

ptet_promoters = read.fasta("/ptetPCP/ptetPCP_fasta/ptet_allPromoters.fasta", as.string = T)
for(compp2 in comps){
  write.fasta(sequences = getSequence(ptet_promoters[which(names(ptet_promoters) %in% preds2_mod[which(preds2_mod$svm.pred == compp2),'Accession'])], as.string = T), 
              names = getName(ptet_promoters[which(names(ptet_promoters) %in% preds2_mod[which(preds2_mod$svm.pred == compp2),'Accession'])]), 
              file.out = paste("/ParameciumPCP/sequences/", compp2, "_predictions_upstream.fasta", sep=""))}

# Do this after getting properties
write.fasta(sequences = getSequence(ptProts[svmProps[which(svmProps$SignalPeptide), 'Accession']], as.string = T), 
            names = getName(ptProts[svmProps[which(svmProps$SignalPeptide), 'Accession']]), 
            file.out = "./ptetPCP_fasta/SignalPeptide.faa")
write.fasta(sequences = getSequence(ptProts[svmProps[which(svmProps$TargetP == "mTP"), 'Accession']], as.string = T), 
            names = getName(ptProts[svmProps[which(svmProps$TargetP == "mTP"), 'Accession']]), 
            file.out = "./ptetPCP_fasta/mitoTarget.faa")

# and MAC/Sup only
macOnly = read.csv("/ParameciumPCP/tables/ptetPCP_MAConly.csv", header=T)
write.fasta(sequences = getSequence(ptProts[macOnly$x], as.string = T), 
            names = getName(ptProts[macOnly$x]), 
            file.out = "/ParameciumPCP/sequences/ptetPCP_MAConly.fasta")

supOnly = read.csv("/ParameciumPCP/tables/ptetPCP_SUPonly.csv", header=T)
write.fasta(sequences = getSequence(ptProts[supOnly$x], as.string = T), 
            names = getName(ptProts[supOnly$x]), 
            file.out = "/ParameciumPCP/sequences/ptetPCP_SUPonly.fasta")

# no imputation
write.fasta(sequences = getSequence(ptProts[rownames(fData(ptetPCP_raw_all_norm_NA))[grep(pattern = "PTET*", rownames(fData(ptetPCP_raw_all_norm_NA)))]],
                                    as.string = T), 
            names = getName(ptProts[rownames(fData(ptetPCP_raw_all_norm_NA))[grep(pattern = "PTET*", rownames(fData(ptetPCP_raw_all_norm_NA)))]]), 
            file.out = "/ParameciumPCP/sequences/ptetPCP_noImpute.fasta")

#-------------------------------------------------------------------------------------------------------------------------
##### Properties
write.csv(fData(preds2), "/ParameciumPCP/tables/properties/ptetPCP_FData.csv", quote = F, row.names = F)

  # Proteins/sequences were plugged into the following software:
    # ParameciumDB Biomart (now Sherlock)
    # TargetP
    # SignalP
    # eggNOG
    # ghostKOALA
    # NLStradamus

svmProps = read.csv("/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", header=T)
targetP = read.csv("/ParameciumPCP/motifs/TargetP/ptetPCP_TargetP.csv", header=T)
svmProps = merge(svmProps, targetP[,c('Accession', 'Prediction_TargetP')], by='Accession')
#write.csv(svmProps, file = "/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", quote=F, row.names = F)
koala = read.csv("/ParameciumPCP/tables/properties/ghostKOALA.csv", header=T)[,c(1:3)]
write.csv(merge(svmProps, koala, by.x="Accession", by.y="Query"), "/ParameciumPCP/tables/properties/ptetPCP_GHprops.csv", quote = F, row.names = F)
svmProps$svm.pred = gsub(pattern = " ", replacement = "\n", svmProps$svm.pred)


# Quick summary of all  features
allSummaries = by(svmProps[grep(pattern = "PTET*", x = svmProps$Accession),], factor(svmProps[grep(pattern = "PTET*", x = svmProps$Accession),'svm.pred']), summary)
allSummariesDF = do.call(rbind.data.frame, allSummaries)
write.csv(allSummariesDF[grep("Median*", allSummariesDF$Freq),], 
          file = "/ParameciumPCP/tables/properties/ptetPCP_Properties-Summary-Median.csv", quote = F)
write.csv(merge(merge(data.frame(svmProps[grep(pattern = "PTET*", x = svmProps$Accession),] %>% group_by(svm.pred) %>% summarise(TMdomain = sum(TMD==T)/(sum(TMD==T) + sum(TMD==F)))),
                      data.frame(svmProps[grep(pattern = "PTET*", x = svmProps$Accession),] %>% group_by(svm.pred) %>% summarise(SignalP = sum(SP==TRUE)/(sum(SP==TRUE) + sum(SP==FALSE)))), by="svm.pred"), 
                data.frame(svmProps[grep(pattern = "PTET*", x = svmProps$Accession),] %>% group_by(svm.pred) %>% summarise(mTP = sum(TargetP=="mTP")/(sum(TargetP=="mTP") + sum(TargetP=="SP") + sum(TargetP=="OTHER")))), by="svm.pred"), 
          file = "/ParameciumPCP/tables/properties/ptetPCP_Properties-Summary-TrueFalse.csv", quote = F, row.names = F)

##### Enrichment of Features
#  Genomic features
# ggplot(svmProps) + geom_count(aes(x=svm.pred, y=Chromosome))
# sort(table(svmProps$Chromosome), decreasing = T)
#   # Uninteresting

ggplot(svmProps) + geom_violin(aes(x=svm.pred, y=mRNA)) +   scale_y_continuous(trans = "log10") + 
  xlab("Predicted Compartment") + ylab("Log mRNA Expression (VEG Growth)")
summary(aov(mRNA ~ svm.pred, data = svmProps))
TukeyHSD(aov(mRNA ~ svm.pred, data = svmProps))
write.csv(data.frame(TukeyHSD(aov(mRNA ~ svm.pred, svmProps))$'svm.pred')[grep("unknown", rownames(data.frame(TukeyHSD(aov(mRNA ~ svm.pred, svmProps))$'svm.pred'))),], 
          file = "/ParameciumPCP/tables/properties/ptetPCP_mRNA_Tukey.csv", quote = F, row.names = T)  # csv for all compartment vs unknown comparisons

# ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=nExons)) + ylim(c(0,10)) + 
#   xlab("Predicted Compartment") + ylab("Exons/gene")
# # meh
# 
# ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=nIntrons)) + ylim(c(0,10)) + 
#   xlab("Predicted Compartment") + ylab("Introns/gene")
# # meh

svmProps = read.csv("/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", header=T)
#svmProps$svm.pred = gsub(pattern = " ", replacement = "\n", svmProps$svm.pred)

tmpA = data.frame(svmProps[which(complete.cases(svmProps[,'DE_Reciliation'])),] %>% group_by(svm.pred) %>% summarise(DErecilliation = sum(DE_Reciliation==TRUE)/(sum(DE_Reciliation==TRUE) + sum(DE_Reciliation==FALSE))))
tmpA$DErecilliation = as.numeric(tmpA$DErecilliation)*100
ggplot(tmpA,aes(x=svm.pred,y=DErecilliation)) + geom_bar(stat='identity') + 
  xlab("Organellar Compartment") + ylab("% of Genes DE after Recilliation")
write.csv(data.frame(svmProps[which(complete.cases(svmProps[,'DE_Reciliation'])),] %>% group_by(svm.pred) %>% summarise(DErecilliation = sum(DE_Reciliation==TRUE)/(sum(DE_Reciliation==TRUE) + sum(DE_Reciliation==FALSE)))),
          file="/ParameciumPCP/tables/properties/ptetPCP_Properties_DE-recilia.csv", quote = F, row.names = F)

tmpB = data.frame(svmProps[which(complete.cases(svmProps[,'DE_TrichocystDischarge'])),] %>% group_by(svm.pred) %>% summarise(DE_TrichocystDischarge = sum(DE_TrichocystDischarge==TRUE)/(sum(DE_TrichocystDischarge==TRUE) + sum(DE_TrichocystDischarge==FALSE))))
tmpB$DE_TrichocystDischarge = as.numeric(tmpB$DE_TrichocystDischarge)*100
ggplot(tmpB,aes(x=svm.pred,y=DE_TrichocystDischarge)) + geom_bar(stat='identity') + 
  xlab("Organellar Compartment") + ylab("% of Genes DE after Trichocyst Discharge")

ggplot(svmProps) + geom_count(aes(x=svm.pred, y=DE_TrichocystDischarge))
write.csv(data.frame(svmProps[which(complete.cases(svmProps[,'DE_TrichocystDischarge'])),] %>% group_by(svm.pred) %>% summarise(DEtrich = sum(DE_TrichocystDischarge==TRUE)/(sum(DE_TrichocystDischarge==TRUE) + sum(DE_TrichocystDischarge==FALSE)))),
          file="/ParameciumPCP/tables/properties/ptetPCP_Properties_DE-trich.csv", quote = F, row.names = F)
# ~2.4% of unknown: 6% of Trichocyst Matrix

# ggplot(svmProps) + geom_count(aes(x=svm.pred, y=DE_Autogamy))
# data.frame(svmProps[which(complete.cases(svmProps[,'DE_Autogamy'])),] %>% group_by(svm.pred) %>% summarise(DE_Auto = sum(DE_Autogamy==TRUE)/(sum(DE_Autogamy==TRUE) + sum(DE_Autogamy==FALSE))))
# # meh
# 
# ggplot(svmProps) + geom_count(aes(x=svm.pred, y=ExpressionGroup))
# data.frame(svmProps[which(complete.cases(svmProps[,'ExpressionGroup'])),] %>% group_by(svm.pred) %>% summarise(DE_Auto = sum(DE_Autogamy==TRUE)/(sum(DE_Autogamy==TRUE) + sum(DE_Autogamy==FALSE))))
# # meh... but
# write.csv(table(svmProps$svm.pred, svmProps$ExpressionGroup), file = "/ParameciumPCP/tables/properties/ptetPCP_expressionGroupSummary.csv")
# # still meh

#  Proteomic features
ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=nAA))  + ylim(0,3500) +
  xlab("Predicted Compartment")  + ylab("Protein Length (N amino acids)")
summary(aov(nAA ~ svm.pred, data = svmProps))
TukeyHSD(aov(nAA ~ svm.pred, data = svmProps))
write.csv(data.frame(TukeyHSD(aov(nAA ~ svm.pred, data = svmProps))$'svm.pred'), 
          file = "/ParameciumPCP/tables/properties/ptetPCP_nAA.csv", quote = F, row.names = T)

ggplot(svmProps) + geom_boxplot(aes(x=svm.pred, y=pI))  + ylim(0,12) +
  xlab("Predicted Compartment")  + ylab("Isoelectric Point")
summary(aov(pI ~ svm.pred, svmProps))
TukeyHSD(aov(pI ~ svm.pred, svmProps))
write.csv(data.frame(TukeyHSD(aov(pI ~ svm.pred, svmProps))$'svm.pred'), 
          file = "/ParameciumPCP/tables/properties/ptetPCP_pI_Tukey.csv", quote = F, row.names = T)

ggplot(svmProps) + geom_violin(aes(x=svm.pred, y=nPSM)) + scale_y_continuous(trans = "log") + 
  xlab("Predicted Compartment")  + ylab("Log PSMs/protein")
summary(aov(nPSM ~ svm.pred, svmProps))
TukeyHSD(aov(nPSM ~ svm.pred, svmProps))
write.csv(data.frame(TukeyHSD(aov(nPSM ~ svm.pred, svmProps))$'svm.pred'), 
          file = "/ParameciumPCP/tables/properties/ptetPCP_nPSM_Tukey.csv", quote = F, row.names = T)

tmpC = data.frame(svmProps[which(complete.cases(svmProps[,'TMD'])),] %>% group_by(svm.pred) %>% summarise(TransmembraneHelix = sum(TMD==TRUE)/(sum(TMD==TRUE) + sum(TMD==FALSE))))
tmpC$TransmembraneHelix = as.numeric(tmpC$TransmembraneHelix)*100
ggplot(tmpC,aes(x=svm.pred,y=TransmembraneHelix)) + geom_bar(stat='identity') + 
  xlab("Organellar Compartment") + ylab("% of Proteins with Predicted Transmembrane Domain")

ggplot(svmProps[which(!is.na(svmProps$TMD)),]) + geom_count(aes(x=svm.pred, y=TMD, color = ..n..)) + scale_size(range=c(0,10)) + 
  xlab("Predicted Compartment") + ylab("Predicted Transmembrane Domain")
table(svmProps$TMD)

tmpD = data.frame(svmProps[which(complete.cases(svmProps[,'SP'])),] %>% group_by(svm.pred) %>% summarise(SignalPeptide = sum(SP==TRUE)/(sum(SP==TRUE) + sum(SP==FALSE))))
tmpD$SignalPeptide = as.numeric(tmpD$SignalPeptide)*100
ggplot(tmpD,aes(x=svm.pred,y=SignalPeptide)) + geom_bar(stat='identity') + 
  xlab("Organellar Compartment") + ylab("% of Proteins with Predicted Signal Peptide")

ggplot(svmProps[which(!is.na(svmProps$SP)),]) + geom_count(aes(x=svm.pred, y=SP, color = ..n..)) + scale_size(range=c(0,10)) + 
  xlab("Predicted Compartment") + ylab("Predicted Signal Peptide")
table(svmProps$SP)[1] / (table(svmProps$SP)[1] + table(svmProps$SP)[2])

tmpE = data.frame(svmProps[which(complete.cases(svmProps[,'TargetP'])),] %>% group_by(svm.pred) %>% summarise(TargetP = sum(TargetP=='mTP')/(sum(TargetP=='mTP') + sum(TargetP=='OTHER'))))
tmpE$TargetP = as.numeric(tmpE$TargetP)*100
ggplot(tmpE,aes(x=svm.pred,y=TargetP)) + geom_bar(stat='identity') + 
  xlab("Organellar Compartment") + ylab("% of Proteins with Predicted Mitochondrial Targeting Sequence")

ggplot(svmProps[grep("PTET*", svmProps$Accession),]) + geom_count(aes(x=svm.pred, y=TargetP, color = ..n..)) + 
  scale_size(range=c(0,10)) + 
  xlab("Predicted Compartment") + ylab("Signal Peptide, Mitochondrial Targetting Sequence, or Neither")

# stats- chi^2 tests
source("/ParameciumPCP/scripts/fun_chiSquare_compartments.R")
tmp = svmProps[, c("svm.pred", "DE_TrichocystDischarge", "DE_Reciliation", "TMD", "SP", "TargetP")]
tmp$DE_TrichocystDischarge = ifelse(tmp$DE_TrichocystDischarge == TRUE, 1, 0)
tmp$DE_Reciliation = ifelse(tmp$DE_Reciliation == TRUE, 1, 0)
tmp$TMD = ifelse(tmp$TMD == TRUE, 1, 0)
tmp$SP = ifelse(tmp$SP == TRUE, 1, 0)
tmp$TargetP = ifelse(tmp$TargetP == 'mTP', 1, 0)

# Reciliation
lCilia = fun_chiSquare_compartments(tmp, "DE_Reciliation")
which(lCilia < 0.0025)  
    # p > 0.0025: Membrane Trafficking Soluble, Axoneme(relevant)
# Trich Discharge
lTrich = fun_chiSquare_compartments(tmp, "DE_TrichocystDischarge")
which(lTrich < 0.0025)  
  # p > 0.0025: ER
# TMD
lTMD = fun_chiSquare_compartments(tmp, "TMD")
which(lTMD < 0.0025)  
  # p > 0.0025: Both MTs (insoluble is relevant), Ribo, ER (Relevant), Nuc Sol, BB core, Cyto, TM, Axo, Protea
# SP
lSP = fun_chiSquare_compartments(tmp, "SP")
which(lSP < 0.0025)  
   # p > 0.0025: Mito, lyso (r), MT sol, Ribo, ER (r), nuc sol, BB core, cyto, TM (r), SA (r)
# mTP
lMTP = fun_chiSquare_compartments(tmp, "TargetP")
which(lMTP < 0.0025)  
 # p > 0.0025: Mito (r),


#-------------------------------------------------------------------------------------------------------------------------
##### Ribosomes
ribos = svmProps[grep(pattern = "Ribosomal protein*", x = svmProps$Description),]
table(ribos$svm.pred)
nucleoli = ribos[ribos$svm.pred == "Nuclei Soluble" | ribos$svm.pred == "Nuclei Insoluble",]
mitRibos = ribos[ribos$svm.pred == "Mitochondria",]
cytoRibos = ribos[ribos$svm.pred == "Ribosome",]

#-------------------------------------------------------------------------------------------------------------------------
# Peroxisome
faa = c("PTET.51.1.P1040049", "PTET.51.1.P1310070", "PTET.51.1.P1040088")
pxa = "PTET.51.1.P1840030"
racemase = "PTET.51.1.P0050127"
reductase = "PTET.51.1.P1590076"
lyase = "PTET.51.1.P0270267"


#-------------------------------------------------------------------------------------------------------------------------
##### NUPs
nup98 = c("PTET.51.1.P0010578", "PTET.51.1.P0950099", "PTET.51.1.P0750016")
nup98 = c("PTET.51.1.P0010578", "PTET.51.1.P0750016")

plotDist(preds2[nup98], pcol = "red", ylim = c(0,1))
fData(preds2[nup98])

nup210 = c("PTET.51.1.P1240016", "PTET.51.1.P0880039")

plotDist(preds2[nup210], pcol = "red", ylim = c(0,1))
fData(preds2[nup210])

plotDist(preds2[c(nup210, nup98)], pcol = "red", ylim = c(0,1))

nup37 = "PTET.51.1.P1740028"
nup42 = "PTET.51.1.P0780138"
nup43 = "PTET.51.1.P0110417"
plotDist(preds2[c(nup37, nup42, nup43)], pcol = "red", ylim = c(0,1))
fData(preds2[c(nup37, nup42, nup43)])

plotDist(preds2[c(nup98, nup210, nup37, nup42, nup43)], pcol = "red", ylim = c(0,1))
fData(preds2[c(nup98, nup210, nup37, nup42, nup43)])

#-------------------------------------------------------------------------------------------------------------------------
##### Glycolysis
phosphoglucoseIsomerase = c("PTET.51.1.P0300257")
phosphofructokinase = c("PTET.51.1.P0480063", "PTET.51.1.P0670026")
aldoase = c("PTET.51.1.P0770173", "PTET.51.1.P0940159", "PTET.51.1.P0810050", "PTET.51.1.P0900034")
triosephosphateIsomerase = c("PTET.51.1.P1550028", "PTET.51.1.P1070138", "PTET.51.1.P1550027", "PTET.51.1.P1370031")
glycerol3phosphateDehydrogenase = c("PTET.51.1.P0380195", "PTET.51.1.P0500184")
phosphoglycerateKinase = c("PTET.51.1.P0700046", "PTET.51.1.P1180061")
phosphoglyceromutase = c("PTET.51.1.P0890070", "PTET.51.1.P1190051")
enolase = c("PTET.51.1.P0870049", "PTET.51.1.P0100278", "PTET.51.1.P0590214", "PTET.51.1.P0040146")
pyruvateKinse = c("PTET.51.1.P0110153",  "PTET.51.1.P0360217", "PTET.51.1.P0070160", "PTET.51.1.P0210069", "PTET.51.1.P0100415")

for(glyco in list(phosphoglucoseIsomerase, phosphofructokinase, aldoase, triosephosphateIsomerase, glycerol3phosphateDehydrogenase, phosphoglycerateKinase, phosphoglyceromutase, enolase, pyruvateKinse)){
  svg(paste("/ParameciumPCP/plots/dist-Glycolysis-", glyco, ".svg"))
  plotDist(preds2[as.character(glyco)], pcol = "black", ylim = c(0,1))
  dev.off()
}

#-------------------------------------------------------------------------------------------------------------------------
##### Motifs
# Nuclear Localization Sequence
fun_nlsCalc = function(nls){length(grep("Posterior @ 0.6\t[0-9]", nls)) / length(grep("PTET.51.1.P*", nls))*100}
nls1 = read_lines(file = "/ParameciumPCP/test/testNLS1.txt")  # read_lines is a wimp
nls2 = read_lines(file = "/ParameciumPCP/test/testNLS2.txt")
nls = c(nls1, nls2)
nlsPresent = as.character(sapply(sapply(nls[grep("Viterbi Path\t[0-9]", nls)-1], strsplit, split = "Predictions for "), "[", 2))
nlsProts = data.frame(nlsPresent, rep(TRUE, length(nlsPresent)))
colnames(nlsProts) = c("Accession", "NLS")  
nlsAbsent = svmProps$Accession[which(!svmProps$Accession %in% nlsProts$Accession)]
nlsProts2 = data.frame(nlsAbsent, rep(FALSE, length(nlsAbsent)))
colnames(nlsProts2) = c("Accession", "NLS")  
nlsProts = rbind(nlsProts, nlsProts2)
nlsProts = merge(nlsProts, svmProps[,c('Accession', 'svm.pred')])
tmpf = nlsProts %>% group_by(svm.pred) %>% summarise(NLS = sum(NLS==TRUE)/(sum(NLS==TRUE) + sum(NLS==FALSE)))
tmpf$NLS = tmpf$NLS*100
write.csv(tmpf, "/ParameciumPCP/tables/properties/ptetPCPnls.csv", quote = F, row.names = F)

ggplot(tmpf,aes(x=svm.pred,y=NLS)) + geom_bar(stat='identity') + 
  xlab("Organellar Compartment") + ylab("% of Proteins with Predicted Nuclear Localization Signals")


# stats

lNLS = fun_chiSquare_compartments(nlsProts[,c(3,2)], "NLS")
which(lNLS < 0.0025)  
 # both nuclear (r), ribosome (r), rest of sig depleted

# # all interesting plots in one
# nlsProts$svm.pred = gsub(pattern = " ", replacement = "\n", nlsProts$svm.pred)
# tmpG = melt(as.data.frame(merge(merge(merge(merge(merge(tmpA, tmpB),tmpC),tmpD),tmpE),tmpf)))
# 
# ggplot(tmpG, aes(x=svm.pred,y=value,fill=variable)) + geom_bar(stat = "identity", position = 'dodge') +  
#   xlab("Organellar Compartment") + ylab("Percentage of Genes/Proteins with a Given Feature") + ylim(c(0,70)) + 
#   scale_fill_discrete(labels=c('DE- Reciliation', 'DE- Trichocyst Discharge', 'Predicted TMD', 'Predicted SP', 'Predicted mTS', 'predicted NLS'))
#   
# ggplot(tmpG, aes(x=svm.pred,y=value,fill=variable)) + geom_bar(stat = "identity", position = 'dodge') +  
#   xlab("Organellar Compartment") + ylab("Percentage of Genes/Proteins with a Given Feature") + ylim(c(0,70)) + 
#   scale_fill_discrete(labels=c('DE- Reciliation', 'DE- Trichocyst Discharge', 'Predicted TMD', 'Predicted SP', 'Predicted mTS', 'predicted NLS')) + 
#   geom_hline(yintercept=length(which(svmProps$DE_Reciliation == T))/ (length(which(svmProps$DE_Reciliation == T)) + length(which(svmProps$DE_Reciliation == F)))*100, color = "indianred", linetype = "dashed") + 
#   geom_hline(yintercept=length(which(svmProps$DE_TrichocystDischarge == T))/ (length(which(svmProps$DE_TrichocystDischarge == T)) + length(which(svmProps$DE_TrichocystDischarge == F)))*100, color = "darkgoldenrod", linetype = "dashed") + 
#   geom_hline(yintercept=length(which(svmProps$TMD == T))/ (length(which(svmProps$TMD == T)) + length(which(svmProps$TMD == F)))*100, color = "green3", linetype = "dashed") + 
#   geom_hline(yintercept=length(which(svmProps$SP == T))/ (length(which(svmProps$SP == T)) + length(which(svmProps$SP == F)))*100, color = "cyan3", linetype = "dashed") + 
#   geom_hline(yintercept=length(which(svmProps$TargetP == 'mTP'))/ (length(which(svmProps$TargetP == 'mTP')) + length(which(svmProps$TargetP == 'OTHER')))*100, color = "steelblue3", linetype = "dashed") + 
#   geom_hline(yintercept=length(which(nlsProts$NLS == T))/ (length(which(nlsProts$NLS == T)) + length(which(nlsProts$NLS == F)))*100, color = "magenta2", linetype = "dashed")


##### Signal Peptides
svmProps = read.csv("/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", header=T)
#svmProps$svm.pred = gsub(pattern = " ", replacement = "\n", svmProps$svm.pred)

spsMeme = read.table("/ParameciumPCP/motifs/SignalP/meme-mod.txt", header = F)
table(spsMeme$V2)

sps = read.table("/ParameciumPCP/motifs/SignalP/fimo-10e4.tsv", header = T)[,3:5]  # motif discovered via MEME and mapped via FIMP (e04)
spsSVM = merge(x = sps, svmProps[,c("Accession", "nAA", "svm.pred", "TMD", "pI")], by.x='sequence_name', by.y="Accession")
spsSVM$percentile = (spsSVM$start/spsSVM$nAA)*100
spsSVM$svm.pred = gsub(" ", "\n", spsSVM$svm.pred)
spsSVM = spsSVM[which(spsSVM$svm.pred %in% c('ER', 'Membrane\nTrafficking\nInsoluble', 'Membrane\nTrafficking\nSoluble', 'Surface\nAntigen', 'Trichocyst\nMatrix', 'Lysosome')),]
table(spsSVM$svm.pred)

# hydrophobicity
library(Peptides)
library(seqinr)
ptProts = read.fasta("/ParameciumPCP/sequences/ptet-proteome-UtoC.fa", seqtype = "AA", as.string = T)
hydros = lapply(getSequence(ptProts[spsSVM$sequence_name], as.string = T), hydrophobicity, scale = "KyteDoolittle")
spsSVM$hydrophobic = as.numeric(hydros)

#ggplot(spsSVM) + geom_point(aes(x=hydrophobic, y=percentile))
#ggplot(spsSVM) + geom_point(aes(x=svm.pred, y=hydrophobic))
#ggplot(spsSVM) + geom_point(aes(x=pI, y=percentile)) + scale_y_log10() + facet_wrap(~ svm.pred)

ggplot(spsSVM, aes(x=svm.pred, y=percentile)) + 
  geom_count(position=position_jitter(h=0.1, w=0.1), alpha = 0.75, shape = 21, size=2) + 
  theme(legend.position = "none") + theme_classic() + 
  xlab("Organellar Compartment") + ylab("Motif Positition (percentile; N to C)") 


# ggplot(melt(spsSVM[,c("svm.pred", "TMD", "percentile")]), aes(x=TMD, y=value)) +
#   geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3) + 
#   theme_classic() + facet_wrap(~ svm.pred) + 
#   xlab("Transmembrane Domain Presence") + ylab("Motif Positition (percentile; N to C)") 

# ggplot(melt(spsSVM[,c("svm.pred", "TMD", "hydrophobic")]), aes(x=TMD, y=value)) +
#   geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3) + 
#   theme_classic() + facet_wrap(~ svm.pred) + 
#   xlab("Transmembrane Domain Presence") + ylab("Hydrophobicity (Kyte-Doolittle Scale") 
# 

ggplot(spsSVM, aes(x=TMD, y=percentile)) + geom_violin() + geom_jitter() + facet_wrap(~ svm.pred) + 
  xlab("Transmembrane Domain Presence") + ylab("Motif Position (percentile: N to C") + theme_classic()
  
  
# TMD position (c-term?)
tmhmm = read.table("/ParameciumPCP/tables/properties/tmhmm.txt")[,c(1,2,5,6)]
tmhmm$V5 = as.numeric(unlist(lapply(sapply(tmhmm$V5, strsplit, "PredHel="), "[", 2)))
tmhmm$V6 = as.character(lapply(sapply(tmhmm$V6, strsplit, "Topology="), "[", 2))

tmhmm2 = merge(spsSVM, tmhmm, by.x="sequence_name", by.y="V1")
tmhmm2 = tmhmm2[,c('sequence_name', 'nAA', 'percentile', 'V6')]

tmhmm3 = tmhmm2[-which(tmhmm2$V6 == "o" | tmhmm2$V6 == "i"),]
tmhmm3$CtermTMD = NA

for(z in 1:nrow(tmhmm3)){  # gotta do it
  if(length(which(((as.numeric(unlist(regmatches(tmhmm3$V6[z], gregexpr("[[:digit:]]+", tmhmm3$V6[z])))) / as.numeric(tmhmm3$nAA[z]))*100) > 66.67))){
    tmhmm3$CtermTMD[z] = T
  } else{
    tmhmm3$CtermTMD[z] = F
  }
}

ggplot(tmhmm3, aes(x=CtermTMD, y=percentile)) + geom_count()

# HDEL or KDEL?
library(seqinr)
ptProts = read.fasta("/Paramecium_Proteomes-Predicted/ptetraurelia_mac_51_annotation_v2.0.protein.fa", seqtype = "AA", as.string = T)

hdel = getName(ptProts[grep("HDEL\\*", getSequence(ptProts, as.string = T))])
kdel = getName(ptProts[grep("KDEL\\*", getSequence(ptProts, as.string = T))])
del = getName(ptProts[grep("DEL\\*", getSequence(ptProts, as.string = T))])

fData(preds2)[which(fData(preds2)$'Accession' %in% hdel),]  # 3 ER
fData(preds2)[which(fData(preds2)$'Accession' %in% kdel),]  # 2 ER
fData(preds2)[which(fData(preds2)$'Accession' %in% del),]  # 6 ER + 1 nuclear

# sortilins
sorts = read.csv("/ParameciumPCP/tables/FOI/sortilins.csv", header=T)

plotDist(preds2[sorts$Accession[which(sorts$Accession %in% fData(preds2)$'Accession')]], pcol = "red")
fData(preds2[sorts$Accession[which(sorts$Accession %in% fData(preds2)$'Accession')]])$'svm.pred'
fData(preds2[sorts$Accession[which(sorts$Accession %in% fData(preds2)$'Accession')]])$'svm'

#-------------------------------------------------------------------------------------------------------------------------
##### mRNA vs PSMs
svmProps$mRNA_VEG_log = as.numeric(log(svmProps$mRNA))
svmProps$nPSM_log = as.numeric(log(svmProps$nPSM))
ggplot(svmProps, aes(x=mRNA_VEG, y=nPSM)) + geom_point() + geom_smooth(method = "lm")
ggplot(svmProps, aes(x=mRNA_VEG_log, nPSM_log)) + geom_point() + geom_smooth(method = "lm") + 
  xlab("Log mRNA Expression (VEG Growth)") + ylab("Log PSMs/protein")

summary(lm(nPSM ~ mRNA, svmProps))  # R^2 0.1098, p-value < 2.2e-16
summary(lm(log(nPSM) ~ log(mRNA), svmProps[,c("mRNA", "nPSM")]+0.00000001 ))  # R^2 0.3214, p-value < 2.2e-16
      # adding extremely small value to handle mRNA = 0

#-------------------------------------------------------------------------------------------------------------------------
##### Word Clouds
library("tm")
library('wordcloud')
library(tidyverse)
source("/ParameciumPCP/scripts/fun_makeWordCloud.R")
svmProps = read.csv("/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", header=T)

commonWords = c("coil", " coil", "coil ", 
                "coiled",  " coiled", "coiled ", "  coiled", "coiled  ", " coiled ", 
                "  coiled  ", "  coiled  ", "Coiled", "Coiled ", " Coiled",
                "fold", "subunit", "protein", "domain", "family", "function", "containing", "profile", "factor", 
                " protein", "protein ", "protein  ", "  protein", "Protein", "unknown", " unknown", "unknown ")
  # makeWordCloud() make word cloud given dataframe and vector of common words

for(orgs in unique(svmProps$svm.pred)){  # gene descriptions
  svg(file=paste("/ParameciumPCP/plots/wordcloud/description_", orgs, ".svg", sep=""))
  fun_makeWordCloud(svmProps[which(svmProps$svm.pred == orgs),], commonWords)
  graphics.off()
}


#-------------------------------------------------------------------------------------------------------------------------
#####
##### Population Genomics
svmProps = read.csv("/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", header=T)
svmProps$svm = gsub(pattern = " ", replacement = "\n", svmProps$svm)
svmProps$svm.pred = gsub(pattern = " ", replacement = "\n", svmProps$svm.pred)

ptPPI = read.table("/ParameciumPCP/tables/popgenomics/all_genes_filtered_tetraurelia.char", header=T)
ptetConversion = read.csv("/ParameciumPCP/tables/popgenomics/ptetIDconversion.csv", header=T)

ptPopGen = merge(merge(ptetConversion, ptPPI, by.x="GeneB", by.y="gene"), 
                 fData(preds2), by.x="GeneA", by.y="Accession")
ptPopGen = ptPopGen[!duplicated(ptPopGen), ]
ptPopGen$svm.pred = gsub(pattern = " ", replacement = "\n", ptPopGen$svm.pred)

ggplot(ptPopGen) + geom_boxplot(aes(x=svm.pred, y=piNpiS)) + ylim(c(0,1)) + 
  xlab("Predicted Compartment") + ylab("piN / piS")

ggplot(ptPopGen) + geom_boxplot(aes(x=svm.pred, y=upstream_pi)) + ylim(c(0,1)) +
  xlab("Predicted Compartment") + ylab("Upstream pi") + ylim(c(0,0.1))
ggplot(ptPopGen) + geom_violin(aes(x=svm.pred, y=upstream_indel)) + ylim(c(0,0.05)) +
  xlab("Predicted Compartment") + ylab("Upstream Indel")  # 11 and 15
ggplot(ptPopGen) + geom_count(aes(x=svm.pred, y=state)) + scale_size(range=c(0,10))

TukeyHSD(aov(piNpiS ~ svm.pred, data = ptPopGen)) # Surface Antigen
data.frame(TukeyHSD(aov(piNpiS ~ svm.pred, data = ptPopGen))$'svm.pred')[grep("unknown", rownames(TukeyHSD(aov(piNpiS ~ svm.pred, data = ptPopGen))$'svm.pred')),]
    # Surface antigen

# slightly different way
ptPopGen$state = ifelse(test = ptPopGen$state == "intact", yes = "intact", "LOF")
ggplot(ptPopGen) + geom_count(aes(x=svm.pred, y=state))

lof = ptPopGen[!is.na(ptPopGen$state),] %>% group_by(svm.pred) %>% summarise(LOF_Intact = sum(sum(state == "LOF")/sum(state == "intact")))
ggplot(data.frame(lof)) + geom_bar(aes(x=svm.pred, y=LOF_Intact), stat="identity") + theme_classic() + 
  xlab("Predicted Compartment") + ylab("")


lof_all = colSums(table(ptPopGen$svm.pred, ptPopGen$state))
chiScores = list()
for(comp in as.character(sort(unique(ptPopGen$svm.pred)))){
  tmpLOF = ptPopGen[which(ptPopGen$svm.pred == comp),]
  tmpDF = rbind(lof_all, table(tmpLOF$state))
  chiScores[[comp]] = chisq.test(tmpDF)$'p.value'
}

which(chiScores < 0.05/17) # Ax, BB-ass, MOM, Perox, Proteasome, SA




#-------------------------------------------------------------------------------------------------------------------------
# ##### Orthos
library(tidyverse)
library(reshape2)
library(RColorBrewer)
orthoTimeline = read.csv("/ParameciumPCP/tables/evolution/ptet-presence-euks.csv", header=T)
orthoTimeline = merge(orthoTimeline, fData(preds2)[,c('Accession', 'svm')], by='Accession')
orthSVMtimeline = data.frame(orthoTimeline %>% group_by(svm) %>% summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))
orthSVMtimeline$Accession = NULL
orthSVMtimeline$svm = gsub(pattern = " ", replacement = "\n", orthSVMtimeline$svm)
orthSVMspp = orthSVMtimeline[,colnames(orthSVMtimeline)[c(1,7:ncol(orthSVMtimeline))]]
orthSVMtimeline = orthSVMtimeline[,colnames(orthSVMtimeline)[1:6]]
mOrthoSVMtimeline = melt(orthSVMtimeline)

ggplot(mOrthoSVMtimeline, aes(x=svm, y=value, fill=variable)) + 
  geom_histogram(stat = "identity", position = "dodge") + theme_classic() + 
  xlab("Predicted Organellar Compartment") 

allOrtho = data.frame(orthoTimeline %>% group_by(svm) %>% summarise_all(mean))
allOrtho$Accession = NULL
allOrtho[,2:ncol(allOrtho)] = allOrtho[,2:ncol(allOrtho)]*100

write.csv(rbind(allOrtho, data.frame(orthoTimeline[,2:ncol(orthoTimeline)] %>% summarise_all(mean))*100), 
          file = "/ParameciumPCP/tables/evolution/ptetPCP-orthoPresence.csv", quote = F, row.names = F)

mOrthoSVMtimeline2 = matrix(ncol = 3, nrow = 1)
colnames(mOrthoSVMtimeline2) = c("svm", "variable", "value")

compss = unique(svmProps$svm)
for(compp in compss){
  tmpSVMTimeline = mOrthoSVMtimeline[which(mOrthoSVMtimeline$svm == compp),]
  tmpSVMTimeline$value = tmpSVMTimeline$value - min(tmpSVMTimeline$value)
  mOrthoSVMtimeline2 = rbind(mOrthoSVMtimeline2, tmpSVMTimeline)
}
mOrthoSVMtimeline2 = mOrthoSVMtimeline2[-1,]
mOrthoSVMtimeline2$variable <- factor(mOrthoSVMtimeline2$variable, 
                                      levels =c("Paramecium", "Ciliates", "Avleolates", "TSR", "BroadEuks"))

ggplot(mOrthoSVMtimeline, aes(x=svm, y=value, fill=variable)) + 
  geom_histogram(stat = "identity", position = "dodge") + theme_classic() + 
  xlab("Predicted Organellar Compartment") 

ggplot(mOrthoSVMtimeline2, aes(x=svm, y=value, fill=variable)) + 
  geom_histogram(stat = "identity", position = "dodge") + theme_classic() + 
  xlab("Predicted Organellar Compartment") 


# Stats
allOrthoCiliate = colSums(table(orthoTimeline$svm, orthoTimeline$Ciliates))
chiScores = list()
for(comp in as.character(sort(unique(orthoTimeline$svm)))){
  tmpOrthoTimeLine = orthoTimeline[which(orthoTimeline$svm == comp),]
  tmpDF = rbind(allOrthoCiliate, table(tmpOrthoTimeLine$Ciliates))
  chiScores[[comp]] = chisq.test(tmpDF)$'p.value'
}
which(chiScores < 0.05/17)  # BB core, Cyto, Lyso, MT-sol, Mito, Ribo, SA

allOrthoAlveolate = colSums(table(orthoTimeline$svm, orthoTimeline$Avleolates))
chiScores = list()
for(comp in as.character(sort(unique(orthoTimeline$svm)))){
  tmpOrthoTimeLine = orthoTimeline[which(orthoTimeline$svm == comp),]
  tmpDF = rbind(allOrthoAlveolate, table(tmpOrthoTimeLine$Avleolates))
  chiScores[[comp]] = chisq.test(tmpDF)$'p.value'
}
which(chiScores < 0.05/17)  # BB core, Cyto, Lyso, MT-sol, Mito, Prot, Ribo, SA

allOrthoTSR = colSums(table(orthoTimeline$svm, orthoTimeline$TSR))
chiScores = list()
for(comp in as.character(sort(unique(orthoTimeline$svm)))){
  tmpOrthoTimeLine = orthoTimeline[which(orthoTimeline$svm == comp),]
  tmpDF = rbind(allOrthoTSR, table(tmpOrthoTimeLine$TSR))
  chiScores[[comp]] = chisq.test(tmpDF)$'p.value'
}
which(chiScores < 0.05/17)  # BB core, Cyto, Lyso, MT-sol, Mito, Prot, Ribo, SA

allOrthoBroadEuks = colSums(table(orthoTimeline$svm, orthoTimeline$BroadEuks))
chiScores = list()
for(comp in as.character(sort(unique(orthoTimeline$svm)))){
  tmpOrthoTimeLine = orthoTimeline[which(orthoTimeline$svm == comp),]
  tmpDF = rbind(allOrthoBroadEuks, table(tmpOrthoTimeLine$BroadEuks))
  chiScores[[comp]] = chisq.test(tmpDF)$'p.value'
}
which(chiScores < 0.05/17)  # BB core, Cyto, Lyso, MT-sol, Mito, Prot, Ribo, SA, TM

# #groupings
eukProtCounts = read.csv("/ParameciumPCP/tables/evolution/ptet-presence-euks-oneByone.csv", header=T)
eukProtKey = melt(read.csv("/ParameciumPCP/tables/evolution/ptetPCP-presence-euks-supergroups.csv", header=T))
eukProtCounts = melt(cbind(eukProtCounts[,1], eukProtCounts[,colnames(eukProtCounts)[which(colnames(eukProtCounts) %in% eukProtKey$Name_to_Use)]]))
colnames(eukProtCounts)[1] = "Accession"

eukProt = merge(eukProtKey, eukProtCounts, by.x="Name_to_Use", by.y = "variable")
eukProt = merge(eukProt, svmProps[,c("svm", "Accession")], by="Accession")
eukProt$Supergroup_UniEuk = gsub(" ", "\n", eukProt$Supergroup_UniEuk)
eukProt$svm = gsub(" ", "\n", eukProt$svm)

tmpEukProt = eukProt[which(eukProt$Supergroup_UniEuk %in% unique(eukProt$Supergroup_UniEuk)[c(1,5, 9, 11)]),]

tmpEukProt = data.frame(tmpEukProt %>%
                          group_by(Accession, Supergroup_UniEuk) %>%
                          summarise(hasOrtho = any(value != 0), SVM = svm))
tmpEukProt = tmpEukProt[-which(duplicated(tmpEukProt)),]

#tmpEukProt$hasOrtho = ifelse(tmpEukProt$hasOrtho == TRUE, 1, 0)

ggplot(tmpEukProt, aes(x=SVM, y=hasOrtho, fill=Supergroup_UniEuk)) +
  geom_bar(stat = "identity", position = 'dodge')

data.frame(tmpEukProt %>% group_by(SVM, Supergroup_UniEuk) %>% summarise(sum(hasOrtho==TRUE)/(sum(hasOrtho==TRUE) + sum(hasOrtho==FALSE))))

#-------------------------------------------------------------------------------------------------------------------------
##### Toxo
library(pRolocdata)
data("Barylyuk2020ToxoLopit")
ptetPreds = fData(preds2)[,c("Accession", "svm.pred")]
toxoPreds = fData(Barylyuk2020ToxoLopit)[,c("Accession", "tagm.map.allocation.pred")]
toxoeuk = read.csv("/ParameciumPCP/tables/evolution/toxo-eukprot.csv", header=T)
ptetToxoEuk = read.csv("/ParameciumPCP/tables/evolution/ptet-toxo-eukprot.csv", header=T)

toxoeuk = toxoeuk[-which(toxoeuk$EukProt == "-"),]
ptetToxoEuk = ptetToxoEuk[-which(ptetToxoEuk$Toxoplasma.gondii == "-"),]
toxoPtet = merge(toxoeuk, ptetToxoEuk, by.x="EukProt", by.y="Toxoplasma.gondii")
toxoPtet = toxoPtet[which(!duplicated(toxoPtet)),]
toxoPtet$GeneID = gsub(pattern = "-t26_1-p1", replacement = "", toxoPtet$GeneID)
toxoPtet$EukProt = NULL
colnames(toxoPtet) = c("toxo", "ptet")
toxoPtet = merge(merge(toxoPtet, toxoPreds, by.x="toxo", by.y="Accession"), ptetPreds, by.x = "ptet", by.y="Accession")
# write.csv(data.frame(table(toxoPtet$tagm.map.allocation.pred, toxoPtet$svm.pred)), file = "/ParameciumPCP/tables/evolution/toxoPtet.csv", 
#           quote = F, row.names = F)
# write.csv(toxoPtet, file = "/ParameciumPCP/tables/evolution/toxoPtet2.csv", 
#           quote = F, row.names = F)

toxoPtet$svm.pred = gsub(" ", "\n", toxoPtet$svm.pred)

ggplot(toxoPtet, aes(svm.pred, tagm.map.allocation.pred)) + geom_count() + scale_size(range=c(0,10)) + 
  theme_bw() + theme(legend.position = "none") + 
  xlab("Predicted P. tetraurelia Compartment") + ylab("Predicted T. gondii Compartment")

### made a version of files collapsing all organelles split for reasons of clustering (e.g, nucleus- chromatin/non-chromatin == nucleus)
toxoPtet3 = read.csv("/ParameciumPCP/tables/evolution/toxoPtet3.csv", header=T)

ggplot(toxoPtet3, aes(svm.pred, tagm.map.allocation.pred)) + geom_count() + scale_size(range=c(0,10)) + 
  theme_bw() + theme(legend.position = "none") + 
  xlab("Predicted P. tetraurelia Compartment") + ylab("Predicted T. gondii Compartment")


# get seqs
library(seqinr)
toxoFa = read.fasta("/ParameciumPCP/sequences/ToxoDB-57_TgondiiME49_AnnotatedProteins.fasta", as.string = T, seqtype = "AA")
write.fasta(file.out = "/ParameciumPCP/sequences/toxoLOPITprots.fasta",
  sequences = getSequence(toxoFa[which(as.character(lapply(sapply(names(toxoFa), strsplit, "-"), "[", 1)) %in% 
        rownames(fData(Barylyuk2020ToxoLopit)))], as.string = T), 
  names = getName(toxoFa[which(as.character(lapply(sapply(names(toxoFa), strsplit, "-"), "[", 1)) %in% 
                       rownames(fData(Barylyuk2020ToxoLopit)))]))
  # for BUSCO... remember to load hmmer/3.3.1 first

#-------------------------------------------------------------------------------------------------------------------------
# Dinos
dino = read.csv("/ParameciumPCP/tables/evolution/ptet-presence-euks-oneByone.csv", header=T)[,c("Paramecium_tetrareulia_ID", "Oxyrrhis_marina")]
merge(svmProps[,c("Accession", "svm.pred")], dino, by.x = "Accession", by.y = "Paramecium_tetrareulia_ID") %>% 
  group_by(svm.pred) %>% summarise(tot = sum(Oxyrrhis_marina))


#############################################################################
# Project predictions onto non-imputed data

preds2_nonImpute_tab = cbind(fData(preds2)[which(rownames(fData(preds2)) %in% rownames(fData(ptetPCP_raw_all_norm_NA))),], 
                             exprs(preds2)[which(rownames(exprs(preds2)) %in% rownames(exprs(ptetPCP_raw_all_norm_NA))),])

write.csv(preds2_nonImpute_tab, file = "/ParameciumPCP/tables/ptetPCP-preds2-noImpute.csv", quote = F)

preds2_nonImpute = filterNA(readMSnSet2("/ParameciumPCP/tables/ptetPCP-preds2-noImpute.csv",
                     ecol = c("MAC.1", "X300g.1", "X1K.1", "X3K.1", "X5K.1", "X9K.1", 
                              "X12K.1", "X15K.1", "X30K.1", "X79K.1", "X120K.1", "SUP.1", 
                              "MAC.2", "X300g.2", "X1K.2", "X3K.2", "X5K.2", "X9K.2", 
                              "X12K.2", "X15K.2", "X30K.2", "X79K.2", "X120K.2", "SUP.2",
                              "MAC.3", "X300g.3", "X1K.3", "X3K.3", "X5K.3", "X9K.3", 
                              "X12K.3", "X15K.3", "X30K.3", "X79K.3", "X120K.3", "SUP.3"),                                           
                     fnames = "Accession"))
set.seed(42)
plot2D(preds2_nonImpute, fcol = "svm.pred", method = "t-SNE")
plot2D(preds2_nonImpute, fcol = "svm.pred")


#------------------------------------------------------------------------------------------------------------------------------------
# Ribosomes

ribos = as.character(unlist(sapply(X = read_lines("/ParameciumPCP/tables/Ribosomes/human-ptet-ribo-blast.txt"), FUN = strsplit, " ")))
ribos = as.character(unlist(sapply(ribos, strsplit, "\\|")))
ribos = ribos[grep(pattern = "*_HUMAN|GSPATT*", ribos)]


finalDF = data.frame(matrix(ncol=2, nrow=1))
colnames(finalDF) = c("HUMAN", "PTET")

for(riboN in 1:length(ribos)){
  riboM = riboN+1
  riboNM = ribos[c(riboN, riboM)]
  
  riboNname = ribos[riboN]
  riboMname = ribos[riboM]
  
  tmpDF = data.frame(matrix(ncol=2, nrow=1))
  colnames(tmpDF) = c("HUMAN", "PTET")
  
  if(length(grep("*_HUMAN", riboNM)) == 2){
    tmpDF$HUMAN = riboNname
    tmpDF$PTET = "absent"
    riboMnameSave = riboMname
  } 
  
  if(length(grep("*_HUMAN", riboNM)) == 1){
    if(grep("GSPATT*", riboNM) == 1){
      tmpDF$HUMAN = riboMnameSave
      tmpDF$PTET = riboNname
    }
    if(grep("GSPATT*", riboNM) == 2){
      tmpDF$HUMAN = riboNname
      tmpDF$PTET = riboMname
      
      riboMnameSave = riboNname
    }
  } 
  
  if(length(grep("*_HUMAN", riboNM)) == 0){
    tmpDF$HUMAN = riboMnameSave
    tmpDF$PTET = riboMname
  } 
  
  finalDF = rbind(finalDF, tmpDF)
}
write.csv(finalDF, file = "/ParameciumPCP/tables/Ribosomes/ribo-present.csv", quote = F, row.names = F)

riboConvert = read.csv("/ParameciumPCP/tables/Ribosomes/ptet-riboConvert.csv", header=T)  # change gene names from genic to protein
riboPtet = merge(finalDF, riboConvert, by.x="PTET", by.y="ID")
riboPtet = merge(riboPtet, fData(preds2), by.x="Related.ID", by.y=0)
riboPtet = riboPtet[!duplicated(riboPtet),]

missingRibos = finalDF[which(finalDF$PTET == "absent"),]
colnames(missingRibos) = c("HUMAN", "Accession")
missingRibos$svm = NA
missingRibos$svm.pred = NA

write.csv(rbind(riboPtet[,c("HUMAN", "Accession", "svm", "svm.pred")], 
      missingRibos), file = "/ParameciumPCP/tables/Ribosomes/ptetPCP-human-ptetp-RibosIDed.csv", quote = F, row.names = F)

#--------------------------------------------------
# Ohnolog families
ptetProps = read.csv("/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", header= T)
ptetProps$ohnologFamily = x

ptetWGD = read.table("/Paramecium_Paraorthologs/ptetraurelia_mac_51_annotation_v2.0.WGD.tree", header=T) 
ptetWGD$NB = NULL
ptetWGD$comb <- apply( ptetWGD[ , colnames(ptetWGD) ] , 1 , paste , collapse = "-" )
ptetWGDsearch = ptetWGD[,"comb"]


for(ptetGene in ptetProps$Accession){
  ptetProps$ohnologFamily[which(ptetProps$Accession == ptetGene)] = grep(ptetGene, ptetWGDsearch)
}

# Remove duplicates
sort(table(ptetProps$ohnologFamily), decreasing = T)

nrow(ptetProps)  # 8995
length(unique(ptetProps$ohnologFamily))  # 6186
ptetProps_reduced = data.frame(ptetProps %>% group_by(ohnologFamily) %>% filter(nPSM == max(nPSM)))
nrow(ptetProps_reduced)  # 6218 ... a few ties
length(unique(ptetProps_reduced$ohnologFamily)) # v


ptetPCP_noIDdups_tab = merge(ptetProps_reduced, exprs(preds2), by.x="Accession", by.y=0)
table(ptetPCP_noIDdups_tab$svm.pred)
write.csv(ptetPCP_noIDdups_tab, file = "/ParameciumPCP/tables/ptetPCP-preds2-noIDdups.csv", quote = F)

preds2_noIDdups = filterNA(readMSnSet2("/ParameciumPCP/tables/ptetPCP-preds2-noIDdups.csv",
                                       ecol = c("MAC.1", "X300g.1", "X1K.1", "X3K.1", "X5K.1", "X9K.1", 
                                                "X12K.1", "X15K.1", "X30K.1", "X79K.1", "X120K.1", "SUP.1", 
                                                "MAC.2", "X300g.2", "X1K.2", "X3K.2", "X5K.2", "X9K.2", 
                                                "X12K.2", "X15K.2", "X30K.2", "X79K.2", "X120K.2", "SUP.2",
                                                "MAC.3", "X300g.3", "X1K.3", "X3K.3", "X5K.3", "X9K.3", 
                                                "X12K.3", "X15K.3", "X30K.3", "X79K.3", "X120K.3", "SUP.3"),                                           
                                       fnames = "Accession"))
set.seed(42)
plot2D(preds2_noIDdups, fcol = "svm.pred", method = "t-SNE")
plot2D(preds2_noIDdups, fcol = "markers", method = "t-SNE")

table(fData(preds2_noIDdups)$'markers')
table(fData(preds2_noIDdups)$'svm.pred')

#get sequences
ptProts = read.fasta("/Paramecium_Proteomes-Predicted/ptetraurelia_mac_51_annotation_v2.0.protein.fa", seqtype = "AA", as.string = T)

write.fasta(sequences = getSequence(ptProts[fData(preds2_noIDdups)$'Accession']), 
            names = getName(ptProts[fData(preds2_noIDdups)$'Accession']), 
            file.out = "/ParameciumPCP/sequences/ptetPCP-noDups.fasta")


# piN/piS of dominant vs nondominant copies?
ptPopGenDom = merge(ptetProps[,c("Accession", "ohnologFamily")], ptPopGen, by.x = "Accession", by.y="GeneA")

table(ptPopGenDom$state)

# state is intact
ptPopGenDomSummary = ptPopGenDom %>% group_by(ohnologFamily) %>%
  summarize(minPSM = min(nPSM), minPSMstate = state[which.min(nPSM)], 
            maxPSM = max(nPSM), maxPSMstate = state[which.max(nPSM)]) %>% 
  mutate(maxDom = ifelse(maxPSMstate == "intact" & minPSMstate != "intact", TRUE, FALSE))

ptPopGenDomYes = data.frame(ptPopGenDomSummary[which(ptPopGenDomSummary$maxDom == TRUE),])  # only a few interesting ones... 
ptPopGenDomYes
# ohnologFamily minPSM    minPSMstate maxPSM maxPSMstate maxDom
# 1              7    112            PTC    203      intact   TRUE
# 2            860    512            PTC    809      intact   TRUE
# 3           1485    196            PTC    309      intact   TRUE
# 4           2073    761 PTC;frameshift    868      intact   TRUE
# 5           2541    552 PTC;frameshift    752      intact   TRUE
# 6           2654     40 PTC;frameshift     47      intact   TRUE
# 7           2692     45 PTC;frameshift     52      intact   TRUE
# 8           3146    245 PTC;frameshift    292      intact   TRUE
# 9           3405     27            PTC     33      intact   TRUE
# 10          3550     10 PTC;frameshift    350      intact   TRUE
# 11          3831    337            PTC    462      intact   TRUE
# 12          3832     45 PTC;frameshift    191      intact   TRUE
# 13          7138    397 PTC;frameshift    445      intact   TRUE
# 14          7323     63 PTC;frameshift     66      intact   TRUE
# 15          7472    185     frameshift    280      intact   TRUE
# 16          8885    131 PTC;frameshift    161      intact   TRUE
# 17          8987   1694     frameshift   2470      intact   TRUE
# 18          9135     32 PTC;frameshift     59      intact   TRUE
# 19          9981    641 PTC;frameshift    783      intact   TRUE
# 20         10802    103 PTC;frameshift    520      intact   TRUE
# 21         10917     99            PTC    219      intact   TRUE
# 22         11164     25 PTC;frameshift     43      intact   TRUE
# 23         11265     75 PTC;frameshift     88      intact   TRUE
# 24         11386     34 PTC;frameshift     90      intact   TRUE
# 25         11463     75            PTC    110      intact   TRUE
# 26         11864    623 PTC;frameshift    626      intact   TRUE
# 27         17136     10 PTC;frameshift     13      intact   TRUE
# 28         17224    226 PTC;frameshift    502      intact   TRUE
# 29         17327     21 PTC;frameshift    101      intact   TRUE

ptPopGenDomYes$diff = ptPopGenDomYes$maxPSM - ptPopGenDomYes$minPSM
ptPopGenDomYes$diffP = ptPopGenDomYes$diff/sum(ptPopGenDomYes$minPSM+ptPopGenDomYes$maxPSM)
ptPopGenDomYes$diffP = ptPopGenDomYes$diffP*100

ptPopGenDomYes %>% arrange(diffP)


#-----------------------------------------------------------------------------------------
# Trichocyst Matrix proteins
trichy = ptetProps[grep("^TMP*", ptetProps$Aliases),]

trichyTMP = trichy %>% filter(svm == "Trichocyst Matrix") %>% pull(Accession)
trichyOTHER = trichy %>% filter(svm != "Trichocyst Matrix") %>% pull(Accession)

plotDist(preds2[trichyTMP], pcol = "black", ylim = c(0,1))
plotDist(preds2[trichyOTHER], pcol = "black", ylim = c(0,1))

#-----------------------------------------------------------------------------------------
# MAC/MIC

# RPB1
plotDist(preds2["PTET.51.1.P1370127"], pcol = "black", ylim = c(0,1))

# RPB2
plotDist(preds2["PTET.51.1.P0480005"], pcol = "black", ylim = c(0,1))

polFOI = read.table("/ParameciumPCP/tables/FOI/ptet_RNApol.tsv", sep = "\t")

plotDist(preds2[polFOI[which(polFOI$Protein.name %in% fData(preds2)$'Accession'),'Protein.name']], 
         pcol = "black", ylim = c(0,1))


#-----------------------------------------------------------------------------------------
# Paralogs

# casein kinase
CSNK2A = read.table("/ParameciumPCP/tables/FOI/ptet_CSNK2A.txt")
plotDist(preds2[CSNK2A$V1], pcol = "black", ylim = c(0,1))

CSNK2B = read.table("/ParameciumPCP/tables/FOI/ptet_CSNK2B.txt")
plotDist(preds2[CSNK2B$V1], pcol = "black", ylim = c(0,1))

CSNK2 = addMarkers(preds2, markers = "/ParameciumPCP/tables/FOI/ptet_CSNK2.csv", 
                   mcol = "Paralog", fcol = "Accession", verbose = T)

set.seed(42)
plot2D(CSNK2, fcol = "Paralog", method = "t-SNE", col = c("red", "blue"))

plotDist(preds2[c(CSNK2A$V1, CSNK2B$V1)], pcol = "black", ylim = c(0,1))


# isocitrate dehydrogenase
IDH = addMarkers(preds2, markers = "/ParameciumPCP/tables/FOI/ptet_IDH.csv", 
                   mcol = "Paralog2", fcol = "Accession", verbose = T)
set.seed(42)
plot2D(IDH, fcol = "Paralog2", method = "t-SNE", col = c("green4", "purple"))

#-----------------------------------------------------------------------------------------


