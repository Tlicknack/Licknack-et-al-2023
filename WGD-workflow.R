# PCP-WGD.R

library(pRoloc)
library(tidyverse)
library(reshape2)

### Important objects:
svmProps = read.csv("/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", header=T)
load("/ParameciumPCP/RData/WGD-All.Rdata")
###



# Workflow
source("/ParameciumPCP/scripts/source_wgdPairs.R")
load("/ParameciumPCP/RData/ptetPCP_superANDunsuper.RData")

ptetPCP = read.csv(file = "/ParameciumPCP/tables/ptetPCP.csv", header=T) 
ptetPCP = ptetPCP[grep(pattern = "PTET*", x = ptetPCP$Accession),]   
dim(ptetPCP)  # 11856 proteins X 248 columns (34 proteins removed)

ptetWGD = read.table("/Paramecium_Paraorthologs/ptetraurelia_mac_51_annotation_v2.0.WGD.tree", header=T) 
ptetWGD$NB = NULL

for(q in 1:ncol(ptetWGD)){  # make dots into NA
  ptetWGD[,q] = sub("^[.]$", NA, x = ptetWGD[,q])
}

# Single copy genes vs multicopy
ptet_singleCopy = ptetWGD[which(rowSums(ifelse(apply(X = ptetWGD, MARGIN = 2, FUN = is.na) == T, yes = 1, no = 0)) == 7),]
ptet_singleCopy = as.character(unlist(ptet_singleCopy)[which(complete.cases(unlist(ptet_singleCopy)))])
length(ptet_singleCopy)  # 12457 single copies
length(which(ptet_singleCopy %in% ptetPCP$Accession))  # 3041 single copy genes in data
  # 3041/12457 = ~24%
length(which(ptet_singleCopy %in% svmProps$Accession)) # 2102 single copy with classification
  # 2102/8995

ptetWGD = ptetWGD[-which(rowSums(ifelse(apply(X = ptetWGD, MARGIN = 2, FUN = is.na) == T, yes = 1, no = 0)) == 7),]
length(unique(as.character(unlist(ptetWGD))[which(complete.cases(unlist(ptetWGD)))]))  # 28003 multi copy
length(which(unique(as.character(unlist(ptetWGD))[which(complete.cases(unlist(ptetWGD)))]) %in% ptetPCP$Accession))  # 8815
  # 8815/28003 = ~31.5%

# Possible comparisons
# tmp
ptetWGD8 = ptetWGD[which(rowSums(ifelse(apply(X = ptetWGD, MARGIN = 2, FUN = is.na) == T, yes = 1, no = 0)) == 0),]
ptetWGD7 = ptetWGD[which(rowSums(ifelse(apply(X = ptetWGD, MARGIN = 2, FUN = is.na) == T, yes = 1, no = 0)) == 1),]
ptetWGD6 = ptetWGD[which(rowSums(ifelse(apply(X = ptetWGD, MARGIN = 2, FUN = is.na) == T, yes = 1, no = 0)) == 2),]
ptetWGD5 = ptetWGD[which(rowSums(ifelse(apply(X = ptetWGD, MARGIN = 2, FUN = is.na) == T, yes = 1, no = 0)) == 3),]
ptetWGD4 = ptetWGD[which(rowSums(ifelse(apply(X = ptetWGD, MARGIN = 2, FUN = is.na) == T, yes = 1, no = 0)) == 4),]
ptetWGD3 = ptetWGD[which(rowSums(ifelse(apply(X = ptetWGD, MARGIN = 2, FUN = is.na) == T, yes = 1, no = 0)) == 5),]
ptetWGD2 = ptetWGD[which(rowSums(ifelse(apply(X = ptetWGD, MARGIN = 2, FUN = is.na) == T, yes = 1, no = 0)) == 6),]

length(unlist(ptetWGD2)[which(complete.cases(unlist(ptetWGD2)))])*1
length(unlist(ptetWGD3)[which(complete.cases(unlist(ptetWGD3)))])*3
length(unlist(ptetWGD4)[which(complete.cases(unlist(ptetWGD4)))])*6
length(unlist(ptetWGD5)[which(complete.cases(unlist(ptetWGD5)))])*10
length(unlist(ptetWGD6)[which(complete.cases(unlist(ptetWGD6)))])*15
length(unlist(ptetWGD7)[which(complete.cases(unlist(ptetWGD7)))])*21
length(unlist(ptetWGD8)[which(complete.cases(unlist(ptetWGD8)))])*28


# Determine presence/absence- calculate euclidean distance if both present
ptetPCP_Chi = ptetPCP[,grep("Found.in.Sample*", colnames(ptetPCP))]
ptetPCP_Chi = cbind(ptetPCP[,2], ptetPCP_Chi)

for(q in 1:ncol(ptetPCP_Chi)){  # code values
  ptetPCP_Chi[,q] = sub("Not Found", 0, x = ptetPCP_Chi[,q])
  ptetPCP_Chi[,q] = sub("Peak Found", 1, x = ptetPCP_Chi[,q])
  ptetPCP_Chi[,q] = sub("High", 2, x = ptetPCP_Chi[,q])
}

chiScores = data.frame(matrix(ncol=5, nrow=1))
colnames(chiScores) = c("Type", "Status", "GeneA", "GeneB", "EuclideanDistance")

euclid = function(a,b) sqrt(sum((a-b)^2))


for(i in 1:nrow(ptetWGD)){
  ptetWGDi = ptetWGD[i,which(!is.na(ptetWGD[i,]))]
  
  for(wgdCombo in combn(colnames(ptetWGDi), 2, simplify = F)){
    ptetWGDic = ptetWGDi[,unlist(wgdCombo)]
    
    tmpChiDF = data.frame(matrix(ncol=5, nrow=1))
    colnames(tmpChiDF) = c("Type", "Status", "GeneA", "GeneB", "EuclideanDistance")
    tmpChiDF$GeneA = ptetWGDic[,1]
    tmpChiDF$GeneB = ptetWGDic[,2]
    
    whichWGD = as.numeric(which(duplicated(rbind(t(data.frame(all_wgd)), 
                                                 t(data.frame(unlist(list(colnames(ptetWGDic)))))), fromLast=T)))
    if(whichWGD <= 4){
      tmpChiDF$Type = "WGD1"
    }
    if(whichWGD > 4 & whichWGD <= 12){
      tmpChiDF$Type = "WGD2"
    }
    if(whichWGD > 12){
      tmpChiDF$Type = "WGD3"
    }
    
    wgdInData = which(ptetPCP_Chi$`ptetPCP[, 2]` %in% as.character(ptetWGDic))
    
    if(length(wgdInData) == 0){ # both copies are missing from the data
      tmpChiDF$Status = "Both Absent"
      tmpChiDF$EuclideanDistance = NA
    }
    if(length(wgdInData) == 1){ # one copy is present in the data
      tmpChiDF$Status = "One Absent"
      tmpChiDF$EuclideanDistance = NA
    }
    if(length(wgdInData) == 2){ # both copies are present in the data
      tmpChiDF$Status = "Both Present"
      tmpChiDF$EuclideanDistance = euclid(as.numeric(ptetPCP_Chi[wgdInData[1], 2:ncol(ptetPCP_Chi)]), 
                                          as.numeric(ptetPCP_Chi[wgdInData[2], 2:ncol(ptetPCP_Chi)]))
    }
    chiScores = rbind(chiScores, tmpChiDF)
  }
}
chiScores = chiScores[-1,]
nrow(chiScores)
summary(chiScores$EuclideanDistance)

write.table(unique(c(chiScores$GeneA, chiScores$GeneB)), "/ParameciumPCP/tables/WGD-gene-list.txt", quote = F, row.names = F)

# Simulate tens of thousands of 
set.seed(42)
ran1 = sample(1:nrow(ptetPCP_Chi), 100000, replace=TRUE)
ran2 = sample(1:nrow(ptetPCP_Chi), 100000, replace=TRUE)

ranchiScores = data.frame(matrix(ncol=5, nrow=100000))
colnames(ranchiScores) = c("Type", "Status", "GeneA", "GeneB", "EuclideanDistance")
ranchiScores$Type = rep("Random", 100000)
ranchiScores$Status = rep("Random", 100000)

for(r in 1:100000){
  ranchiScores$GeneA[r] = ptetPCP_Chi[ran1[r], 1]
  ranchiScores$GeneB[r] = ptetPCP_Chi[ran2[r], 1]
  ranchiScores$EuclideanDistance[r] = euclid(as.numeric(ptetPCP_Chi[ran1[r], 2:ncol(ptetPCP_Chi)]), 
                                             as.numeric(ptetPCP_Chi[ran2[r], 2:ncol(ptetPCP_Chi)]))
}

# We want the lowest 5 percentile  
length(which(complete.cases(chiScores$EuclideanDistance)))  # 4214 complete pairs
hist(ranchiScores$EuclideanDistance, breaks=100)
hist(chiScores$EuclideanDistance, breaks=100)

wgdSimilar = chiScores[which(chiScores$EuclideanDistance < as.numeric(quantile(ranchiScores$EuclideanDistance, 0.05))),]
wgdDifferent = chiScores[-which(chiScores$EuclideanDistance < as.numeric(quantile(ranchiScores$EuclideanDistance, 0.05))),]

write.csv(wgdSimilar, file = "/ParameciumPCP/tables/ptetPCP_WGD-similar.csv", quote = F, row.names = F)
write.csv(wgdDifferent, file = "/ParameciumPCP/tables/ptetPCP_WGD-different.csv", quote = F, row.names = F)
ranchiScores$zscore = NA

# Summary
nrow(chiScores)
table(chiScores$Type)
table(chiScores$Status)
nrow(wgdDifferent[which(complete.cases(wgdDifferent)),])
nrow(wgdSimilar[which(complete.cases(wgdSimilar)),])
table(wgdDifferent[which(complete.cases(wgdDifferent)),'Type'])
table(wgdSimilar[which(complete.cases(wgdSimilar)),'Type'])

ggplot(rbind(ranchiScores, chiScores), aes(x=EuclideanDistance, fill = Status)) + 
  geom_histogram(data=subset(rbind(ranchiScores, chiScores),Status == 'Random'), fill = "red", alpha = 0.5, bins = 100) + 
  geom_histogram(data=subset(rbind(ranchiScores, chiScores),Status == 'Both Present'), fill = "blue", alpha = 0.5, bins = 100) + 
  geom_vline(xintercept = as.numeric(quantile(ranchiScores$EuclideanDistance, 0.05)), color = "red", linetype = "dotted", size = 1) + 
  theme_classic() + xlab("Euclidean Distance Between Protein Pairs") + ylab("Number of Protein Pairs") + 
  geom_vline(xintercept = 17.86057, color = "darkred", linetype = "dotted", size = 1)

# Calculate Z score
chiScores = chiScores %>% 
  mutate(zscore = ((EuclideanDistance - mean(ranchiScores$EuclideanDistance)) / sd(ranchiScores$EuclideanDistance))) %>%
  arrange(zscore)

chiScores[which(chiScores$zscore > 1.65),]  # 95% confidence --> 48 pairs
chiScores[which(chiScores$zscore > 1.96),]  # 97.5% confidence --> 31 pairs
chiScores[which(chiScores$zscore > 2.34),]  # 99% confidence --> 18 pairs
chiScores[which(chiScores$zscore > 2.80),]  # 99.75% confidence --> 4 pairs

# make fasta (align and measure Pdist in mega)
library(seqinr)
ptProts = read.fasta("/Paramecium_Proteomes-Predicted/ptetraurelia_mac_51_annotation_v2.0.protein.fa", seqtype = "AA", as.string = T)

for(f in 1:nrow(chiScores)){
  chiRow = chiScores[f,]
  write.fasta(sequences = getSequence(ptProts[c(chiRow$GeneA, chiRow$GeneB)], as.string = T), 
              names = getName(ptProts[c(chiRow$GeneA, chiRow$GeneB)]), 
              file.out = paste("/ParameciumPCP/sequences/WGD-all/", chiRow$GeneA, "_", chiRow$GeneB, ".fasta", sep = ""))
}

# Compare evolutionary divergence values in incomplete pairs vs complete
# retrieve from mega
megFiles = list.files("/ParameciumPCP/sequences/WGD-all/", full.names = T)
megFiles = megFiles[grep(pattern = "*[0-9].meg", x = megFiles)]
FmegTab = data.frame(matrix(ncol = 3, nrow = 1))
colnames(FmegTab) = c("GeneA", "GeneB", "Sequence_Distance")

for(meg in megFiles){
  megTab = read.table(meg, skip = 37)[2]  # .meg files with 2 seqs have 37 lines of unimportance
  colnames(megTab) = "Sequence_Distance"
  megTab$GeneA = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[1]
  megTab$GeneB = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[2]
  FmegTab = rbind(FmegTab, megTab)
}
FmegTab = FmegTab[2:nrow(FmegTab),]
FmegTab = merge(FmegTab, chiScores)


ggplot(FmegTab[FmegTab$Status != "Both Absent",], aes(x=Status, y=Sequence_Distance)) + geom_boxplot() + 
  theme_classic() + xlab("Status of Ohnologs") + ylab("Evolutionary Distance (Poisson Corrected)") + facet_wrap(~ Type)

summary(aov(Sequence_Distance ~ Status, data = FmegTab[FmegTab$Status != "Both Absent",]))
TukeyHSD(aov(Sequence_Distance ~ Status, data = FmegTab[FmegTab$Status != "Both Absent",]))
summary(aov(Sequence_Distance ~ Type, data = FmegTab[FmegTab$Status != "Both Absent",]))
TukeyHSD(aov(Sequence_Distance ~ Type, data = FmegTab[FmegTab$Status != "Both Absent",]))
summary(aov(Sequence_Distance ~ Status * Type, data = FmegTab[FmegTab$Status != "Both Absent",]))
TukeyHSD(aov(Sequence_Distance ~ Status * Type, data = FmegTab[FmegTab$Status != "Both Absent",]))

# add expression level
svmProps = read.csv("/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", header=T)
FmegRNAtab = merge(merge(FmegTab, svmProps[,c("Accession", "mRNA")], by.x = "GeneA", by.y = "Accession"), 
      svmProps[,c("Accession", "mRNA")], by.x = "GeneB", by.y = "Accession")

ggplot(FmegRNAtab) + geom_point(aes(x=mRNA.x, y=mRNA.y)) +   
  scale_y_log10() + scale_x_log10() + xlab("mRNA Expression (GeneA)") + ylab("mRNA Expression (GeneB)") + theme_classic()

summary(FmegRNAtab$mRNA.x)
summary(FmegRNAtab$mRNA.y)
FmegRNAtab[which(FmegRNAtab$mRNA.y < 1),]

length(FmegRNAtab[which(FmegRNAtab$mRNA.x < as.numeric(quantile(svmProps$mRNA, 0.05))),])
length(FmegRNAtab[which(FmegRNAtab$mRNA.y < as.numeric(quantile(svmProps$mRNA, 0.05))),])


#save.image("/ParameciumPCP/RData/WGD-All.Rdata")
load("/ParameciumPCP/RData/WGD-All.Rdata")

ggplot(FmegTab, aes(x=Sequence_Distance, y=EuclideanDistance)) + geom_point() + geom_smooth(method = "lm") +
  xlab("Sequence Divergence") + ylab("Euclidean Distance") + theme_classic()

ggplot(FmegTab, aes(x=Sequence_Distance, y=EuclideanDistance)) + geom_point() + geom_smooth(method = "lm") +
  scale_y_log10() + scale_x_log10() + xlab("Sequence Divergence") + ylab("Euclidean Distance") + theme_classic()

 
# ggplot(FmegTab, aes(x=Sequence_Distance, y=EuclideanDistance)) + geom_point() + geom_smooth(method = "lm") +
#   scale_y_log10() + scale_x_log10() + facet_wrap(~ Type)
# 
# ggplot(FmegTab[which(FmegTab$EuclideanDistance < as.numeric(quantile(ranchiScores$EuclideanDistance, 0.05))),], 
#        aes(x=Sequence_Distance, y=EuclideanDistance)) +
#   geom_point() + geom_smooth(method = "lm") +
#   scale_y_log10() + scale_x_log10()
# 
# ggplot(FmegTab[which(FmegTab$EuclideanDistance > as.numeric(quantile(ranchiScores$EuclideanDistance, 0.05))),], 
#        aes(x=Sequence_Distance, y=EuclideanDistance)) +
#   geom_point() + geom_smooth(method = "lm") +
#   scale_y_log10() + scale_x_log10()


# load localization data: redo analysis on differential WGD retention and organelle
# WGD Compartment enrichment
library(data.table)
library(RColorBrewer)
svmProps = read.csv("/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", header=T)
#svmProps$svm.pred = gsub(pattern = " ", replacement = "\n", svmProps$svm.pred)  # only for ggplot

tmp = svmProps[, c("svm.pred", "WGD1", "WGD2", "WGD3")]
tmp$WGD1 = ifelse(tmp$WGD1 == TRUE, 1, 0)
tmp$WGD2 = ifelse(tmp$WGD2 == TRUE, 1, 0)
tmp$WGD3 = ifelse(tmp$WGD3 == TRUE, 1, 0)

df <- melt(data.table(tmp), id.vars = "svm.pred", variable.name = "WGD", value.name="T_F")
ggplot(df, aes(x = svm.pred , y= T_F, fill = WGD)) +
  geom_bar(position="dodge", stat = "summary") + scale_fill_manual(values=brewer.pal(n=3, name="Set1")) + 
  ylim(c(0,1)) + xlab("Predicted Compartment") + ylab("WGD1/2/3 Presence") + theme(legend.position = "none") +
  geom_hline(yintercept=0.6988, color = "red", linetype = "dashed") + 
  geom_hline(yintercept=0.3103, color = "blue3", linetype = "dashed") + 
  geom_hline(yintercept=0.06281, color = "forestgreen", linetype = "dashed")

# stats
# WGD1
allProts = tmp %>% summarise(true = sum(WGD1 == 1, na.rm = T), false = sum(WGD1 == 0, na.rm = T))
rownames(allProts) = "all"
groupedProts = tmp %>% group_by(svm.pred) %>% 
  summarise(true = sum(WGD1 == 1, na.rm = T),
            false = sum(WGD1 == 0, na.rm = T))
comps = unique(fData(preds2)$'svm.pred')
chiScores = list()
for(compp in comps){
  tmp1 = data.frame(groupedProts[which(groupedProts$svm.pred == compp),c(2:3)])
  rownames(tmp1) = compp
  tmp2 = rbind(tmp1, allProts)
  chii = chisq.test(tmp2)
  chiScores[[compp]] = chii$p.value
} 
which(chiScores < 0.0025)  # p > 0.0025: lysosome, insol MT, Ribosome, ER, BB-core, Trich, SA, Proteasome


# WGD2
allProts = tmp %>% summarise(true = sum(WGD2 == 1, na.rm = T), false = sum(WGD2 == 0, na.rm = T))
rownames(allProts) = "all"
groupedProts = tmp %>% group_by(svm.pred) %>% 
  summarise(true = sum(WGD2 == 1, na.rm = T),
            false = sum(WGD2 == 0, na.rm = T))
comps = unique(fData(preds2)$'svm.pred')
chiScores = list()
for(compp in comps){
  tmp1 = data.frame(groupedProts[which(groupedProts$svm.pred == compp),c(2:3)])
  rownames(tmp1) = compp
  tmp2 = rbind(tmp1, allProts)
  chii = chisq.test(tmp2)
  chiScores[[compp]] = chii$p.value
}
which(chiScores < 0.0025)  # p > 0.0025: Ribosome, sol nuclei, cytosol, Trich matrix, Proteasome

# WGD3
allProts = tmp %>% summarise(true = sum(WGD3 == 1, na.rm = T), false = sum(WGD3 == 0, na.rm = T))
rownames(allProts) = "all"
groupedProts = tmp %>% group_by(svm.pred) %>% 
  summarise(true = sum(WGD3 == 1, na.rm = T),
            false = sum(WGD3 == 0, na.rm = T))
comps = unique(fData(preds2)$'svm.pred')
chiScores = list()
for(compp in comps){
  tmp1 = data.frame(groupedProts[which(groupedProts$svm.pred == compp),c(2:3)])
  rownames(tmp1) = compp
  tmp2 = rbind(tmp1, allProts)
  chii = chisq.test(tmp2)
  chiScores[[compp]] = chii$p.value
}
which(chiScores < 0.0025)  # p > 0.0025: mito, Ribosome, sol nuclei, TM, BB-ass


# Evol Translocation
svmProps = read.csv("/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", header=T)
#svmProps$svm.pred = gsub(pattern = " ", replacement = "\n", svmProps$svm.pred)  # only for ggplot
#svmProps$svm = gsub(pattern = " ", replacement = "\n", svmProps$svm)  # only for ggplot

wgdDF = read.csv("/ParameciumPCP/tables/ptetPCP-WGD-PPSS-Pdist.csv", header=T)
table(wgdDF$sameClass)
table(wgdDF[wgdDF$Type == "WGD1",]$'sameClass')

summary(aov(sameClass ~ Sequence_Distance, wgdDF))
ggplot(wgdDF) + geom_boxplot(aes(x=sameClass, y=Sequence_Distance))
ggplot(wgdDF) + geom_boxplot(aes(x=sameClass, y=Sequence_Distance)) + facet_wrap(~ Type) + 
  xlab("Ohnolog in Same Class") + ylab("Sequence Distance (Poisson Model)") + ggtitle("Whole Protein")

summary(aov(sameClass ~ Sequence_Distance, wgdDF[which(wgdDF$Type == "WGD1"),]))
summary(aov(sameClass ~ Sequence_Distance, wgdDF[which(wgdDF$Type == "WGD2"),]))
summary(aov(sameClass ~ Sequence_Distance, wgdDF[which(wgdDF$Type == "WGD3"),]))


mostDivgernt_all = merge(rbind(
  wgdDF[which(wgdDF$PPSS %in% sort(wgdDF[which(wgdDF$Type == "WGD1"),'PPSS'])),],
  wgdDF[which(wgdDF$PPSS %in% sort(wgdDF[which(wgdDF$Type == "WGD2"),'PPSS'])),],
  wgdDF[which(wgdDF$PPSS %in% sort(wgdDF[which(wgdDF$Type == "WGD3"),'PPSS'])),]), 
  svmProps[,c('Accession', 'svm', 'svm.pred')], by.x="GeneA", by.y = "Accession")
mostDivgernt_all = merge(mostDivgernt_all, svmProps[,c('Accession', 'svm', 'svm.pred')], by.x="GeneB", by.y="Accession")
colnames(mostDivgernt_all)[5:8] = c("GeneA_svm", "GeneA_svm.pred", "GeneB_svm", "GeneB_svm.pred")

mostDivgernt_all$GeneA_svm = gsub(pattern = " ", replacement = "\n", mostDivgernt_all$GeneA_svm)  # only for ggplot
mostDivgernt_all$GeneB_svm.pred = gsub(pattern = " ", replacement = "\n", mostDivgernt_all$GeneB_svm.pred)  # only for ggplot
mostDivgernt_all$GeneA_svm.pred = gsub(pattern = " ", replacement = "\n", mostDivgernt_all$GeneA_svm.pred)  # only for ggplot
mostDivgernt_all$GeneB_svm = gsub(pattern = " ", replacement = "\n", mostDivgernt_all$GeneB_svm)  # only for ggplot

ggplot(mostDivgernt_all) + geom_count(aes(x=GeneA_svm.pred, y=GeneB_svm.pred, color = ..n..)) + scale_size(range=c(0,10)) + 
  xlab("Predicted Compartment (Gene A)") + ylab("Predicted Compartment (Gene B)") + theme(legend.position = "none")
ggplot(mostDivgernt_all) + geom_count(aes(x=GeneA_svm, y=GeneB_svm, color = ..n..)) + scale_size(range=c(0,10)) + 
  xlab("Classified Compartment (Gene A)") + ylab("Classified Compartment (Gene B)") + theme(legend.position = "none")


### Gene Duplicate Paralog Analysis
svmProps = read.csv("/ParameciumPCP/tables/properties/ptetPCP_Properties.csv", header=T)
ghostKO = read.csv("/ParameciumPCP/tables/properties/ghostKOALA.csv", header=T)
paralogs = merge(svmProps, 
                 ghostKO, 
                 by.y = "Query", 
                 by.x ="Accession")[,c("Accession", "PutativeGeneName", "svm", "mRNA", "DE_Reciliation", "DE_TrichocystDischarge", "TMD", "SP", "TargetP")]
nrow(paralogs)
length(unique(paralogs$PutativeGeneName)) #1921 genes
length(which(paralogs$PutativeGeneName == "")) # 4825

arrange(paralogs %>% group_by(PutativeGeneName) %>% summarise(size = length(PutativeGeneName)), 
        desc(size), group_by = PutativeGeneName)


ggplot(paralogs) + geom_count(aes(x=svm, y=PutativeGeneName))

diffParalogs = paralogs %>% group_by(PutativeGeneName) %>% summarise(nDiff = length(unique(svm)))
table(diffParalogs$nDiff)

diffParalogs[which(diffParalogs$nDiff > 5),]



######## N-terminal, C-terminal distances
library(seqinr)
fastas = list.files("/ParameciumPCP/sequences/WGD", full.names = T)

for(fasta in fastas){
  tmpFasta = read.fasta(fasta, seqtype = "AA", as.string = T)
  tmpFasta1 = tmpFasta[1]
  tmpFasta2 = tmpFasta[2]
  
  tmpFasta1_N = substr(unlist(getSequence(tmpFasta1, as.string = T)), 
                       start = 1, stop = nchar(getSequence(tmpFasta1, as.string = T))/2)
  tmpFasta1_C = substr(unlist(getSequence(tmpFasta1, as.string = T)), 
                       start = nchar(getSequence(tmpFasta1, as.string = T))/2, stop = nchar(getSequence(tmpFasta1, as.string = T)))
  tmpFasta2_N = substr(unlist(getSequence(tmpFasta2, as.string = T)), 
                       start = 1, stop = nchar(getSequence(tmpFasta2, as.string = T))/2)
  tmpFasta2_C = substr(unlist(getSequence(tmpFasta2, as.string = T)), 
                       start = nchar(getSequence(tmpFasta2, as.string = T))/2, stop = nchar(getSequence(tmpFasta2, as.string = T)))
  
  nameN = paste("/ParameciumPCP/sequences/WGD-Nterm/N_", getName(tmpFasta1), "_", getName(tmpFasta2), ".fasta", sep = "")
  write.fasta(sequences = list(tmpFasta1_N, tmpFasta2_N), names = getName(tmpFasta), file.out = nameN)
  nameC = paste("/ParameciumPCP/sequences/WGD-Cterm/C_", getName(tmpFasta1), "_", getName(tmpFasta2), ".fasta", sep = "")
  write.fasta(sequences = list(tmpFasta1_C, tmpFasta2_C), names = getName(tmpFasta), file.out = nameC)

}

# Nterm
megFilesN = list.files("/ParameciumPCP/sequences/WGD-Nterm-dists/", full.names = T)
megFilesN = megFilesN[grep(pattern = "*[0-9].meg", x = megFilesN)]
FmegTabN = data.frame(matrix(ncol = 3, nrow = 1))
colnames(FmegTabN) = c("GeneA", "GeneB", "Sequence_Distance")

for(meg in megFilesN){
  megTab = read.table(meg, skip = 37)[2]  
  colnames(megTab) = "Sequence_Distance"
  meg = gsub("N_", "", meg)
  megTab$GeneA = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[1]
  megTab$GeneB = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[2]
  FmegTabN = rbind(FmegTabN, megTab)
}
FmegTabN = FmegTabN[2:nrow(FmegTabN),]
wgdDFn = wgdDF[,c(1:8,10:11)]
FmegTabN = FmegTabN %>% right_join(wgdDFn, by=c("GeneA","GeneB"))

summary(aov(sameClass ~ Sequence_Distance, FmegTabN))
ggplot(FmegTabN) + geom_boxplot(aes(x=sameClass, y=Sequence_Distance))
ggplot(FmegTabN) + geom_boxplot(aes(x=sameClass, y=Sequence_Distance)) + facet_wrap(~ Type) + 
  xlab("Ohnolog in Same Class") + ylab("Sequence Distance (Poisson Model)") + ggtitle("N-Terminus")

summary(aov(sameClass ~ Sequence_Distance, FmegTabN[which(FmegTabN$Type == "WGD1"),])) # 0.301
summary(aov(sameClass ~ Sequence_Distance, FmegTabN[which(FmegTabN$Type == "WGD2"),])) # 0.000198
summary(aov(sameClass ~ Sequence_Distance, FmegTabN[which(FmegTabN$Type == "WGD3"),])) # 0.000439


# Cterm
megFilesC = list.files("/ParameciumPCP/sequences/WGD-Cterm-dists/", full.names = T)
megFilesC = megFilesC[grep(pattern = "*[0-9].meg", x = megFilesC)]
FmegTabC = data.frame(matrix(ncol = 3, nrow = 1))
colnames(FmegTabC) = c("GeneA", "GeneB", "Sequence_Distance")

for(meg in megFilesC){
  megTab = read.table(meg, skip = 37)[2]  
  colnames(megTab) = "Sequence_Distance"
  meg = gsub("C_", "", meg)
  megTab$GeneA = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[1]
  megTab$GeneB = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[2]
  FmegTabC = rbind(FmegTabC, megTab)
}
FmegTabC = FmegTabC[2:nrow(FmegTabC),]
FmegTabC = merge(FmegTabC, chiScores)
wgdDFc = wgdDF[,c(1:8,10:11)]
FmegTabC = FmegTabC %>% right_join(wgdDFc, by=c("GeneA","GeneB"))

summary(aov(sameClass ~ Sequence_Distance, FmegTabC))
ggplot(FmegTabC) + geom_boxplot(aes(x=sameClass, y=Sequence_Distance))
ggplot(FmegTabC) + geom_boxplot(aes(x=sameClass, y=Sequence_Distance)) + facet_wrap(~ Type.x) + 
  xlab("Ohnolog in Same Class") + ylab("Sequence Distance (Poisson Model)") + ggtitle("C-Terminus")
 
summary(aov(sameClass ~ Sequence_Distance, FmegTabC[which(FmegTabC$Type.x == "WGD1"),])) # 0.761
summary(aov(sameClass ~ Sequence_Distance, FmegTabC[which(FmegTabC$Type.x == "WGD2"),])) # 2.22e-05
summary(aov(sameClass ~ Sequence_Distance, FmegTabC[which(FmegTabC$Type.x == "WGD3"),])) # 0.00101

# How about first 20aas?

for(fasta in fastas){
  tmpFasta = read.fasta(fasta, seqtype = "AA", as.string = T)
  tmpFasta1 = tmpFasta[1]
  tmpFasta2 = tmpFasta[2]
  
  tmpFasta1seq = substr(unlist(getSequence(tmpFasta1, as.string = T)), 
                       start = 1, stop = 20)
  tmpFasta2seq = substr(unlist(getSequence(tmpFasta2, as.string = T)), 
                       start = 1, stop = 20)

  name = paste("/ParameciumPCP/sequences/WGD-Nterm-20/", getName(tmpFasta1), "_", getName(tmpFasta2), ".fasta", sep = "")
  write.fasta(sequences = list(tmpFasta1seq, tmpFasta2seq), names = getName(tmpFasta), file.out = name)
}

megFilesNtwenty = list.files("/ParameciumPCP/sequences/WGD-Nterm-dists/", full.names = T)
megFilesNtwenty = megFilesNtwenty[grep(pattern = "*[0-9].meg", x = megFilesNtwenty)]
FmegTabN = data.frame(matrix(ncol = 3, nrow = 1))
colnames(FmegTabN) = c("GeneA", "GeneB", "Sequence_Distance")

for(meg in megFilesNtwenty){
  megTab = read.table(meg, skip = 37)[2]  
  colnames(megTab) = "Sequence_Distance"
  meg = gsub("N_", "", meg)
  megTab$GeneA = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[1]
  megTab$GeneB = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[2]
  FmegTabN = rbind(FmegTabN, megTab)
}
FmegTabN = FmegTabN[2:nrow(FmegTabN),]
wgdDFn = wgdDF[,c(1:8,10:11)]
FmegTabN = FmegTabN %>% right_join(wgdDFn, by=c("GeneA","GeneB"))


#-----------------------------------------------------------------------------------------
### Diff Identification (gene family)
load("/ParameciumPCP/RData/ptetPCP_superANDunsuper.RData")
library(stringr)

ptetPCP = read.csv(file = "/ParameciumPCP/tables/ptetPCP.csv", header=T) 
ptetPCP = ptetPCP[grep(pattern = "PTET*", x = ptetPCP$Accession),]   
dim(ptetPCP)  # 11856 proteins X 248 columns (34 proteins removed)

cdhit = read.csv("/ParameciumPCP/tables/CDhit/CDHIT_clusters.csv", header=T)
ptetWGD = read.table("/Paramecium_Paraorthologs/ptetraurelia_mac_51_annotation_v2.0.WGD.tree", header=T) 
ptetWGD$NB = NULL

for(q in 1:ncol(ptetWGD)){  # make dots into NA
  ptetWGD[,q] = sub("^[.]$", NA, x = ptetWGD[,q])
}

ptetWGDfam = data.frame(matrix(ncol=13, nrow=1))
colnames(ptetWGDfam) = c("OhnologFam", "Ngenes", "NgenesIDed", "NgenesClass", "NdiffClass", colnames(ptetWGD))

for(i in 1:nrow(ptetWGD)){
  ptetWGDi = ptetWGD[i,]
  
  if(as.numeric(rowSums(is.na(ptetWGDi))) == 0){
    ptetWGDiGenes = as.character(ptetWGDi)
  } else{
    ptetWGDiGenes = as.character(ptetWGDi)[-which(is.na(as.character(ptetWGDi)))]
  }

  
  tmpWGDi = data.frame(matrix(ncol=13, nrow=1))
  colnames(tmpWGDi) = c("OhnologFam", "Ngenes", "NgenesIDed", "NgenesClass", "NdiffClass", colnames(ptetWGD))
  
  tmpWGDi$OhnologFam = rownames(ptetWGDi)
  tmpWGDi$Ngenes = length(ptetWGDiGenes)
  
  ptetPCP[which(ptetPCP$Accession %in% ptetWGDiGenes),1:3]
  cdhit[grep(paste(ptetWGDiGenes,collapse = "|"), cdhit$Accession.s.),]
  nrow(cdhit[grep(paste(ptetWGDiGenes,collapse = "|"), cdhit$Accession.s.),])
  
  if(nrow(cdhit[grep(paste(ptetWGDiGenes,collapse = "|"), cdhit$Accession.s.),]) == length(ptetWGDiGenes)){  # X genes in X cdhit rows
    tmpWGDi$NgenesIDed = nrow(ptetPCP[which(ptetPCP$Accession %in% ptetWGDiGenes),])
  } else{  # X genes in non-X rows... complicated. count nrow(ptetPCP) for IDed genes. then add genes in CDhit clusters only if they're ohnologs
    tmpWGDi$NgenesIDed = length(grep(gsub("; ", "|", paste(cdhit[grep(paste(ptetWGDiGenes,collapse = "|"), cdhit$Accession.s.),'Accession.s.'], collapse = "|")), 
                                     ptetWGDiGenes))
  }
  tmpWGDi$NgenesClass = nrow(fData(preds2)[which(fData(preds2)$'Accession' %in% ptetPCP[which(ptetPCP$Accession %in% ptetWGDiGenes),'Accession']),])
  tmpWGDi$NdiffClass = length(unique(fData(preds2)[which(fData(preds2)$'Accession' %in% ptetPCP[which(ptetPCP$Accession %in% ptetWGDiGenes),'Accession']),'svm']))
  tmpWGDi[,colnames(ptetWGD)] = as.character(ptetWGDi)
  ptetWGDfam = rbind(ptetWGDfam, tmpWGDi)
}

ptetWGDfam = ptetWGDfam[c(2:nrow(ptetWGDfam)),]
ptetWGDfam$IDratio = ptetWGDfam$NgenesIDed/ptetWGDfam$Ngenes
ptetWGDfam$ClassProp = ifelse(test = ptetWGDfam$NdiffClass == 1, 
                              yes = 1, 
                              no = (ptetWGDfam$NdiffClass)/(ptetWGDfam$NgenesClass))

ptet_singleCopy
ptetWGDfamMulti = ptetWGDfam[which(ptetWGDfam$Ngenes > 1),]

nrow(ptetWGDfamMulti)
length(which(ptetWGDfamMulti$NgenesIDed > 0))
summary(ptetWGDfamMulti[which(ptetWGDfamMulti$IDratio > 0),'IDratio'])


length(which(ptetWGDfamMulti$IDratio == 1))  # how many have all IDed
length(which(ptetWGDfamMulti$IDratio == 0))  # how many have all IDed

FmegTab$Genes = apply(FmegTab[,c("GeneA", "GeneB")], 1, paste, collapse = "|")  # for easier grepping

#save.image("/ParameciumPCP/RData/WGD-All.Rdata")
load("/ParameciumPCP/RData/WGD-All.Rdata")

# Make ptet gene family fastas
library(seqinr)
ptProts

ptetWGDfamMulti$OhnologFam = as.numeric(ptetWGDfamMulti$OhnologFam)

for(rr in ptetWGDfamMulti$OhnologFam){
  # write.fasta(sequences = getSequence(ptProts[ptetWGDfamMulti[rr,5:12][!is.na(ptetWGDfamMulti[rr,5:12])]]), 
  #             names = getName(ptProts[ptetWGDfamMulti[rr,5:12][!is.na(ptetWGDfamMulti[rr,5:12])]]), 
  #             file.out = paste("/ParameciumPCP/sequences/WGD-families/", rr, ".fasta", sep=""))
  # 
  write.fasta(sequences = getSequence(ptProts[ptetWGDfamMulti[which(ptetWGDfamMulti$OhnologFam == rr),5:12][!is.na(ptetWGDfamMulti[which(ptetWGDfamMulti$OhnologFam == rr),5:12])]], as.string = T), 
    names = getName(ptProts[ptetWGDfamMulti[which(ptetWGDfamMulti$OhnologFam == rr),5:12][!is.na(ptetWGDfamMulti[which(ptetWGDfamMulti$OhnologFam == rr),5:12])]]), 
    file.out = paste("/ParameciumPCP/sequences/WGD-families/", rr, ".fasta", sep=""))
  
}
  # align and pdist(poisson) in bash
  # take Sum-of-pairs for MSA score

megFilesNfam = list.files("/ParameciumPCP/sequences/WGD-families-distances/", full.names = T)
megFilesNfam = megFilesNfam[grep(pattern = "*[0-9].csv", x = megFilesNfam)]
FmegTabF = data.frame(matrix(ncol = 2, nrow = 1))
colnames(FmegTabF) = c("OhnologFam", "Sequence_Distance")

for(meg in megFilesNfam){
  megTab = read.csv(meg, header = T)
  colnames(megTab) = "Sequence_Distance"
  megTab$OhnologFam = gsub(".csv", "", basename(meg))
  FmegTabF = rbind(FmegTabF, megTab)
}
FmegTabF = FmegTabF[2:nrow(FmegTabF),]
FmegTabF = merge(FmegTabF, ptetWGDfamMulti, by="OhnologFam")

plot1 = ggplot(data = FmegTabF, aes(x = IDratio, y = Sequence_Distance, fill = IDratio)) +
  geom_jitter() +
  geom_smooth(method = "lm") + xlab("Proportion of Ohnologs Identified") + ylab("Mean Sequence Distance (Poisson Model)") + 
  labs(title = "Sequence Distance vs. Identification Ratio", x = "ID Ratio", y = "Sequence Distance")

plot2 = ggplot(data = FmegTabF, aes(x = ClassProp, y = Sequence_Distance, fill = ClassProp)) +
  geom_jitter() +
  geom_smooth(method = "lm") + xlab("Unique Compartments per Identified Protein") + ylab("Mean Sequence Distance (Poisson Model)") + 
  labs(title = "Sequence Distance vs. Classification Ratio", x = "Classification Ratio", y = "Sequence Distance")

svg("/ParameciumPCP/plots/WGD-IDvsClass.svg")
grid.arrange(plot1, plot2)
dev.off()

summary(lm(IDratio ~ Sequence_Distance, FmegTabF))  #  <2e-16
summary(lm(ClassProp ~ Sequence_Distance, FmegTabF))  #  0.0341

summary(aov(IDratio ~ Sequence_Distance, FmegTabF))  #  <2e-16
summary(aov(ClassProp ~ Sequence_Distance, FmegTabF))  # 0.0341

cor(FmegTabF$IDratio, FmegTabF$Sequence_Distance, method = "spearman")  # -0.3112867
cor(FmegTabF$ClassProp, FmegTabF$Sequence_Distance, use="complete.obs", method = "spearman")  # -0.09158769



#save.image("/ParameciumPCP/RData/WGD-All.Rdata")
load("/ParameciumPCP/RData/WGD-All.Rdata")

# most diverse class
FmegTabF[which(FmegTabF$ClassProp > 0.25),]
# OhnologFam Sequence_Distance Ngenes NgenesIDed NgenesClass NdiffClass   ID1                ID2
# 10003        0.14763125      3          3           3          1        PTET.51.1.P0600183               <NA>
# 10004        0.04091135      4          3           3          1        PTET.51.1.P0600185 PTET.51.1.P0520082
# ID3                  ID4                 ID5  ID6  ID7  ID8    IDratio ClassProp
# PTET.51.1.P1370066   PTET.51.1.P0270092  <NA> <NA> <NA> <NA>    1.00   0.3333333
# PTET.51.1.P1370069   PTET.51.1.P0270090  <NA> <NA> <NA> <NA>    0.75   0.3333333






#--------------------------------------------------------------------
#Yeast data
library(pRolocdata)
data("yeast2018")
WGDyeast = read.csv("/ParameciumPCP/tables/WGD/yeastWGD.csv", header=T)
WGDyeast = WGDyeast[which(complete.cases(WGDyeast)),]
fdY = fData(yeast2018)

WGDyeastSVM = merge(merge(WGDyeast, fData(yeast2018)[which(as.character(fData(yeast2018)$'Accession') %in% WGDyeast$GeneA_Uniprot),c('Accession', 'svm', 'predicted.location')],
      by.x = "GeneA_Uniprot", by.y = "Accession"),
      fData(yeast2018)[which(as.character(fData(yeast2018)$'Accession') %in% WGDyeast$GeneB_Uniprot),c('Accession', 'svm', 'predicted.location')], 
      by.x = "GeneB_Uniprot", by.y = "Accession")
colnames(WGDyeastSVM)[10:13] = c("SVM_GeneA_Class", "SVM_GeneA_Pred", "SVM_GeneB_Class", "SVM_GeneB_Pred")

YchiScores = data.frame(matrix(ncol=5, nrow=1))
colnames(YchiScores) = c("Type", "Status", "GeneA", "GeneB", "EuclideanDistance")

euclid = function(a,b) sqrt(sum((a-b)^2))

for(q in 1:nrow(WGDyeast)){
  WGDyeastQ = WGDyeast[q,]
  tmpChiDF = data.frame(matrix(ncol=5, nrow=1))
  colnames(tmpChiDF) = c("Type", "Status", "GeneA", "GeneB", "EuclideanDistance")
  tmpChiDF$GeneA = WGDyeastQ$GeneA_Uniprot
  tmpChiDF$GeneB = WGDyeastQ$GeneB_Uniprot
  tmpRows = which(rownames(exprs(yeast2018)) %in% c(WGDyeastQ$GeneA_Uniprot, WGDyeastQ$GeneB_Uniprot))
  
  if(length(tmpRows) == 0){ # both copies are missing from the data
    tmpChiDF$Status = "Both Absent"
    tmpChiDF$EuclideanDistance = NA
  }
  if(length(tmpRows) == 1){ # one copy is present in the data
    tmpChiDF$Status = "One Absent"
    tmpChiDF$EuclideanDistance = NA
  }
  if(length(tmpRows) == 2){ # both copies are present in the data
    tmpChiDF$Status = "Both Present"
    tmpChiDF$EuclideanDistance = euclid(2^(exprs(yeast2018)[tmpRows,1]),
                                        2^(exprs(yeast2018)[tmpRows,2]))
  }
  YchiScores = rbind(YchiScores, tmpChiDF)
}
YchiScores = YchiScores[-1,]
YchiScores$Type = rep("WGD", nrow(YchiScores))


# one present one absent
YchiScores[YchiScores$Status == "One Absent",]

# Most divergent
YchiScores[which(YchiScores$EuclideanDistance == max(YchiScores$EuclideanDistance, na.rm = T)),]
plotDist(yeast2018[as.character(YchiScores[which(YchiScores$EuclideanDistance == max(YchiScores$EuclideanDistance, na.rm = T)),c('GeneA', 'GeneB')])], 
         pcol = c("red", "black"))
YchiScores[which(YchiScores$EuclideanDistance >  0.2 ),]

YchiScores[which(YchiScores$EuclideanDistance == min(YchiScores$EuclideanDistance, na.rm = T)),c('GeneA', 'GeneB')]
plotDist(yeast2018[as.character(YchiScores[which(YchiScores$EuclideanDistance == min(YchiScores$EuclideanDistance, na.rm = T)),c('GeneA', 'GeneB')])], 
         pcol = c("red", "black"))

# relocalization
WGDyeastSVM = cbind(WGDyeastSVM, WGDyeastSVM %>% summarise(relocalClass = SVM_GeneA_Class != SVM_GeneB_Class, 
                                                           relocalPred = SVM_GeneA_Pred != SVM_GeneB_Pred))

write.csv(WGDyeastSVM, 
          "/ParameciumPCP/tables/WGD/relocalization-yeast.csv", quote = F, row.names = F)
write.csv(data.frame(table(WGDyeastSVM$SVM_GeneA, WGDyeastSVM$SVM_GeneB))[which(data.frame(table(WGDyeastSVM$SVM_GeneA, WGDyeastSVM$SVM_GeneB))[,3] != 0),], 
          "/ParameciumPCP/tables/WGD/relocalization-yeastSummary.csv", quote = F, row.names = F)

WGDyeastSVM$SVM_GeneA_Class = gsub(" ", "\n", WGDyeastSVM$SVM_GeneA_Class)
WGDyeastSVM$SVM_GeneA_Pred = gsub(" ", "\n", WGDyeastSVM$SVM_GeneA_Pred)

WGDyeastSVM$SVM_GeneB_Class = gsub(" ", "\n", WGDyeastSVM$SVM_GeneB_Class)
WGDyeastSVM$SVM_GeneB_Pred = gsub(" ", "\n", WGDyeastSVM$SVM_GeneB_Pred)

ggplot(WGDyeastSVM) + geom_count(aes(x=SVM_GeneA_Class, y=SVM_GeneB_Class, color = ..n..)) + scale_size(range=c(0,10)) + 
  xlab("Classified Compartment (Gene A)") + ylab("Classified Compartment (Gene B)") + theme(legend.position = "none")

ggplot(WGDyeastSVM) + geom_boxplot(aes(x=relocalClass, y=AAidRatio))
ggplot(WGDyeastSVM) + geom_boxplot(aes(x=relocalClass, y=LengthRatio))
summary(aov(relocalClass ~ AAidRatio + LengthRatio, WGDyeastSVM))

WGDyeastSame = WGDyeastSVM[which(WGDyeastSVM$SVM_GeneA == WGDyeastSVM$SVM_GeneB),]
nrow(WGDyeastSame)
WGDyeastDiff = WGDyeastSVM[which(WGDyeastSVM$SVM_GeneA != WGDyeastSVM$SVM_GeneB),]
nrow(WGDyeastDiff)

sort(table(fData(yeast2018)$'svm'), decreasing = T)

# Get yeast seqs
library(seqinr)
yeastProt = read.fasta("/ParameciumPCP/sequences/yeast-proteome.fa", as.string = T, seqtype = "AA")
names(yeastProt) = as.character(lapply(sapply(getName(yeastProt), strsplit, "\\|"), "[", 2))

which(names(yeastProt) %in% c(WGDyeastSVM$GeneB_Uniprot, WGDyeastSVM$GeneA_Uniprot))

for(y in 1:nrow(WGDyeastSVM)){
  write.fasta(sequences = getSequence(yeastProt[c(WGDyeastSVM$GeneB_Uniprot[y], WGDyeastSVM$GeneA_Uniprot[y])], as.string = T), 
              names = names(yeastProt[c(WGDyeastSVM$GeneB_Uniprot[y], WGDyeastSVM$GeneA_Uniprot[y])]), 
              file.out = paste("/ParameciumPCP/sequences/WGD-yeast/", 
                               WGDyeastSVM$GeneB_Uniprot[y], "_", WGDyeastSVM$GeneA_Uniprot[y], ".fasta", sep = ""))
}


















#
#
#
#
#
#





##
##
## Old
#-------------------------------------------------------------------------------------------------------------------------
##### Gene Duplication
load("/ParameciumPCP/RData/PCpulldown-data.RData")
ptetWGD = read.table("/Paramecium_Paraorthologs/ptetraurelia_mac_51_annotation_v2.0.WGD.tree", header=T) 
ptetWGD$NB = NULL

for(q in 1:ncol(ptetWGD)){  # make dots into NA
  ptetWGD[,q] = sub("^[.]$", NA, x = ptetWGD[,q])
}
head(ptetWGD)

source("/ParameciumPCP/scripts/source_wgdPairs.R")
source("/ParameciumPCP/scripts/fun_WGDtestMakeBind.R")

wgdDF = data.frame(matrix(nrow=1, ncol=4))  # initiate final DF
wgdDF[1,] = c("GeneA", "GeneB", "Type", "PPSS")

for(i in 1:nrow(ptetWGD)){
  ptetWGD_i = ptetWGD[i,]
  nProts = length(which(is.na(ptetWGD_i) == F))
  
  if(nProts > 1){
    wgdDF = fun_WGDtestMakeBind(wgd1_1, ptetWGD_i, lPPSS_noImpute, wgd = "WGD1")
    wgdDF = fun_WGDtestMakeBind(wgd1_2, ptetWGD_i, lPPSS_noImpute, wgd = "WGD1")
    wgdDF = fun_WGDtestMakeBind(wgd1_3, ptetWGD_i, lPPSS_noImpute, wgd = "WGD1")
    wgdDF = fun_WGDtestMakeBind(wgd1_4, ptetWGD_i, lPPSS_noImpute, wgd = "WGD1")
    
    wgdDF = fun_WGDtestMakeBind(wgd2_1, ptetWGD_i, lPPSS_noImpute, wgd = "WGD2")
    wgdDF = fun_WGDtestMakeBind(wgd2_2, ptetWGD_i, lPPSS_noImpute, wgd = "WGD2")
    wgdDF = fun_WGDtestMakeBind(wgd2_3, ptetWGD_i, lPPSS_noImpute, wgd = "WGD2")
    wgdDF = fun_WGDtestMakeBind(wgd2_4, ptetWGD_i, lPPSS_noImpute, wgd = "WGD2")
    wgdDF = fun_WGDtestMakeBind(wgd2_5, ptetWGD_i, lPPSS_noImpute, wgd = "WGD2")
    wgdDF = fun_WGDtestMakeBind(wgd2_6, ptetWGD_i, lPPSS_noImpute, wgd = "WGD2")
    wgdDF = fun_WGDtestMakeBind(wgd2_7, ptetWGD_i, lPPSS_noImpute, wgd = "WGD2")
    wgdDF = fun_WGDtestMakeBind(wgd2_8, ptetWGD_i, lPPSS_noImpute, wgd = "WGD2")
    
    wgdDF = fun_WGDtestMakeBind(wgd3_1, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_2, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_3, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_4, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_5, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_6, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_7, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_8, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_9, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_10, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_11, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_12, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_13, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_14, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_15, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
    wgdDF = fun_WGDtestMakeBind(wgd3_16, ptetWGD_i, lPPSS_noImpute, wgd = "WGD3")
  }
}
wgdDF = wgdDF[complete.cases(wgdDF),] # Clean it up
colnames(wgdDF) = wgdDF[1,]
wgdDF = wgdDF[2:nrow(wgdDF),]  # 440 rows

table(wgdDF$Type)  # 1: 195 pairs; 2: 187 pairs; 3: 58 pairs
head(wgdDF[order(wgdDF$PPSS),])
tail(wgdDF[order(wgdDF$PPSS),])

write.csv(wgdDF, "/ParameciumPCP/tables/ptetPCP_WGD-PPSS.csv", quote = F, row.names = F)

# Make random DF to compare
ranDF = data.frame(matrix(nrow=1, ncol=4))
colnames(ranDF) = c("GeneA", "GeneB", "Type", "PPSS")
set.seed(42)
ranRows1 = round(runif(58, min = 1, length(lPPSS_noImpute)))  # smallest group is 58 (WGD3)
set.seed(666)
ranRows2 = round(runif(58, min = 1, length(lPPSS_noImpute)))

for(f in 1:58){
  tmpRanDF = data.frame(matrix(nrow=1, ncol=4))
  colnames(tmpRanDF) = c("GeneA", "GeneB", "Type", "PPSS")
  tmpRanDF$GeneA = names(lPPSS_noImpute)[f]
  tmpRanDF$GeneB = names(lPPSS_noImpute[[f]])[f]
  tmpRanDF$Type = "Random"
  tmpRanDF$PPSS = lPPSS_noImpute[[ranRows1[f]]][ranRows2[f]]
  ranDF = rbind(ranDF, tmpRanDF)
}
ranDF = ranDF[2:nrow(ranDF),]
pssWGDran = rbind(wgdDF, ranDF)

# plot
pssWGDran$PPSS = as.numeric(pssWGDran$PPSS)
ggplot(pssWGDran) + geom_violin(aes(x=Type, y=PPSS)) + xlab("Ohnolog Type (or Random Pair)") + ylab("Protein Similarity Score")

for(q in 1:nrow(wgdDF)){                               # all pairwise distance plots- both SVG and PNG (switch out 2 things)
  if(length(which(row.names(fData(preds2)) %in% as.character(wgdDF[q,1:2]))) == 2){
    png(paste("/ParameciumPCP/plots/distProf-WGD/", as.character(wgdDF[q,1]), "_", 
              as.character(wgdDF[q,2]), "_hist.png", sep=""))
    plotDist(preds2[as.character(pssWGDran[q,1:2])], pcol = c("black", "red"), ylim = c(0,1))
    dev.off()  # shut off plot saving
  }
}

# stats- different from random but not each other
wilcox.test(x=pssWGDran[which(pssWGDran$Type == "WGD1"),"PPSS"],
            y=pssWGDran[which(pssWGDran$Type == "Random"),"PPSS"])  # p-value = 2.2e-16
wilcox.test(x=pssWGDran[which(pssWGDran$Type == "WGD1"),"PPSS"],
            y=pssWGDran[which(pssWGDran$Type == "WGD2"),"PPSS"])  # p-value = .6018
wilcox.test(x=pssWGDran[which(pssWGDran$Type == "WGD1"),"PPSS"],
            y=pssWGDran[which(pssWGDran$Type == "WGD3"),"PPSS"])  # p-value = .9405
wilcox.test(x=pssWGDran[which(pssWGDran$Type == "WGD2"),"PPSS"],
            y=pssWGDran[which(pssWGDran$Type == "WGD3"),"PPSS"])  # p-value = .7593



# Sequence vs PPSS
megFiles = list.files("/ParameciumPCP/megacc/ptetPCP_WGD-pdist/", full.names = T)
megFiles = megFiles[grep(pattern = "*[0-9].meg", x = megFiles)]
FmegTab = data.frame(matrix(ncol = 3, nrow = 1))
colnames(FmegTab) = c("GeneA", "GeneB", "Sequence_Distance")

for(meg in megFiles){
  megTab = read.table(meg, skip = 37)[2]  # .meg files with 2 seqs have 37 lines of unimportance
  colnames(megTab) = "Sequence_Distance"
  megTab$GeneA = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[1]
  megTab$GeneB = unlist(strsplit(x = substr(x = basename(meg), start = 1, 37), "_"))[2]
  FmegTab = rbind(FmegTab, megTab)
}
FmegTab = FmegTab[2:nrow(FmegTab),]
pDistPPSSsvm = mostDivgernt_all %>% right_join(FmegTab, by=c("GeneA","GeneB"))
pDistPPSSsvm = pDistPPSSsvm[which(complete.cases(pDistPPSSsvm)),]

write.csv(pDistPPSSsvm, "/ParameciumPCP/tables/ptetPCP-WGD-PPSS-Pdist.csv", quote = F, row.names = F)

pDistPPSSsvm %>% group_by(Type) %>% summarise(range = range(Sequence_Distance))
pDistPPSSsvm %>% group_by(Type) %>% summarise(mean = mean(Sequence_Distance))
pDistPPSSsvm %>% group_by(Type) %>% summarise(median = median(Sequence_Distance))
pDistPPSSsvm %>% group_by(Type) %>% summarise(b = sum(Sequence_Distance > 0.1))

ggplot(pDistPPSSsvm, aes(x=Sequence_Distance, y=PPSS)) + 
  geom_point() + geom_smooth(method = "lm", col = "black") + 
  xlim(0,1.5) + ylim(0,1) + xlab("Sequence Divergence") + ylab("Protein Profile Similarity Score") + 
  scale_x_log10() + scale_y_log10()

ggplot(pDistPPSSsvm, aes(x=Sequence_Distance, y=PPSS)) + 
  geom_point() + geom_smooth(method = "lm", col = "black") + 
  xlim(0,1.5) + ylim(0,1) + xlab("Sequence Divergence") + ylab("Protein Profile Similarity Score") + 
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ Type)

WGD1.lm_log = lm(log(PPSS) ~ log(Sequence_Distance), data=pDistPPSSsvm[which(pDistPPSSsvm$Type == "WGD1"),]) 
summary(WGD1.lm_log) # p-value: 0.015; Adjusted R-squared:  0.0248 
WGD2.lm_log = lm(log(PPSS) ~ log(Sequence_Distance), data=pDistPPSSsvm[which(pDistPPSSsvm$Type == "WGD2"),]) 
summary(WGD2.lm_log) # p-value: p-value: 2.01e-05; Adjusted R-squared:  0.08893 
WGD3.lm_log = lm(log(PPSS) ~ log(Sequence_Distance), data=pDistPPSSsvm[which(pDistPPSSsvm$Type == "WGD3"),]) 
summary(WGD3.lm_log) # p-value: p-value: 0.002684; Adjusted R-squared:  0.1347

plot(fitted(WGD1.lm_log), resid(WGD1.lm_log))  # residuals
abline(0,0)
plot(fitted(WGD2.lm_log), resid(WGD2.lm_log))
abline(0,0)
plot(fitted(WGD3.lm_log), resid(WGD3.lm_log))
abline(0,0)

WGD.lm = lm(log(PPSS) ~ log(Sequence_Distance) , data=pDistPPSSsvm)
summary(WGD.lm)
