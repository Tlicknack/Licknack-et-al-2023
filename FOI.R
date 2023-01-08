load("/ParameciumPCP/RData/ptetPCP_superAND unsuper.RData")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###  Contractile Vacuole Complex
cvc = read.csv("/ParameciumPCP/tables/FOI/contractileVacuoleComplex.csv", header=T)
cvc = cvc[as.character(cvc$ID) %in% fData(preds2)$'Accession',]
foi_cvc = FeaturesOfInterest(fnames = cvc$ID, 
                            description = "Contractile Vacuole Proteins", object = preds2)

plot2D(preds2, fcol = "markers")
highlightOnPlot(preds2, foi = foi_cvc, pch = 24, col="black", bg="gray")
plot2D(preds2, fcol = "markers", dims = c(1,3))
highlightOnPlot(preds2, foi = foi_cvc, pch = 24, col="black", bg="gray")


plotDist(preds2[cvc$ID], pcol = "red")

merge(cvc, fData(preds2)[,c('Accession', 'svm.pred')], by.x="ID", by.y="Accession")
merge(cvc, fData(preds2)[,c('Accession', 'svm')], by.x="ID", by.y="Accession")

#### Phagosomes
phago = read.csv("/ParameciumPCP/tables/FOI/phagosome.csv", header=T)
phago = phago[as.character(phago$ID) %in% fData(preds2)$'Accession',]
foi_phago = FeaturesOfInterest(fnames = phago$ID, 
                             description = "Food Vacuole Proteins", object = preds2)
plot2D(preds2, fcol = "markers")
highlightOnPlot(preds2, foi = foi_phago, pch = 24, col="black", bg="gray")

plotDist(preds2[phago$ID], pcol = "red")

merge(phago, fData(preds2)[,c('Accession', 'svm.pred')], by.x="ID", by.y="Accession")
merge(phago, fData(preds2)[,c('Accession', 'svm')], by.x="ID", by.y="Accession")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Sortilins

sorts = read.csv("/ParameciumPCP/tables/FOI/sortilins.csv", header=T)
sorts = sorts[as.character(sorts$Accession) %in% fData(preds2)$'Accession',]
foi_sorts = FeaturesOfInterest(fnames = sorts$Accession, 
                               description = "SOR1-16", object = preds2)
plot2D(preds2, fcol = "markers")
highlightOnPlot(preds2, foi = foi_sorts, pch = 24, col="black", bg="gray")

plotDist(preds2[sorts$Accession], pcol = "red")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###  Apicoplast

library(seqinr)
ptProts = read.fasta("/ParameciumPCP/sequences/ptet-proteome-UtoC.fa", seqtype = "AA", as.string = T)

apis = read.table("/ParameciumPCP/test/apicoplast.txt", header=F)

write.fasta(sequences = getSequence(ptProts[apis$V1], as.string = T), ### all protein sequences
            names = getName(ptProts[apis$V1]), 
            file.out = "/ParameciumPCP/sequences/apicoplast-pred-orthos.fasta")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### DNARepair

dnaRepair = read.csv("./ptetPCP_csv/FOI/ptetPCP_DNArepair.csv", header=T)
dnaRepair = dnaRepair[dnaRepair$Accession %in% fData(preds2)$'Accession',]
foi_dnaRepair = FeaturesOfInterest(fnames = dnaRepair$Accession, 
                                   description = "Annotated DNA Repair Proteins", object = preds2)
  # Plots
set.seed(42)
.tSne = plot2D(preds2, fcol = "svm.pred", method = "t-SNE")

highlightOnPlot(.tSne, foi = foi_dnaRepair, pch = 24, col="black", bg="yellow")
addLegend(preds2, fcol = "svm.pred", where = "topright", bty = "n", cex = .6)

plotDist(preds2[dnaRepair$Accession], fcol = NULL, pcol = "yellow2")

  # Properties
svmProps 

svmProps[which(svmProps$Accession %in% dnaRepair$Accession),c(1,13:14)]

plotDist(preds2["PTET.51.1.P0250220"], pcol = "red")  # ku70-2; MT3
plotDist(preds2["PTET.51.1.P0540024"], pcol = "red") # PtLigIV_1; Nuclei
plotDist(preds2["PTET.51.1.P1460025"], pcol = "red") # PtKu80-1; MT3
plotDist(preds2["PTET.51.1.P1510135"], pcol = "red") # PtKu80-2; Cytosol

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###  Rab.R
rabs = read.csv("./ptetPCP_csv/FOI/ptetPCP_FOI_Rab.csv", header = T)
rabs = rabs[which(rabs$Protein.Name %in% fData(preds2)$'Accession'),]
foi_rab = FeaturesOfInterest(fnames = rabs$Protein.Name, 
                             description = "Annotated Rab Proteins", object = preds2)

# Plots
set.seed(42)
.tSne = plot2D(preds2, fcol = "svm.pred", method = "t-SNE")
highlightOnPlot(.tSne, foi = foi_rab, pch = 24, col="black", bg="orange")
addLegend(preds2, fcol = "svm.pred", where = "topright", bty = "n", cex = .6)

plotDist(preds2[rabs$Protein.Name], fcol = NULL, pcol = "orange")

# get fasta
library(seqinr)
ptProts = read.fasta("/Paramecium_Proteomes-Predicted/ptetraurelia_mac_51_annotation_v2.0.protein.fa", seqtype = "AA", as.string = T)
write.fasta(sequences = getSequence(ptProts[rabs$Protein.Name], as.string = T), names = getName(ptProts[rabs$Protein.Name]), file.out = "./ptetPCP_fasta/ptetPCP_rab.faa")

# properties
svmProps = read.csv("./ptetPCP_csv/ptetPCP_Properties.csv", header= T)
svmProps[which(svmProps$Accession %in% rabs$Protein.Name),c(1,5,8,9,11,13:15)]

table(svmProps[which(svmProps$Accession %in% rabs$Protein.Name),c(13)])
table(svmProps[which(svmProps$Accession %in% rabs$Protein.Name),c(11)])
svmProps[which(svmProps$Accession %in% rabs$Protein.Name),c(1,13)]
  # Outer Mitochondrial Membrane?
    # PTET.51.1.P0120313  PTET.51.1.P0870070  
plotDist(preds2["PTET.51.1.P0120313"], fcol = NULL, pcol = "red")
plotDist(preds2["PTET.51.1.P0870070"], fcol = NULL, pcol = "red")

  # Nuclear?
    # PTET.51.1.P0190040  PTET.51.1.P0620243
plotDist(preds2["PTET.51.1.P0190040"], fcol = NULL, pcol = "red")
plotDist(preds2["PTET.51.1.P0620243"], fcol = NULL, pcol = "red")

# make all plots
setwd("/ptetPCP/ptetPCP_plots/Rab/")
for(w in 1:nrow(rabs)){
  png(paste("dist_Rab_", as.character(rabs[w,1]), ".png", sep = "")) 
  plotDist(preds2[as.character(rabs[w,1])], fcol = NULL, pcol = "red", ylim = c(0,1))
  dev.off()  # shut off plot saving
}

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Cilia
marked_cilia = addMarkers(preds2, markers = "/ptetPCP/ptetPCP_csv/markersCilia.csv", 
                    mcol = "markersCilia", fcol = "Accession", verbose = T)
marked_ciliaM = addMarkers(preds2, markers = "/ptetPCP/ptetPCP_csv/markersCiliarMembrane.csv", 
                          mcol = "markersCiliarMembrane", fcol = "Accession", verbose = T)

# Plot projection
set.seed(42)
plot2D(marked_cilia, fcol = "markersCilia", method = "t-SNE")
set.seed(42)
plot2D(marked_ciliaM, fcol = "markersCiliarMembrane", method = "t-SNE")


# Plot distribution
plotDist(marked_cilia[which(fData(marked_cilia)$'markersCilia' == 'cilia')], pcol = "red")
plotDist(marked_ciliaM[which(fData(marked_ciliaM)$'markersCiliarMembrane' == 'Cilia-Membrane')], pcol = "red")

# Export Summary Table
ciliaTable = merge(data.frame(sort(table(fData(marked_cilia)[which(fData(marked_cilia)$'markersCilia' == "cilia"),'svm.pred']), 
                                   decreasing = T)), 
                   data.frame(sort(table(fData(marked_ciliaM)[which(fData(marked_ciliaM)$'markersCiliarMembrane' == "Cilia-Membrane"),'svm.pred']), 
                                   decreasing = T)), 
                   by = "Var1")
colnames(ciliaTable) = c("Predicted Compartment", "Cilia MS Predicted", "Cilia Membrane MS Predicted")

ciliaTable[nrow(ciliaTable)+1,] = c(".", colSums(ciliaTable[,2:3]))

write.csv(ciliaTable, "/ptetPCP/ptetPCP_csv/ciliaTable.csv", quote = F, row.names = F)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### Actin
fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0130204"])  # Act1_1 (comets around phagosomes)
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0130204"])
fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0850133"])  # Act1_5 ()
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0850133"])
fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P1470101"])  # Act4_1 (phagosomes, oral cavity, CLEAVAGE FURROW)
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P1470101"])
fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P1280118"])  # Act4_2 (phagosomes, oral cavity, CLEAVAGE FURROW)
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P1280118"])
fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P1560086"])  # Act5_1 (phagosomes, oral cavity)
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P1560086"])
fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0600063"])  # Act6_2 ()
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0600063"])
fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0620211"])  # Arp3_1 ()
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0620211"])
fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0240238"])  # Arp3_2 ()
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0240238"])

fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0190343"])  # Alp1_1 ()
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0190343"])

fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0630141"])  # Unnamed Arp (Arp8?)
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0630141"])

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
### PFUS
fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P1290118"])
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P1290118"])

fData(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0020067"])
exprs(preds2_ptetLOPITv3_scaled_all_impute_norm_zero_marked_ftmp4["PTET.51.1.P0020067"])
