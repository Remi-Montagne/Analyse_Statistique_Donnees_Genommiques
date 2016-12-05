rm(list=ls())
setwd("~/M2S1/Analyse_Stats_Donnees_Genomiques")

# DESeq package is used to identify differentially expressed genes
library(DESeq) 
library(DESeq2) 

help(DESeq)
#-------------------------------------------------------
# Importation des données O. tauri 
#-------------------------------------------------------
rawCountTable = read.table("mapping_rawdata_allGenes.txt", header = T, row.names = 1)[,-1]
geneLength    = read.table("mapping_rawdata_allGenes.txt", header = T, row.names = 1)[,1]
# --> Attention à garder à part l'information de longueur des gènes
sampleInfo = read.table("sample_info.txt", header = T, row.names = 1)

cond1 <- sort(unique(sampleInfo[sampleInfo$Condition == "Cond1", 5]))
cond2 <- sort(as.character(unique(sampleInfo[sampleInfo$Condition == "Cond2", 5])))
l <- length(cond2)
cond2 <- c(cond2[l], cond2[-l])
upreg <- c()
downreg <- c()

#calcul des genes up-regules et downregules en condition 1
for(i in seq(1, length(cond1) - 1, 2)) {
  c1 <- cond1[i]
  c2 <- cond1[i+1]
replicates1 <- row.names(sampleInfo[sampleInfo$Name == c2,])
replicates2 <- row.names(sampleInfo[sampleInfo$Name == c1,])

# combine data
diffParam = cbind(rawCountTable[,replicates1], rawCountTable[,replicates2])

# for DESeq
conds = c(rep("Cond1", length(replicates1)), rep("Cond2", length(replicates2)))

# create metadata
expDesign = data.frame(row.names = colnames(rawCountTable[,c(replicates1, replicates2)]), 
                       condition = conds,
                       libType   = rep("RNAseq", length(c(replicates1, replicates2))))

# create a countDataSet
cds = newCountDataSet(rawCountTable[,c(replicates1, replicates2)], expDesign$condition)

# Estimate normalization factors
cds = estimateSizeFactors(cds)

# Inspect size factors
sizeFactors(cds)

# To estimate the gene dispersions
cds = estimateDispersions(cds)
# To inspect the intermediate steps, a fitInfo object is stored, which contains the per 	genes estimate, the fitted curve
# and the values that will subsequently be used for inference.
str(fitInfo(cds))

# differential analysis (ratio A/B is performed) 
res = nbinomTest(cds, condA = "Cond1", condB = "Cond2")
res <- subset(res, ! is.na(res$padj))
str(res)
subs <- res[res$padj < 0.05,]


upreg <- c(upreg, length(subs[subs$log2FoldChange >= 1,1]))
downreg <- c(downreg, length(subs[subs$log2FoldChange <= -1,1]))

}

#calcul des genes up-regules et downregules en condition2
for(i in seq(1, length(cond2) - 1, 2)) {
  c1 <- cond2[i]
  c2 <- cond2[i+1]
  replicates1 <- row.names(sampleInfo[sampleInfo$Name == c2,])
  replicates2 <- row.names(sampleInfo[sampleInfo$Name == c1,])
  
  # combine data
  diffParam = cbind(rawCountTable[,replicates1], rawCountTable[,replicates2])
  
  # for DESeq
  conds = c(rep("Cond1", length(replicates1)), rep("Cond2", length(replicates2)))
  
  # create metadata
  expDesign = data.frame(row.names = colnames(rawCountTable[,c(replicates1, replicates2)]), 
                         condition = conds,
                         libType   = rep("RNAseq", length(c(replicates1, replicates2))))
  
  # create a countDataSet
  cds = newCountDataSet(rawCountTable[,c(replicates1, replicates2)], expDesign$condition)
  
  # Estimate normalization factors
  cds = estimateSizeFactors(cds)
  
  # Inspect size factors
  sizeFactors(cds)
  
  # To estimate the gene dispersions
  cds = estimateDispersions(cds)
  # To inspect the intermediate steps, a fitInfo object is stored, which contains the per 	genes estimate, the fitted curve
  # and the values that will subsequently be used for inference.
  str(fitInfo(cds))
  
  # differential analysis (ratio A/B is performed) 
  res = nbinomTest(cds, condA = "Cond1", condB = "Cond2")
  res <- subset(res, ! is.na(res$padj))
  str(res)
  subs <- res[res$padj < 0.05,]
  
  
  upreg <- c(upreg, length(subs[subs$log2FoldChange >= 1,1]))
  downreg <- c(downreg, length(subs[subs$log2FoldChange <= -1,1]))
  
}

comptes <- rbind(upreg, downreg)

barplot(comptes, col = c(2,3), beside = TRUE, names.arg = c("3h-light", "6h-light", "3h-dark", "6h-dark", "3h", "9h", "15h", "22h"))
abline(v = 12.5, lty = 2, lwd = 2)
legend("left", legend = c("up-regulation", "down-regulation"), fill = c(2, 3), cex = 1.1)
legend(x = 5, y = 600, legend = "condition1", box.lty = 0)
legend(x = 17, y = 600, legend = "condition2", box.lty = 0)
legend("left", legend = c("up-regulation", "down-regulation"), fill = c(2, 3), cex = 1.1)