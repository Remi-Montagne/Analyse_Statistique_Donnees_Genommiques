rm(list=ls())

# DESeq package is used to identify differentially expressed genes
library(DESeq) 

help(DESeq)
#-------------------------------------------------------
# Importation des données O. tauri 
#-------------------------------------------------------
rawCountTable = read.table("mapping_rawdata_allGenes.txt", header = T, row.names = 1)[,-1]
geneLength    = read.table("mapping_rawdata_allGenes.txt", header = T, row.names = 1)[,1]
# --> Attention à garder à part l'information de longueur des gènes
sampleInfo    = read.table("sample_info.txt", header = T, row.names = 1)

replicates1 <- c("HCA.38","HCA.42","HCA.46")
replicates2 <- c("HCA.26","HCA.30","HCA.34")

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
res = nbinomTest(cds, condA = "Cond2", condB = "Cond1")
str(res)

subset <- res[res$pval <= 0.05,]
fc <- subset$foldChange

upreg <- fc[fc > 0]
downreg <- fc[fc < 0]
comptes <- c(length(upreg), length(downreg))
hist(comptes, col = c(2,3))
