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


#-------------------------------------------------------
# Sélection des réplicats associés aux conditions à comparer
#-------------------------------------------------------
# Je compare S1 vs S2
replicates1 <- c("HCA.3","HCA.4","HCA.5")
replicates2 <- c("HCA.9","HCA.10","HCA.11")

# J'extrait les colonnes de mon tableau de comptage 
# qui correspondent aux conditions S1 vs S2
diffParam = cbind(rawCountTable[,replicates1], rawCountTable[,replicates2])


#-------------------------------------------------------
# Formatage des données pour l'utilisation de DESeq
#-------------------------------------------------------
# for DESeq
conds = c(rep("Cond1", length(replicates1)), rep("Cond2", length(replicates2)))

# create metadata
expDesign = data.frame(row.names = colnames(rawCountTable[,c(replicates1, replicates2)]), 
                       condition = conds,
                       libType   = rep("RNAseq", length(c(replicates1, replicates2))))

# create a countDataSet
cds = newCountDataSet(rawCountTable[,c(replicates1, replicates2)], expDesign$condition)


#-------------------------------------------------------
# Analyse différentielle par la méthode DESeq
#-------------------------------------------------------
## le détail de la méthode est disponible en tapant help(DESeq)

# Estimation des facteurs de normalisations
cds = estimateSizeFactors(cds)
# Inspect size factors
sizeFactors(cds)

# Estimation du paramètres de dispersion de la loi BN
cds = estimateDispersions(cds)
str(fitInfo(cds))

# differential analysis (ratio A/B is performed) 
res = nbinomTest(cds, condA = "Cond2", condB = "Cond1")

#-------------------------------------------------------
# Analyse des résultats de l'analyse différentielle
#-------------------------------------------------------

str(res)

## ENONCE TP 3
## A rendre pour la séance 4 (reproduire la figure "supplémentary material numéro 3" de l'article)
## Pour chaque comparaison +Fe / -Fe: représenter le nombre de gènes différentiellement sur-exprimés 
## et le nombre de gènes différentiellement sous-exprimés pour chacune des comparaisons



