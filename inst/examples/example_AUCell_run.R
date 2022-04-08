# This example is run using a fake expression matrix.
# Therefore, the output will be meaningless.

############# Fake expression matrix #############
set.seed(123)
exprMatrix <- matrix(data=sample(c(rep(0, 5000), sample(1:3, 5000, replace=TRUE))),
                     nrow=20, 
                     dimnames=list(paste("Gene", 1:20, sep=""), 
                                   paste("Cell", 1:500, sep="")))
exprMatrix <- as(exprMatrix, "dgCMatrix")

# In this example we use two gene sets: 10 and 5 random genes
# (see other formatting examples at the end)
fewGenes <- sample(rownames(exprMatrix), 10)
otherGenes <- sample(rownames(exprMatrix), 5)

geneSets <- list(geneSet1=fewGenes,
                 geneSet2=otherGenes)
geneSets

# Calculate AUCell score for the genes in the sets
# To be able to run this fake example (which contain only 20 genes),
# we use aucMaxRank=5 (top 25% of the genes in the ranking)
cells_AUC <- AUCell_run(exprMatrix, geneSets, aucMaxRank=5)


## To run in paralell:
# cells_AUC <- AUCell_run(exprMatrix, geneSets, aucMaxRank=5, 
#                   BPPARAM=BiocParallel::MulticoreParam(5))



# Format of the output:
cells_AUC

# To subset & access the AUC slot (as matrix):
cells_AUC[1:2,]
cells_AUC[,3:4]
getAUC(cells_AUC)[,1:5]


# These methods are also available:
dim(cells_AUC)
nrow(cells_AUC)
ncol(cells_AUC)
colnames(cells_AUC)[1:4]
rownames(cells_AUC)

#########################################################
# Alternatives for the input of gene sets:

# a) Character vector (i.e. only one gene-set)
# It will take the default name 'geneSet'
fewGenes
test <- AUCell_run(exprMatrix, fewGenes, aucMaxRank=5)

# b) List
geneSets <- list(geneSet1=fewGenes,
                 geneSet2=otherGenes)
geneSets
test <- AUCell_run(exprMatrix, fewGenes, aucMaxRank=5)

# c) GeneSet object (from GSEABase)
library(GSEABase)
geneSetOne <- GeneSet(fewGenes, setName="geneSetOne")
geneSetOne
test <- AUCell_run(exprMatrix, fewGenes, aucMaxRank=5)


# d) GeneSetCollection object (from GSEABase)
geneSetTwo <- GeneSet(otherGenes, setName="geneSetTwo")
geneSets <- GeneSetCollection(geneSetOne, geneSetTwo)
geneSets
test <- AUCell_run(exprMatrix, fewGenes, aucMaxRank=5)

