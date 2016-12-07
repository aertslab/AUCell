# This example is run using a fake expression matrix.
# Therefore, the output will be meaningless.

############# Fake expression matrix #############
set.seed(123)
exprMatrix <- matrix(data=sample(c(rep(0, 5000), sample(1:3, 5000, replace=TRUE))),
                     nrow=20)
rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
colnames(exprMatrix) <- paste("Cell", 1:500, sep="")
dim(exprMatrix)
##################################################

######### Previous step in the workflow ##########
# Step 1.
cells_rankings <- AUCell.buildRankings(exprMatrix, plotStats=FALSE)
##################################################

############## Step 2: Calculate AUC #############

# The gene sets should be provided as a list
# (In this example we use 10 and 5 random genes)
genes <- sample(rownames(exprMatrix), 10)
otherGenes <- sample(rownames(exprMatrix), 5)

geneSets <- list(geneSet1=genes,
                 geneSet2=otherGenes)
geneSets

# Calculate AUC with the rankings from Step 1.
# To be able to run this fake example (which contain only 20 genes),
# we use aucMaxRank=5 (top 25% of the genes in the ranking)
cells_AUC <- AUCell.calcAUC(geneSets, cells_rankings, aucMaxRank=5)

dim(cells_AUC)
head(cells_AUC)
