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

############# Begining of the workflow ###########
# Step 1.
cells_rankings <- AUCell.buildRankings(exprMatrix, plotStats=FALSE)

# Step 2.
# (Gene set: 10 random genes)
genes <- sample(rownames(exprMatrix), 10)
geneSets <- list(geneSet1=genes)
# (aucMaxRank=5 to run with this fake example, it will return 'high' AUC values)
cells_AUC <- AUCell.calcAUC(geneSets, cells_rankings, aucMaxRank=5)
##################################################

# Plot histogram:
AUC.plot(cells_AUC@AUC["geneSet1",], gSetName="A fake gene-set", nBreaks=10)

