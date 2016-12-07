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

cells_rankings <- AUCell.buildRankings(exprMatrix, plotStats=TRUE)
dim(cells_rankings)
