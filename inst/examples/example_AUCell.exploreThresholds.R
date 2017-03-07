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

######### Previous steps in the workflow #########
# Step 1.
cells_rankings <- AUCell.buildRankings(exprMatrix, plotStats=FALSE)

# Step 2.
# (Gene sets: random genes)
geneSets <- list(geneSet1=sample(rownames(exprMatrix), 10),
                 geneSet2=sample(rownames(exprMatrix), 5))
cells_AUC <- AUCell.calcAUC(geneSets, cells_rankings, aucMaxRank=5)
##################################################

############## Step 3: Assign cells ##############

# 1. Plot histograms and obtain some pre-computed thresholds
# (this example is only meant to show the interface/arguments of the function,
# see the vignette for meaningful examples)
par(mfrow=c(1,2)) # Plot is divided into one row and two columns
thresholds <- AUCell.exploreThresholds(cells_AUC, seed=123, plotHist=TRUE)
thresholds$geneSet1$aucThr

# 2. Obtain cells over a given threshold:
names(which(cells_AUC@AUC["geneSet1",] > 0.19))

# Alternative: assign cells according to the 'automatic' threshold
cells_assignment <- AUCell.exploreThresholds(cells_AUC, seed=123,
                                             plotHist=FALSE, assignCells=TRUE)
# Cells assigned:
lapply(cells_assignment, function(x) x$assignment)
# Threshold applied:
sapply(cells_assignment, function(x) x$aucThr$selected)


