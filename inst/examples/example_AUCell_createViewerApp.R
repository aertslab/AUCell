

######
# Fake run of AUCell
set.seed(123)
exprMatrix <- matrix(
  data=sample(c(rep(0, 5000), sample(1:3, 5000, replace=TRUE))),
  nrow=20, 
  dimnames=list(paste("Gene", 1:20, sep=""), 
                paste("Cell", 1:500, sep="")))
geneSets <- list(geneSet1=sample(rownames(exprMatrix), 10),
                 geneSet2=sample(rownames(exprMatrix), 5))

cells_rankings <- AUCell_buildRankings(exprMatrix)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=5, nCores=1)
selectedThresholds <- NULL
cellsTsne<- Rtsne::Rtsne(t(exprMatrix),max_iter = 10)$Y
rownames(cellsTsne) <- colnames(exprMatrix)
######


library(shiny); library(rbokeh)

# Create app 
aucellApp <- AUCell_createViewerApp(auc=cells_AUC, thresholds=selectedThresholds, 
                                    tSNE=cellsTsne, exprMat=exprMatrix)

# The exact commands to lauch the app depend on the R settings 
# For example:
# options(shiny.host="0.0.0.0") 
# savedSelections <- runApp(aucellApp) 
# (see Shiny's doc for help)

