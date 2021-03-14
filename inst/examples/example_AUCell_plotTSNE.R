

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

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=5, nCores=1)
selectedThresholds <- rowMeans(getAUC(cells_AUC))
cellsTsne<- Rtsne::Rtsne(t(exprMatrix),max_iter = 10)$Y
# cellsTsne<- tsne::tsne(t(exprMatrix),max_iter = 10)
rownames(cellsTsne) <- colnames(exprMatrix)
######


par(mfrow=c(2,3))
thrs <- AUCell_plotTSNE(tSNE=cellsTsne, exprMat=NULL,
                        cellsAUC=cells_AUC, thresholds=selectedThresholds, 
                        plots = c("histogram", "binaryAUC", "AUC"))
 


#####
# Color based on the known phenodata:
cellInfo <- data.frame(cellType1=sample(LETTERS[1:3],ncol(exprMatrix), replace=TRUE), 
                       cellType2=sample(letters[5:7],ncol(exprMatrix), replace=TRUE), 
                       nGenes=abs(rnorm(ncol(exprMatrix))), 
                       row.names=colnames(exprMatrix))
colVars <- list(cellType2=setNames(c("skyblue","magenta", "darkorange"),letters[5:7]))
# dev.off()
plotTsne_cellProps(cellsTsne, cellInfo, colVars=colVars)

