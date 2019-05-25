#' @title orderAUC
#' @description Reorder the gene-sets based on AUC similarity 
#' @param auc AUC (as returned by calcAUC)
#' @return gene-set names in the suggested order
#' @importFrom stats as.dist hclust
# @importFrom S4Vectors as.dist
#' @importFrom stats cor
#' @examples 
#' # cellsAUC <- cellsAUC[orderAUC(cellsAUC),]
#' @export
orderAUC <- function(auc){
  variableGs <- names(which(apply(getAUC(auc), 1, sd) > 0))
  gsDist <-as.dist(1-cor(t(getAUC(auc)[variableGs,]), method="spear"))
  gsClust <- hclust(gsDist, method="ward.D2")
  gsClusters <- setNames(dynamicTreeCut::cutreeDynamic(gsClust, distM=as.matrix(gsDist), verbose = FALSE), gsClust$labels)
  gsOrder <- gsClust$labels[gsClust$order]
  gsOrder <- gsOrder[order(gsClusters[gsOrder], decreasing = TRUE)]
  return(gsOrder)
}
