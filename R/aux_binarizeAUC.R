
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- getThresholdSelected(thresholds)
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  ### Convert to matrix (regulons with zero assigned cells are lost)
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}

