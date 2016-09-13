# Runs .auc.assignmnetThreshold on each gene set (AUC)

# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
# @import data.table
#'
#' @title to do
#' @description to do
#' @param to do
#' @details # required?
#' @return to do
#' @examples  #to do
#' @export
AUCell.exploreThresholds <- function(aucMatrix, seed=987, thrP=0.01, nCores=1, plotHist=TRUE, densAdjust=2, assignCells=FALSE, verbose=TRUE) # nDist=NULL,
{
  suppressPackageStartupMessages(library(data.table)) # NEW
  if(nCores > 1 && plotHist) {
    nCores <- 1
    warning("Cannot plot in paralell")
  }

  if(nCores==1)
  {
    assigment <- lapply(colnames(aucMatrix), function(gSetName)
    {
      auc <- aucMatrix[,gSetName]
      aucThr <- .auc.assignmnetThreshold(auc, thrP=thrP, seed=seed, plotHist=plotHist, gSetName=gSetName, densAdjust=densAdjust)
      assignedCells <- NULL
      if(!is.null(aucThr)) assignedCells <- names(auc)[which(auc>aucThr$selected)]

      list(aucThr=aucThr, seed=seed, assignment=assignedCells)

    })
    names(assigment) <- colnames(aucMatrix)

  }else
  {
    # Run each geneSet in parallel
    suppressMessages(require("doMC", quietly=TRUE)); require("doRNG", quietly=TRUE)
    registerDoMC(nCores)
    if(verbose) message(paste("Using", getDoParWorkers(), "cores."))

    assigment <- foreach(gSetName=colnames(aucMatrix)) %dorng%
    {
      auc <- aucMatrix[,gSetName]
      aucThr <- .auc.assignmnetThreshold(auc, thrP=thrP, seed=seed, plotHist=plotHist, gSetName=gSetName)
      assignedCells <- NULL
      if(!is.null(aucThr)) assignedCells <- names(auc)[which(auc>aucThr$selected)]

      tmp <- list(aucThr=aucThr, seed=seed, assignment=assignedCells)
      return(setNames(list(tmp), gSetName))
    }
    attr(assigment, "rng") <- NULL
    assigment <- unlist(assigment, recursive = FALSE)[colnames(aucMatrix)]
  }
  # If cell assignment is not requested, remove from the output (it was initially calculated to know the number of cells)
  if(!assignCells) assigment <- lapply(assigment, function(x) x[which(!names(x) %in% "assignment")])

  assigment
}



#
# library(mixtools)
# # pdf("auc_thresholds_Zeisel.pdf")
# # par(mfrow=c(2,1))
# # cellsAssigned <- list()

#   cellAssignment <- unname(sapply(cellAssignment[,3], function(x) sapply(strsplit(gsub("]","",x), "_[", fixed=TRUE), function(y) c(y[[2]],y[[1]]))))
#   cellsAssigned[[sig]]  <- setNames(cellAssignment[1,], cellAssignment[2,])
#
#
# cellTypeAssigned_mixAUC <- c(lapply(names(cellsAssigned), function(sig) setNames(paste(rep(sig, length(cellsAssigned[[sig]]))), names(cellsAssigned[[sig]]))))
# names(cellTypeAssigned_mixAUC) <- names(cellsAssigned)
# cellTypeAssigned_mixAUC <- sapply(cellTypeAssigned_mixAUC, names)
# # length(unlist(cellTypeAssigned_mixAUC))
# # # 7690
# save(cellTypeAssigned_mixAUC, file="cellTypeAssigned_literatureSigs_v3(AUC threholds).RData")
#
