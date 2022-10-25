# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
#
#' @title AUCell_assignCells
#' @description Assigns whether the gene sets are considered "active" on each cell
#' based on the given thresholds
#' @param cellsAUC AUC object returned by \code{\link{AUCell_calcAUC}}.
#' @param thresholds Thresholds selected for each geneset (named vector).
#' @param nCores Number of cores to use for computation.
#' @return List with the following elements for each gene-set:
#' \itemize{
#' \item 'aucThr' threshold value, in the same format as AUCell_exploreThresholds()
#' \item 'assignment' List of cells that pass the selected AUC threshold
#' }
#' @seealso Previous step in the workflow: \code{\link{AUCell_calcAUC}} 
#' and optionally \code{\link{AUCell_exploreThresholds}} 
#'
#' See the package vignette for examples and more details:
#' \code{vignette("AUCell")}
#' @example inst/examples/example_AUCell_exploreThresholds.R
#' 
# Function copied from AUCell_exploreThresholds and removed the assingment
# (not the best practice: in the future better use common code)
#' @export
AUCell_assignCells <- function(cellsAUC, thresholds, nCores=1)
{
  #####  Extract AUC matrix
  aucMatrix <- NULL
  if(methods::is(cellsAUC,"aucellResults")){
    aucMatrix <- getAUC(cellsAUC)
  }
  
  if(methods::is(cellsAUC,"matrixWrapper"))
  {
    stop('The AUC was calculated with a previous AUCell version. \n",
         "Please update them with updateAucellResults(..., objectType="AUC")')
  }
  
  if(is.matrix(cellsAUC)) {
    aucMatrix <- cellsAUC
  }
  if(!is.matrix(aucMatrix)) stop("cellsAUC should contain the AUC values.")
  
  rowSumAUC <- rowSums(aucMatrix)
  if(any(rowSumAUC==0)) warning("Skipping genesets with all AUC 0: ", 
                                paste(names(rowSumAUC)[which(rowSumAUC==0)], collapse=", "), 
                                immediate. = TRUE)
  aucMatrix <- aucMatrix[rowSumAUC>0,,drop=FALSE]
  #########
  
  if(nCores==1)
  {
    gSetName <- NULL
    assigment <- lapply(rownames(aucMatrix), function(gSetName)
    {
      aucThr <- thresholds[gSetName]
      assignedCells <- NULL
      if(!is.null(aucThr)){
        auc <- aucMatrix[gSetName,]
        assignedCells <- names(auc)[which(auc>aucThr)]
      }

      list(aucThr=list(selected=aucThr), assignment=assignedCells)

    })
    names(assigment) <- rownames(aucMatrix)
    
  }else
  {
    # Run each geneSet in parallel
    doRNG::registerDoRNG(nCores)
    if(verbose)
      message("Using ", foreach::getDoParWorkers(), " cores.")

    "%dopar%"<- foreach::"%dopar%"
    suppressPackageStartupMessages(assigment <-
                                     doRNG::"%dorng%"(foreach::foreach(gSetName=rownames(aucMatrix)),
                                                      {
                                                        aucThr <- thresholds[gSetName]
                                                        assignedCells <- NULL
                                                        if(!is.null(aucThr)){
                                                          auc <- aucMatrix[gSetName,]
                                                          assignedCells <- names(auc)[which(auc>aucThr$selected)]
                                                        }
                                                        tmp <- list(aucThr=list(selected=aucThr), assignment=assignedCells)
                                                        return(setNames(list(tmp), gSetName))
                                                      }))
    attr(assigment, "rng") <- NULL
    assigment <- unlist(assigment, recursive = FALSE)[rownames(aucMatrix)]
  }

  return(assigment)
}





