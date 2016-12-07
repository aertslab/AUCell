# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @import data.table
#'
#' @title Build gene expression rankings for each cell
#' @description Builds the "rankings" for each cell: expression-based ranking for all the genes in each cell.
#'
#' The genes with same expression value are shuffled. Therefore, genes with expression '0' are randomly sorted at the end of the ranking.
#'
#' These "rankings" can be seen as a new representation of the original dataset. Once they are calculated, they can be saved for future analyses.
#' @param exprMat Expression matrix (genes as rows, cells as columns)
#' @param plotStats Should the function plot the expression boxplots/histograms? (TRUE / FALSE). These plots can also be produced with the function \code{\link{plotGeneCount}}.
#' @param nCores Number of cores to use for computation.
#' @param verbose Should the function show progress messages? (TRUE / FALSE)
#' @return data.table of genes (row) by cells (columns) with the ranking of the gene within the cell.
#' @details
#' It is important to check that most cells have at least the number of expressed/detected genes that are going to be used to calculate the AUC (`aucMaxRank` in `calcAUC()`). The histogram provided by `AUCell.buildRankings()` allows to quickly check this distribution. `plotGeneCount(exprMatrix)` allows to obtain only the plot before building the rankings.
#' @seealso Next step in the workflow: \code{\link{AUCell.calcAUC}}.
#'
#' See the package vignette for examples and more details: \code{vignette("AUCell")}
#' @example inst/examples/example_AUCell.buildRankings.R
#' @export
AUCell.buildRankings <- function(exprMat, plotStats=TRUE, nCores=1, verbose=TRUE)
{
  if(!is.data.table(exprMat)) exprMat <- data.table(exprMat, keep.rownames=TRUE)
  setkey(exprMat, "rn") # (reorders rows)

  if(plotStats)
  {
    msg <- tryCatch(plotGeneCount(exprMat[,-"rn", with=FALSE], verbose=verbose),
                                                    error = function(e) {
                                                      return(e)
                                                    })
    if("error" %in% class(msg)) warning(paste("There has been an error in plotGeneCount() [Message: ", msg$message, "]. Proceeding to calculate the rankings...", sep=""))
  }

  colsNam <- colnames(exprMat)[-1] # 1=row names
  if(nCores==1)
  {
    # The rankings are saved in exprMat (i.e. By reference)
    exprMat[, (colsNam):=lapply(-.SD, frank, ties.method="random"), .SDcols=colsNam]

  }else
  {
    suppressMessages(require("doMC", quietly=TRUE)); require("doRNG", quietly=TRUE)
    registerDoMC(nCores)
    if(verbose) message(paste("Using", getDoParWorkers(), "cores."))

    suppressWarnings(colsNamsGroups <- split(colsNam, (1:length(colsNam)) %% nCores)) # Expected warning: Not multiple
    rowNams <- exprMat$rn
    exprMat <- foreach(colsGr=colsNamsGroups, .combine="cbind") %dorng%
    {
      # exprMat[, (colsGr):=lapply(-.SD, frank, ties.method="random"), .SDcols=colsGr]  # Edits by reference: how to make it work in paralell...?
      subMat <- exprMat[,colsGr, with=FALSE]
      subMat[, (colsGr):=lapply(-.SD, frank, ties.method="random"), .SDcols=colsGr]
    }
    exprMat <- data.table(rn=rowNams, exprMat[,colsNam, with=FALSE]) # Keep initial order & recover rownames
    setkey(exprMat, "rn")
  }
  return(exprMat)
}

