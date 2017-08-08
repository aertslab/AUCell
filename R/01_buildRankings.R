# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @import data.table
#' @import SummarizedExperiment
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
#' It is important to check that most cells have at least the number of expressed/detected genes that are going to be used to calculate the AUC (`aucMaxRank` in `calcAUC()`). The histogram provided by `AUCell_buildRankings()` allows to quickly check this distribution. `plotGeneCount(exprMatrix)` allows to obtain only the plot before building the rankings.
#' @seealso Next step in the workflow: \code{\link{AUCell_calcAUC}}.
#'
#' See the package vignette for examples and more details: \code{vignette("AUCell")}
#' @example inst/examples/example_AUCell_buildRankings.R
#' @export
AUCell_buildRankings <- function(exprMat, plotStats=TRUE, nCores=1, verbose=TRUE)
{
  if(!is.data.table(exprMat))
    exprMat <- data.table(exprMat, keep.rownames=TRUE)  # TO DO: Replace by sparse matrix??? (e.g. dgTMatrix)
  setkey(exprMat, "rn") # (reorders rows)

  nGenesDetected <- numeric(0)
  if(plotStats)
  {
    msg <- tryCatch(plotGeneCount(exprMat[,-"rn", with=FALSE], verbose=verbose),
            error = function(e) {
              return(e)
            })
    if("error" %in% class(msg)) {
      warning("There has been an error in plotGeneCount() [Message: ",
              msg$message, "]. Proceeding to calculate the rankings...", sep="")
    }else{
      if(is.numeric(nGenesDetected))
        nGenesDetected <- msg
    }

  }

  colsNam <- colnames(exprMat)[-1] # 1=row names
  if(nCores==1)
  {
    # The rankings are saved in exprMat (i.e. By reference)
    exprMat[, (colsNam):=lapply(-.SD, frank, ties.method="random"), .SDcols=colsNam]

  }else
  {
    # doRNG::registerDoRNG(nCores)
    doParallel::registerDoParallel()
    options(cores=nCores)

    if(verbose)
      message("Using ", foreach::getDoParWorkers(), " cores.")

    suppressWarnings(colsNamsGroups <- split(colsNam, (seq_along(colsNam)) %% nCores)) # Expected warning: Not multiple
    rowNams <- exprMat$rn

    "%dopar%"<- foreach::"%dopar%"
    suppressPackageStartupMessages(exprMat <-  doRNG::"%dorng%"(foreach::foreach(colsGr=colsNamsGroups, .combine=cbind),
    {
      # exprMat[, (colsGr):=lapply(-.SD, frank, ties.method="random"), .SDcols=colsGr]
      # Edits by reference: how to make it work in paralell...?
      subMat <- exprMat[,colsGr, with=FALSE]
      subMat[, (colsGr):=lapply(-.SD, frank, ties.method="random"), .SDcols=colsGr]
    }))
    # Keep initial order & recover rownames
    exprMat <- data.table(rn=rowNams, exprMat[,colsNam, with=FALSE])
    setkey(exprMat, "rn")
  }

  rn <- exprMat$rn
  exprMat <- as.matrix(exprMat[,-1])
  rownames(exprMat) <- rn

  # return(matrixWrapper(matrix=exprMat, rowType="gene", colType="cell",
  #                      matrixType="Ranking", nGenesDetected=nGenesDetected))
  names(dimnames(exprMat)) <- c("genes", "cells")
  new("aucellResults",
      SummarizedExperiment(assays=list(ranking=exprMat)),
      nGenesDetected=nGenesDetected)
}

