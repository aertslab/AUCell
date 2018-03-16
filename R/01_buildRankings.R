# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @import data.table
#' @import SummarizedExperiment
#' @importFrom methods new
#'
#' @title Build gene expression rankings for each cell
#' @description Builds the "rankings" for each cell:
#' expression-based ranking for all the genes in each cell.
#'
#' The genes with same expression value are shuffled.
#' Therefore, genes with expression '0' are randomly sorted at the
#' end of the ranking.
#'
#' These "rankings" can be seen as a new representation of the original dataset.
#' Once they are calculated, they can be saved for future analyses.
#' @param exprMat Expression matrix (genes as rows, cells as columns)
#' The expression matrix can also be provided as one of the Bioconductor classes:
#' \itemize{
#' \item \link[Biobase]{ExpressionSet}:
#' The matrix will be obtained through exprs(exprMatrix)
#' \item \link[SummarizedExperiment]{RangedSummarizedExperiment} and derived classes (e.g. \link[SingleCellExperiment]{SingleCellExperiment} ):
#' The matrix will be obtained through assay(exprMatrix),
#' -which will extract the first assay (usually the counts)-
#' or the assay name given in 'assayName'
#' }
#' @param plotStats Should the function plot the expression boxplots/histograms?
#' (TRUE / FALSE). These plots can also be produced
#' with the function \code{\link{plotGeneCount}}.
#' @param nCores Number of cores to use for computation.
#' @param verbose Should the function show progress messages? (TRUE / FALSE)
#' @param assayName Name of the assay containing the expression matrix (e.g. in \link[SingleCellExperiment]{SingleCellExperiment} objects)
#' @return data.table of genes (row) by cells (columns)
#' with the ranking of the gene within the cell.
#' @details
#' It is important to check that most cells have at least the number of
#' expressed/detected genes that are going to be used to calculate the AUC
#' (`aucMaxRank` in `calcAUC()`).
#' The histogram provided by `AUCell_buildRankings()`
#' allows to quickly check this distribution.
#' `plotGeneCount(exprMatrix)` allows to obtain only the plot before
#' building the rankings.
#' @seealso Next step in the workflow: \code{\link{AUCell_calcAUC}}.
#'
#' See the package vignette for examples and more details:
#' \code{vignette("AUCell")}
#' @example inst/examples/example_AUCell_buildRankings.R
#' @rdname AUCell_buildRankings
#' @export
setGeneric("AUCell_buildRankings", signature="exprMat",
  function(exprMat, plotStats=TRUE, nCores=1, verbose=TRUE)
  {
   standardGeneric("AUCell_buildRankings")
  })

#' @rdname AUCell_buildRankings
#' @aliases AUCell_buildRankings,matrix-method
setMethod("AUCell_buildRankings", "matrix",
  function(exprMat, plotStats=TRUE, nCores=1, verbose=TRUE)
  {
    .AUCell_buildRankings(exprMat=exprMat, plotStats=plotStats, nCores=nCores, verbose=verbose)
  })

#' @rdname AUCell_buildRankings
#' @aliases AUCell_buildRankings,SummarizedExperiment-method
setMethod("AUCell_buildRankings", "SummarizedExperiment",
function(exprMat, plotStats=TRUE, nCores=1, verbose=TRUE, assayName=NULL)
  {
    if(is.null(assayName))
    {
      if(length(SummarizedExperiment::assays(exprMat))>1)
        warning("More than 1 assays are available. Using the first one (", names(assays(exprMat))[1], ").")
      
      exprMat <- SummarizedExperiment::assay(exprMat)
    }else{
      if(!assayName %in% names(assays(exprMat)))
        stop("There isn't an assay named ", assayName, ".")
      
      exprMat <- SummarizedExperiment::assays(exprMat)[[assayName]]
    }
      
    .AUCell_buildRankings(exprMat=exprMat, plotStats=plotStats, nCores=nCores, verbose=verbose)
  })

#' @rdname AUCell_buildRankings
#' @aliases AUCell_buildRankings,ExpressionSet-method
setMethod("AUCell_buildRankings", "ExpressionSet",
  function(exprMat, plotStats=TRUE, nCores=1, verbose=TRUE)
  {
    exprMat <- Biobase::exprs(exprMat)
    .AUCell_buildRankings(exprMat=exprMat, plotStats=plotStats, nCores=nCores, verbose=verbose)
  })

.AUCell_buildRankings <- function(exprMat, plotStats=TRUE, nCores=1, verbose=TRUE)
{
  if(!data.table::is.data.table(exprMat))
    exprMat <- data.table::data.table(exprMat, keep.rownames=TRUE)
    # TO DO: Replace by sparse matrix??? (e.g. dgTMatrix)
  data.table::setkey(exprMat, "rn") # (reorders rows)

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
    exprMat[, (colsNam):=lapply(-.SD, data.table::frank, ties.method="random"),
            .SDcols=colsNam]

  }else
  {
    # doRNG::registerDoRNG(nCores)
    doParallel::registerDoParallel()
    options(cores=nCores)

    if(verbose)
      message("Using ", foreach::getDoParWorkers(), " cores.")

    # Expected warning: Not multiple
    suppressWarnings(colsNamsGroups <- split(colsNam,
                                             (seq_along(colsNam)) %% nCores))
    rowNams <- exprMat$rn

    colsGr <- NULL
    "%dopar%"<- foreach::"%dopar%"
    suppressPackageStartupMessages(exprMat <-
                      doRNG::"%dorng%"(foreach::foreach(colsGr=colsNamsGroups,
                                                        .combine=cbind),
    {
      # Edits by reference: how to make it work in paralell...?
      subMat <- exprMat[,colsGr, with=FALSE]
      subMat[, (colsGr):=lapply(-.SD, data.table::frank, ties.method="random"),
             .SDcols=colsGr]
    }))
    # Keep initial order & recover rownames
    exprMat <- data.table::data.table(rn=rowNams, exprMat[,colsNam, with=FALSE])
    data.table::setkey(exprMat, "rn")
  }

  rn <- exprMat$rn
  exprMat <- as.matrix(exprMat[,-1])
  rownames(exprMat) <- rn

  # return(matrixWrapper(matrix=exprMat, rowType="gene", colType="cell",
  #                      matrixType="Ranking", nGenesDetected=nGenesDetected))
  names(dimnames(exprMat)) <- c("genes", "cells")
  new("aucellResults",
      SummarizedExperiment::SummarizedExperiment(assays=list(ranking=exprMat)),
      nGenesDetected=nGenesDetected)
}

