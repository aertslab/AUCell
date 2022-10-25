# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @rawNamespace import(data.table, except = shift)
# @import SummarizedExperiment
#' @importFrom methods new
#' @importFrom DelayedMatrixStats colRanks
#' @importClassesFrom DelayedArray DelayedArray
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
#' The expression matrix can also be provided as one of the R/Bioconductor classes:
#' \itemize{
#' \item \link[Matrix]{dgCMatrix-class}:
#' Sparse matrix 
#' \item \link{RangedSummarizedExperiment} and derived classes (e.g. \link[SingleCellExperiment]{SingleCellExperiment} ):
#' The matrix will be obtained through assay(exprMatrix),
#' -which will extract the first assay (usually the counts)-
#' or the assay name given in 'assayName'
#' \item \code{ExpressionSet}:
#' The matrix will be obtained through exprs(exprMatrix)
#' }
#' @param featureType Name for the rows (e.g. "genes"). Only for naming the rankings, not used internally.
#' @param plotStats Should the function plot the expression boxplots/histograms?
#' (TRUE / FALSE). These plots can also be produced
#' with the function \code{\link{plotGeneCount}}.
#' @param splitByBlocks Whether to split the matrix by blocks in the ranking calculation. 
#' Allows using multiple cores.
#' FALSE by default. If using sparse matrices it is automatically set to TRUE.
#' @param BPPARAM Set to use multiple cores. Only used if 'splitByBlocks=TRUE'
#' @param verbose Should the function show progress messages? (TRUE / FALSE)
#' @param assayName Name of the assay containing the expression matrix (e.g. in \link[SingleCellExperiment]{SingleCellExperiment} objects)
#' @param keepZeroesAsNA Convert zeroes to NA instead of locating randomly at the end of the ranking.
#' @param ... Other arguments
#' @param nCores Deprecated
#' @param mctype Deprecated
#' @return Ranking of the feature within the cell (features as rows, cells as columns)
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
  function(exprMat, featureType="genes", plotStats=TRUE, splitByBlocks=FALSE,
           BPPARAM=NULL,
           keepZeroesAsNA=FALSE, verbose=TRUE, 
           nCores=NULL, mctype=NULL, ...)
  {
   standardGeneric("AUCell_buildRankings")
  })

# Sparse matrix
#' @rdname AUCell_buildRankings
#' @aliases AUCell_buildRankings,dgCMatrix-method
setMethod("AUCell_buildRankings", "dgCMatrix",
          function(exprMat, featureType="genes", plotStats=TRUE, splitByBlocks=TRUE,
                   BPPARAM=NULL,
                   keepZeroesAsNA=FALSE, verbose=TRUE, nCores=NULL, mctype=NULL)
          {
            .AUCell_buildRankings(exprMat=exprMat, featureType=featureType,
                                  splitByBlocks=TRUE, BPPARAM=BPPARAM,
                                  plotStats=plotStats, keepZeroesAsNA=keepZeroesAsNA, verbose=verbose,
                                  nCores=nCores, mctype=mctype)
          })

# Non-sparse matrix
#' @rdname AUCell_buildRankings
#' @aliases AUCell_buildRankings,matrix-method
setMethod("AUCell_buildRankings", "matrix",
  function(exprMat, featureType="genes", plotStats=TRUE, splitByBlocks=FALSE,
           BPPARAM=NULL,
           keepZeroesAsNA=FALSE, verbose=TRUE, 
           nCores=NULL, mctype=NULL)
  {
    .AUCell_buildRankings(exprMat=exprMat, featureType=featureType, 
                          splitByBlocks=splitByBlocks, BPPARAM=BPPARAM,
                          plotStats=plotStats, keepZeroesAsNA=keepZeroesAsNA, verbose=verbose,
                          nCores=nCores, mctype=mctype)
  })

#' @rdname AUCell_buildRankings
#' @aliases AUCell_buildRankings,ExpressionSet-method
setMethod("AUCell_buildRankings", "ExpressionSet",
          function(exprMat, featureType="genes", plotStats=TRUE, splitByBlocks=FALSE, 
                   BPPARAM=NULL,
                   keepZeroesAsNA=FALSE, verbose=TRUE, 
                   nCores=NULL, mctype=NULL)
          {
            exprMat <- Biobase::exprs(exprMat)
            .AUCell_buildRankings(exprMat=exprMat, featureType=featureType,
                                  splitByBlocks=splitByBlocks, BPPARAM=BPPARAM,
                                  plotStats=plotStats, keepZeroesAsNA=keepZeroesAsNA, verbose=verbose,
                                  nCores=nCores, mctype=mctype)
          })

#' @rdname AUCell_buildRankings
#' @aliases AUCell_buildRankings,SummarizedExperiment-method
setMethod("AUCell_buildRankings", "SummarizedExperiment",
function(exprMat, featureType="genes", plotStats=TRUE, splitByBlocks=FALSE,
         BPPARAM=NULL,
         keepZeroesAsNA=FALSE, verbose=TRUE, assayName=NULL, 
         nCores=NULL, mctype=NULL)
  {
    if(is.null(assayName))
    {
      if(length(SummarizedExperiment::assays(exprMat))>1)
        warning("More than 1 assays are available. Using the first one (", 
                names(assays(exprMat))[1], ").")
      
      exprMat <- SummarizedExperiment::assay(exprMat)
    }else{
      if(!assayName %in% names(assays(exprMat)))
        stop("There isn't an assay named ", assayName, ".")
      
      exprMat <- SummarizedExperiment::assays(exprMat)[[assayName]]
    }
    .AUCell_buildRankings(exprMat=exprMat, featureType=featureType, 
                          splitByBlocks=splitByBlocks, BPPARAM=BPPARAM,
                          plotStats=plotStats, keepZeroesAsNA=keepZeroesAsNA, verbose=verbose,
                          nCores=nCores, mctype=mctype)
  })


### Run without delayed array, it can be added outside if needed
# Accepts regular or sparse matrix ()
.AUCell_buildRankings <- function(exprMat, featureType="genes",
                                  splitByBlocks=FALSE,
                                  keepZeroesAsNA=FALSE, 
                                  BPPARAM=NULL,
                                  plotStats=TRUE, verbose=TRUE,
                                  nCores=NULL, mctype=NULL)
{
  if((!splitByBlocks) && ("dgCMatrix" %in% class(exprMat)))
    splitByBlocks <- TRUE
  
  if(!is.null(nCores))
      warning("nCores is no longer used. It will be deprecated in the next AUCell version.")
    
  ## Still needed? 
  if(keepZeroesAsNA){
    zeroesCoords <- which(exprMat==0, arr.ind=TRUE)
  }
  
  ## Calculate expression stats?
  nGenesDetected <- numeric(0)
  if(plotStats)
  {
    msg <- tryCatch(plotGeneCount(exprMat, plotStats=plotStats, verbose=verbose),
                    error = function(e) {
                      return(e)
                    })
    if(methods::is(msg,"error")) {
      warning("There has been an error in plotGeneCount() [Message: ",
              msg$message, "]. Proceeding to calculate the rankings...", sep="")
    }else{
      if(is.numeric(nGenesDetected))
        nGenesDetected <- msg
    }
  }
  
  ## Rank genes
  rowNames <- rownames(exprMat)
  colNames <- colnames(exprMat)
  exprMat <- -exprMat # ro rank decreasingly
  
  # Would prefer not to use the block inside the function... but
  # the method for sparse matrices does not allow ties.method='random'
  if(splitByBlocks){
    exprMat <- do.call(cbind,
                       blockApply(DelayedArray(exprMat),
                                  FUN=DelayedMatrixStats::colRanks,
                                  ties.method="random", preserveShape=TRUE,
                                  BPPARAM=BPPARAM,
                                  grid=colAutoGrid(exprMat))) 
  }else{
    exprMat <- DelayedMatrixStats::colRanks(exprMat, ties.method="random", preserveShape=TRUE, decreasing=TRUE) 
  }
  
  rownames(exprMat) <- rowNames
  colnames(exprMat) <- colNames
  
  if(keepZeroesAsNA){
    exprMat[which(zeroesCoords==0, arr.ind=TRUE)] <- NA
  }
  
  ## Format & return
  names(dimnames(exprMat)) <- c(featureType, "cells")
  new("aucellResults",
      SummarizedExperiment::SummarizedExperiment(assays=list(ranking=exprMat)),
      nGenesDetected=nGenesDetected)
}

