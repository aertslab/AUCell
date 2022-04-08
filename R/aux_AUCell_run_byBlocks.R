#' @importFrom DelayedArray DelayedArray blockApply colAutoGrid
#' @title Run AUCell 
#' @description Runs AUCell (calculates the ranking + score genesets)
#' @param exprMat Expression matrix (genes/regions as rows, cells as columns)
#' The expression matrix can also be provided as one of the R/Bioconductor classes:
#' \itemize{
#' \item \link[Matrix]{dgCMatrix-class}:
#' Sparse matrix 
#' \item \link{RangedSummarizedExperiment} and derived classes (e.g. \link[SingleCellExperiment]{SingleCellExperiment} ):
#' The matrix will be obtained through assay(exprMatrix),
#' -which will extract the first assay (usually the counts)-
#' or the assay name given in 'assayName'
#' }
#' @param geneSets List of gene-sets (or signatures) to test in the cells.
#' The gene-sets should be provided as \code{\link{GeneSet}},
#' \code{\link{GeneSetCollection}} or character list (see examples).
#' @param featureType Name for the rows (e.g. "genes"). Only for naming the rankings, not used internally.
#' @param keepZeroesAsNA Convert zeroes to NA instead of locating randomly at the end of the ranking.
#' @param normAUC Whether to normalize the maximum possible AUC to 1 (Default: TRUE).
#' @param aucMaxRank Threshold to calculate the AUC (see 'details' section)
#' @param BPPARAM Set to use multiple cores. Only used if 'splitByBlocks=TRUE'
#' @param assayName Name of the assay containing the expression matrix (e.g. in \link[SingleCellExperiment]{SingleCellExperiment} objects)
#' @param ... Other arguments
#' @return Matrix with the AUC values (gene-sets as rows, cells as columns).
#' @details In a simplified way, the AUC value represents the fraction of genes,
#'  within the top X genes in the ranking, that are included in the signature.
#' The parameter 'aucMaxRank' allows to modify the number of genes
#' (maximum ranking) that is used to perform this computation.
#' By default, it is set to 5\% of the total number of genes in the rankings.
#' Common values may range from 1 to 20\%.
#' @seealso Includes \code{\link{AUCell_buildRankings}} and \code{\link{AUCell_calcAUC}}.
#' Next step in the workflow: \code{\link{AUCell_exploreThresholds}}.
#'
#' See the package vignette for examples and more details:
#' \code{vignette("AUCell")}
#' @example inst/examples/example_AUCell_run.R
#'
#' @rdname AUCell_run
#' @export
setGeneric("AUCell_run", signature="exprMat",
           function(exprMat, 
                    geneSets, 
                    featureType='genes',
                    keepZeroesAsNA=FALSE, 
                    normAUC=TRUE, 
                    aucMaxRank=ceiling(0.05*nrow(exprMat)), 
                    BPPARAM=NULL,
                    ...)
           {
             standardGeneric("AUCell_run")
           })

# Sparse matrix
#' @rdname AUCell_run
#' @aliases AUCell_run,dgCMatrix-method
setMethod("AUCell_run", "dgCMatrix", 
          function(exprMat, 
                   geneSets, 
                   featureType='genes',
                   keepZeroesAsNA=FALSE, 
                   normAUC=TRUE, 
                   aucMaxRank=ceiling(0.05*nrow(exprMat)), 
                   BPPARAM=NULL)
          {
                        .AUCell_run(exprMat=exprMat, geneSets=geneSets, featureType=featureType, 
                                    keepZeroesAsNA=keepZeroesAsNA, 
                                    normAUC=normAUC, aucMaxRank=aucMaxRank, 
                                    BPPARAM=BPPARAM)
          })

# Non-sparse matrix
#' @rdname AUCell_run
#' @aliases AUCell_run,matrix-method
setMethod("AUCell_run", "matrix",
          function(exprMat, 
                   geneSets, 
                   featureType="genes", 
                   keepZeroesAsNA=FALSE, 
                   normAUC=TRUE, 
                   aucMaxRank=ceiling(0.05*nrow(exprMat)), 
                   BPPARAM=NULL)
          {
            .AUCell_run(exprMat=exprMat, geneSets=geneSets, featureType=featureType, 
                        keepZeroesAsNA=keepZeroesAsNA, 
                        normAUC=normAUC, aucMaxRank=aucMaxRank, 
                        BPPARAM=BPPARAM)
          })

#' @rdname AUCell_run
#' @aliases AUCell_run,SummarizedExperiment-method
setMethod("AUCell_run", "SummarizedExperiment",
          function(exprMat, 
                   geneSets, 
                   featureType="genes", 
                   keepZeroesAsNA=FALSE, 
                   normAUC=TRUE, 
                   aucMaxRank=ceiling(0.05*nrow(exprMat)), 
                   BPPARAM=NULL, assayName=NULL)
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
            
            .AUCell_run(exprMat=exprMat, geneSets=geneSets, featureType=featureType, 
                        keepZeroesAsNA=keepZeroesAsNA, 
                        normAUC=normAUC, aucMaxRank=aucMaxRank, 
                        BPPARAM=BPPARAM)
          })

.AUCell_run <- function(exprMat, 
                      geneSets, 
                      featureType="genes", 
                      keepZeroesAsNA=FALSE, 
                      normAUC=TRUE, 
                      aucMaxRank=ceiling(0.05*nrow(exprMat)), 
                      BPPARAM=NULL) 
{
  collected <- blockApply(DelayedArray(exprMat), 
                          FUN=.AUCell_run_internal, 
                          geneSets=geneSets,
                          featureType=featureType,
                          keepZeroesAsNA=keepZeroesAsNA,
                          normAUC=normAUC,
                          aucMaxRank=aucMaxRank,
                          BPPARAM=BPPARAM,
                          grid=colAutoGrid(exprMat)
  )
  do.call(cbind, collected)
}

.AUCell_run_internal <- function(exprMat_block, geneSets, featureType, keepZeroesAsNA, normAUC, aucMaxRank) 
{
  ranked <- AUCell_buildRankings(exprMat_block, splitByBlocks=FALSE, 
                                 featureType=featureType,
                                 keepZeroesAsNA=keepZeroesAsNA, 
                                 plotStats=FALSE, verbose=FALSE)
  AUCell_calcAUC(geneSets, ranked, normAUC=normAUC, aucMaxRank=aucMaxRank, verbose=FALSE)
}
