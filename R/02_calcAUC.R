# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
#' @rawNamespace import(data.table, except = shift)
#' @import GSEABase
#' @importFrom stats setNames
#' @importFrom methods new
#' @title Calculate AUC
#' @description Calculates the 'AUC' for each gene-set in each cell.
#' @param geneSets List of gene-sets (or signatures) to test in the cells.
#' The gene-sets should be provided as \code{\link{GeneSet}},
#' \code{\link{GeneSetCollection}} or character list (see examples).
#' @param rankings 'Rankings' created for this dataset with
#' \code{\link{AUCell_buildRankings}}.
#' @param nCores Number of cores to use for computation.
#' @param normAUC Wether to normalize the maximum possible AUC to 1 (Default: TRUE).
#' @param aucMaxRank Threshold to calculate the AUC (see 'details' section)
#' @param verbose Should the function show progress messages? (TRUE / FALSE)
#' @return Matrix with the AUC values (gene-sets as rows, cells as columns).
#' @details In a simplified way, the AUC value represents the fraction of genes,
#'  within the top X genes in the ranking, that are included in the signature.
#' The parameter 'aucMaxRank' allows to modify the number of genes
#' (maximum ranking) that is used to perform this computation.
#' By default, it is set to 5\% of the total number of genes in the rankings.
#' Common values may range from 1 to 20\%.
#' @seealso Previous step in the workflow: \code{\link{AUCell_buildRankings}}.
#' Next step in the workflow: \code{\link{AUCell_exploreThresholds}}.
#'
#' See the package vignette for examples and more details:
#' \code{vignette("AUCell")}
#' @example inst/examples/example_AUCell_calcAUC.R
#'
# @docType methods
#' @rdname AUCell_calcAUC
#' @export
setGeneric("AUCell_calcAUC", signature="geneSets",
  function(geneSets, rankings, nCores=1, normAUC=TRUE,
          aucMaxRank=ceiling(0.05*nrow(rankings)), verbose=TRUE)
  {
   standardGeneric("AUCell_calcAUC")
})

#' @rdname AUCell_calcAUC
#' @aliases AUCell_calcAUC,list-method
setMethod("AUCell_calcAUC", "list",
  function(geneSets, rankings, nCores=1, normAUC=TRUE,
           aucMaxRank=ceiling(0.05*nrow(rankings)), verbose=TRUE)
  {
    .AUCell_calcAUC(geneSets=geneSets, rankings=rankings, nCores=nCores,
                    aucMaxRank=aucMaxRank, verbose=verbose)
  })

#' @rdname AUCell_calcAUC
#' @aliases AUCell_calcAUC,character-method
setMethod("AUCell_calcAUC", "character",
  function(geneSets, rankings, nCores=1, normAUC=TRUE,
           aucMaxRank=ceiling(0.05*nrow(rankings)), verbose=TRUE)
  {
    geneSets <- list(geneSet=geneSets)

    .AUCell_calcAUC(geneSets=geneSets, rankings=rankings, nCores=nCores,
                    aucMaxRank=aucMaxRank, verbose=verbose)
  })

#' @rdname AUCell_calcAUC
#' @aliases AUCell_calcAUC,GeneSet-method
setMethod("AUCell_calcAUC", "GeneSet",
  function(geneSets, rankings, nCores=1, normAUC=TRUE,
           aucMaxRank=ceiling(0.05*nrow(rankings)), verbose=TRUE)
  {
    geneSets <- setNames(list(GSEABase::geneIds(geneSets)),
                         GSEABase::setName(geneSets))

    .AUCell_calcAUC(geneSets=geneSets, rankings=rankings, nCores=nCores,
                    aucMaxRank=aucMaxRank, verbose=verbose)
  })

#' @rdname AUCell_calcAUC
#' @aliases AUCell_calcAUC,GeneSetCollection-method
setMethod("AUCell_calcAUC", "GeneSetCollection",
  function(geneSets, rankings, nCores=1, normAUC=TRUE,
           aucMaxRank=ceiling(0.05*nrow(rankings)), verbose=TRUE)
  {
    geneSets <- GSEABase::geneIds(geneSets)

    .AUCell_calcAUC(geneSets=geneSets, rankings=rankings, nCores=nCores,
                    aucMaxRank=aucMaxRank, verbose=verbose)
  })

# Takes named character list as input
.AUCell_calcAUC <- function(geneSets, rankings, 
                            nCores=1, mctype=c("domc","snow")[1],
                            normAUC=TRUE,
                            aucMaxRank=ceiling(0.05*nrow(rankings)), verbose=TRUE)
{
  if(!is.list(geneSets))
    stop("geneSets should be a named list.")
  if(is.null(names(geneSets)))
    stop("geneSets should be a named list.")
  if(any(lengths(geneSets)<=0))
  {
    warning("Ignoring the following empty sets: ", paste0(names(lengths(geneSets))[which(lengths(geneSets)<=0)], collapse=", "))
    geneSets <- geneSets[which(lengths(geneSets)>0)]
  }
  if(length(geneSets)<=0) stop("No geneSets provided or remaining.")
  
  if(nCores > length(geneSets))
    nCores <- length(geneSets) # No point in using more...
  if((aucMaxRank < 300) && verbose)
    warning("Using only the first ", aucMaxRank,
            " genes (aucMaxRank) to calculate the AUC.", immediate.=TRUE)
  if(aucMaxRank <=0) stop("aucMaxRank should be a positive value.")

  if(!methods::is(rankings,"aucellResults") ||
     SummarizedExperiment::assayNames(rankings)!="ranking"){
    if(methods::is(rankings,"matrixWrapper"))
    {
      stop('These rankings were built with a previous AUCell version. \n",
           "Please update them with updateAucellResults(..., objectType="ranking")')
    }
    stop("Rankings should be a the object returned by AUCell_buildRankings()")
  } else{
    rankings <- getRanking(rankings)
  }
  # if(!key(rankings) == "rn") stop("The rankings key should be 'rn'.")

  if(!normAUC) .AUC.geneSet <- .AUC.geneSet_old()
  if(normAUC) .AUC.geneSet <- .AUC.geneSet_norm
  
  ######################################################################
  #### 1. Calculate the AUC for each gene set
  if(nCores==1)
  {
    gSetName <- NULL
    aucMatrix <- sapply(names(geneSets), function(gSetName)
      .AUC.geneSet(geneSet=geneSets[[gSetName]], rankings=rankings,
                   aucMaxRank=aucMaxRank, gSetName=gSetName))
    aucMatrix <- t(aucMatrix)
  }
  if(nCores>1)
  {
    if(!mctype %in% c("domc")) stop("Valid 'mctype': 'doMC'") # 'snow'  
    if(mctype=="snow")
    {
      # library(parallel)
      # library(doSNOW)
      # library(plyr)
      cl <- parallel::makeCluster(nCores, type = "SOCK")
      doSNOW::registerDoSNOW(cl)
      if(verbose) message("Using ", length(cl), " cores with SNOW.")
      # clusterEvalQ(cl, library(AUCell))
      parallel::clusterExport(cl, c("geneSets", "rankings", "aucMaxRank", ".AUC.geneSet", ".auc", ".AUC.geneSet_norm", ".AUC.geneSet_old"), envir=environment())
      opts <- list(preschedule=TRUE)
      # clusterSetRNGStream(cl, seed)
      aucMatrix <- suppressWarnings(plyr::llply(.data=names(geneSets), 
                               .fun=function(gSetName) 
                                 setNames(list(.AUC.geneSet(geneSet=geneSets[[gSetName]], rankings=rankings, aucMaxRank=aucMaxRank, gSetName=gSetName)), gSetName),
                               .parallel=TRUE, .paropts=list(.options.snow=opts), .inform=FALSE))
      aucMatrix <- do.call(rbind,
                           unlist(aucMatrix, recursive=FALSE))
      parallel::stopCluster(cl)
    }
    if(mctype=="domc")
    {
     doMC::registerDoMC(nCores)
     if(verbose)
       message("Using ", foreach::getDoParWorkers(), " cores with doMC.")
  
     aucMatrix <- foreach::"%dopar%"(foreach::foreach(gSetName=names(geneSets)),
      {
        setNames(list(.AUC.geneSet(geneSet=geneSets[[gSetName]],
                                   rankings=rankings, aucMaxRank=aucMaxRank,
                                   gSetName=gSetName)),
                 gSetName)
      })
      aucMatrix <- do.call(rbind,
                           unlist(aucMatrix, recursive=FALSE))
    }
  }

  aucMatrix <- aucMatrix[intersect(names(geneSets), rownames(aucMatrix)),,drop=FALSE]
  missingSets <- names(geneSets)[which(!names(geneSets) %in% rownames(aucMatrix))]
  if(length(missingSets)>0) warning("The AUC for the following sets was not calculated: ", paste(missingSets, collapse=", "))
  
  ######################################################################
  ##### Messages for missing genes
  missingGenes <- as.matrix(aucMatrix[,c("missing", "nGenes") , drop=FALSE])
  missingPercent <- as.numeric(missingGenes[,"missing"])/
    as.numeric(missingGenes[,"nGenes"])
  missingPercent <- setNames(missingPercent, rownames(missingGenes))

  if(all(missingPercent>=.80))
    stop("Fewer than 20% of the genes in the gene sets are included in the rankings.",
         "Check wether the gene IDs in the 'rankings' and 'geneSets' match.")

  if(any(missingPercent>.80))
  {
    warning("The following gene sets will be excluded from the analysis",
            "(less than 20% of their genes are available):\n",
            paste(names(missingPercent)[which(missingPercent >= .80)],
                  collapse=", "),
            sep="", immediate.=TRUE)
    aucMatrix <- aucMatrix[which(missingPercent < .80),,drop=FALSE]
  }

  missingGenes <- missingGenes[rownames(aucMatrix),,drop=FALSE]
  if(sum(missingGenes[,"missing"])>0)
  {
    msg1 <- "Genes in the gene sets NOT available in the dataset: \n"
    msg2 <-  sapply(rownames(missingGenes)[which(missingGenes[,"missing"]>0.01)],
                    function(gSetName){
                      paste("\t", gSetName, ": \t",
                            missingGenes[gSetName,"missing"],
                            " (",round(missingPercent[gSetName]*100),"% of ",
                            missingGenes[gSetName,"nGenes"],")",sep="")
                    })
    if(verbose)
      message(msg1, paste(msg2, collapse="\n"), sep="")
  }
  # (remove missing genes info from AUC matrix)
  aucMatrix <- aucMatrix[,1:(ncol(aucMatrix)-2), drop=FALSE]

  ######################################################################
  #### End: Return
  names(dimnames(aucMatrix)) <- c("gene sets", "cells")
  new("aucellResults",
      SummarizedExperiment::SummarizedExperiment(assays=list(AUC=aucMatrix)))
}

# add?: AUCThreshold
.AUC.geneSet_old <- function(geneSet, rankings, aucMaxRank, gSetName="")
{
  geneSet <- unique(geneSet)
  nGenes <- length(geneSet)
  geneSet <- geneSet[which(geneSet %in% rownames(rankings))]
  missing <- nGenes-length(geneSet)

  gSetRanks <- rankings[which(rownames(rankings) %in% geneSet),,drop=FALSE]
  rm(rankings)

  aucThreshold <- round(aucMaxRank)
  maxAUC <- aucThreshold * nrow(gSetRanks)  # database.gene_count 

  # Apply by columns (i.e. to each ranking)
  auc <- apply(gSetRanks, 2, .auc, aucThreshold, maxAUC)

  return(c(auc, missing=missing, nGenes=nGenes))
}

.AUC.geneSet_norm <- function(geneSet, rankings, aucMaxRank, gSetName="")
{
  geneSet <- unique(geneSet)
  nGenes <- length(geneSet)
  geneSet <- geneSet[which(geneSet %in% rownames(rankings))]
  missing <- nGenes-length(geneSet)
  
  gSetRanks <- rankings[which(rownames(rankings) %in% geneSet),,drop=FALSE]
  rm(rankings)
  
  aucThreshold <- round(aucMaxRank)
  ########### NEW version:  #######################
  x_th <- 1:nrow(gSetRanks)
  x_th <- sort(x_th[x_th<aucThreshold])
  y_th <- seq_along(x_th)
  maxAUC <- sum(diff(c(x_th, aucThreshold)) * y_th) 
  ############################################
  
  # Apply by columns (i.e. to each ranking)
  auc <- apply(gSetRanks, 2, .auc, aucThreshold, maxAUC)
  
  return(c(auc, missing=missing, nGenes=nGenes))
}


# oneRanking <- gSetRanks[,3, with=FALSE]
.auc <- function(oneRanking, aucThreshold, maxAUC)
{
  x <- unlist(oneRanking)

  x <- sort(x[x<aucThreshold])
  y <- seq_along(x)
  sum(diff(c(x, aucThreshold)) * y)/maxAUC
}
