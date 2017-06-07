# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
#' @import data.table
#'
#' @title AUCell.calcAUC
#' @description Calculates the 'AUC' for each gene-set in each cell.
#' @param geneSets List of gene-sets (or signatures) to test in the cells. The gene-sets should be provided as a 'named list' in which each element is a gene-set (i.e. \code{list(geneSet1=c("gene1", "gene2"))})
#' @param rankings 'Rankings' created for this dataset with \code{\link{AUCell.buildRankings}}.
#' @param nCores Number of cores to use for computation.
#' @param aucMaxRank Threshold to calculate the AUC (see 'details' section).
#' @param verbose Should the function show progress messages? (TRUE / FALSE)
#' @return Matrix with the AUC values (gene-sets as rows, cells as columns).
#' @details In a simplified way, the AUC value represents the fraction of genes, within the top X genes in the ranking, that are included in the signature.
#' The parameter 'aucMaxRank' allows to modify the number of genes (maximum ranking) that is used to perform this computation.
#' By default, it is set to 5\% of the total number of genes in the rankings. Common values may range from 1 to 20\%.
#' @seealso Previous step in the workflow: \code{\link{AUCell.buildRankings}}. Next step in the workflow: \code{\link{AUCell.exploreThresholds}}.
#'
#' See the package vignette for examples and more details: \code{vignette("AUCell")}
#' @example inst/examples/example_AUCell.calcAUC.R
#' @export
AUCell.calcAUC <- function(geneSets, rankings, nCores=1, aucMaxRank=ceiling(0.05*nrow(rankings)), verbose=TRUE) #, seed=123, plotHist=TRUE
{
  if(!is.list(geneSets)) stop("geneSets should be a named list.")
  if(is.null(names(geneSets))) stop("geneSets should be a named list.")
  if(nCores > length(geneSets)) nCores <- length(geneSets) # No point in using more...
  # if(aucMaxRank < 300) warning(paste("Using only the first ", aucMaxRank, " genes (aucMaxRank) to calculate the AUC."))
  if(aucMaxRank < 300) warning(paste("Using only the first ", aucMaxRank, " genes (aucMaxRank) to calculate the AUC."), immediate.=TRUE)

  if(class(rankings)!="matrixWrapper"){
    stop("Rankings should be a the object returned by AUCell.buildRankings()")
  } else{
    rankings <- getRanking(rankings)
  }
  # if(!key(rankings) == "rn") stop("The rankings key should be 'rn'.")

  ######################################################################
  #### 1. Calculate the AUC for each gene set
  if(nCores==1)
  {
    aucMatrix <- sapply(names(geneSets), function(gSetName)
      .AUC.geneSet(geneSet=geneSets[[gSetName]], rankings=rankings, aucMaxRank=aucMaxRank, gSetName=gSetName))
    aucMatrix <- t(aucMatrix)
  }else
  {
    # Run each geneSet in parallel
    suppressMessages(require("doMC", quietly=TRUE))
    doMC::registerDoMC(nCores)
    if(verbose) message(paste("Using", getDoParWorkers(), "cores."))

    aucMatrix <- foreach(gSetName=names(geneSets)) %dopar%
    {
      setNames(list(.AUC.geneSet(geneSet=geneSets[[gSetName]], rankings=rankings, aucMaxRank=aucMaxRank, gSetName=gSetName)), gSetName)
    }
    aucMatrix <- do.call(rbind, unlist(aucMatrix, recursive = FALSE)[names(geneSets)])
  }

  ######################################################################
  ##### Messages for missing genes
  missingGenes <- as.matrix(aucMatrix[,c("missing", "nGenes") , drop=FALSE])
  missingPercent <- as.numeric(missingGenes[,"missing"])/as.numeric(missingGenes[,"nGenes"])
  missingPercent <- setNames(missingPercent, rownames(missingGenes))

  if(all(missingPercent>=.80)) stop("Fewer than 20% of the genes in the gene sets are included in the rankings. Check wether the gene IDs in the 'rankings' and 'geneSets' match.")

  if(any(missingPercent>.80))
  {
    warning(paste("The following gene sets will be excluded from the analysis (less than 20% of their genes are available):\n",
     paste(names(missingPercent)[which(missingPercent >= .80)], collapse=", "), sep=""), immediate.=TRUE)
    aucMatrix <- aucMatrix[which(missingPercent < .80),,drop=FALSE]
  }

  if(sum(missingGenes[,"missing"])>0)
  {
    msg1 <- "Genes in the gene sets NOT available in the dataset: \n"
    msg2 <-  sapply(rownames(missingGenes)[which(missingGenes[,"missing"]>0)], function(gSetName){
        paste("\t", gSetName, ": \t", missingGenes[gSetName,"missing"],
          " (",round(missingPercent[gSetName]*100),"% of ", missingGenes[gSetName,"nGenes"],")",sep="")})
    if(verbose) message(paste(msg1, paste(msg2, collapse="\n"), sep=""))
  }
  # (remove missing genes info from AUC matrix)
  aucMatrix <- aucMatrix[,1:(ncol(aucMatrix)-2), drop=FALSE]


  ######################################################################
  #### End: Return
  matrixWrapper(matrix=aucMatrix, rowType="gene-set", colType="cell", matrixType="AUC")
}

.AUC.geneSet <- function(geneSet, rankings, aucMaxRank, gSetName="")  # add?: AUCThreshold
{
  geneSet <- unique(geneSet)
  nGenes <- length(geneSet)
  geneSet <- geneSet[which(geneSet %in% rownames(rankings))]
  missing <- nGenes-length(geneSet)

  # gSetRanks <- subset(rankings, rn %in% geneSet)[,-"rn", with=FALSE] # gene names are no longer needed
  gSetRanks <- rankings[which(rownames(rankings) %in% geneSet),,drop=FALSE]
  rm(rankings)

  aucThreshold <- round(aucMaxRank)
  maxAUC <- aucThreshold * nrow(gSetRanks)     # database.gene_count  ->  IS THIS CORRECT?
  # maxAUC <- sum(diff(c(1:min(nGenes,(aucThreshold-1)), aucThreshold)) * 1:(aucThreshold-1))

  # Apply by columns (i.e. to each ranking)
  auc <- apply(gSetRanks, 2, .calcAUC, aucThreshold, maxAUC)

  c(auc, missing=missing, nGenes=nGenes)
}

# oneRanking <- gSetRanks[,3, with=FALSE]
.calcAUC <- function(oneRanking, aucThreshold, maxAUC)
{
  x <- unlist(oneRanking)

  x <- sort(x[x<aucThreshold])
  y <- 1:length(x)
  sum(diff(c(x, aucThreshold)) * y)/maxAUC
}
