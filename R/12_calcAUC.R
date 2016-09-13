
# AUCell(rankings, geneSets)
# seed: add to return
# plotHist: or filename?
# AUCell.calcAUC <- function(geneSets, rankings, nCores=4, aucMaxRank=0.05*nrow(rankings), verbose=TRUE)
# {
#   calcAUC(geneSets=geneSets, rankings=rankings, nCores=nCores, aucMaxRank=aucMaxRank, verbose=verbose)
# }


# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
#' @import data.table
#'
#' @title to do
#' @description to do
#' @param to do
#' @return to do
#' @examples  #to do
#' @export
AUCell.calcAUC <- function(geneSets, rankings, nCores=1, aucMaxRank=0.03*nrow(rankings), verbose=TRUE) #, seed=123, plotHist=TRUE
{
  if(!is.list(geneSets)) stop("geneSets should be a named list.")
  if(is.null(names(geneSets))) stop("geneSets should be a named list.")
  if(nCores > length(geneSets)) nCores <- length(geneSets) # No point in using more...

  if(!is.data.table(rankings)) stop("Rankings should be a data.table (i.e. genes x [cells or motifs])")
  # if(!key(rankings) == "rn") stop("The rankings key should be 'rn'.")

  ######################################################################
  #### 1. Calculate the AUC for each gene set
  if(nCores==1)
  {
    aucMatrix <- sapply(names(geneSets), function(gSetName)  .AUC.geneSet(geneSet=geneSets[[gSetName]], rankings=rankings, aucMaxRank=aucMaxRank, gSetName=gSetName))
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
    aucMatrix <- do.call(cbind, unlist(aucMatrix, recursive = FALSE)[names(geneSets)])
  }

  ######################################################################
  ##### Messages for missing genes
  missinGenes <- aucMatrix[c("missing", "nGenes"), , drop=FALSE]
  missingPercent <- missinGenes["missing",, drop=FALSE]/missinGenes["nGenes",, drop=FALSE]
  if(all(missingPercent>=.80)) stop("Fewer than 20% of the genes in the gene sets are included in the rankings. Check wether the gene IDs in the 'rankings' and 'geneSets' match.")

  if(any(missingPercent>.80))
  {
    warning(paste("The following gene sets will be excluded from the analysis (less than 20% of their genes are available):\n",
     paste(names(missingPercent)[which(missingPercent >= .80)], collapse=", "), sep=""), immediate.=TRUE)
    aucMatrix <- aucMatrix[,which(missingPercent < .80),drop=FALSE]
  }

  if(sum(missinGenes["missing",])>0)
  {
    msg1 <- "Genes in the gene sets NOT available in the dataset: \n"
    msg2 <-  sapply(colnames(missinGenes)[which(missinGenes["missing",]>0)], function(gSetName)
        paste("\t", gSetName, ": \t", missinGenes["missing",gSetName],
          " (",round(missingPercent[,gSetName]*100),"% of ", missinGenes["nGenes",gSetName],")",sep=""))
    if(verbose) message(paste(msg1, paste(msg2, collapse="\n"), sep=""))
  }
  # (remove missing genes info from AUC matrix)
  aucMatrix <- aucMatrix[1:(nrow(aucMatrix)-2),, drop=FALSE]


  ######################################################################
  #### End: Return
  aucMatrix
}

.AUC.geneSet <- function(geneSet, rankings, aucMaxRank, gSetName="")  # add?: AUCThreshold
{
  geneSet <- unique(geneSet)
  nGenes <- length(geneSet)
  geneSet <- geneSet[which(geneSet %in% rankings$rn)]
  missing <- nGenes-length(geneSet)

  # stop(paste(c(sum(!geneSet %in% rankings$rn), class(geneSet)), collapse=", "))
  # gSetRanks <- rankings[geneSet][,-"rn", with=FALSE] # gene names are no longer needed
  gSetRanks <- subset(rankings, rn %in% geneSet)[,-"rn", with=FALSE] # gene names are no longer needed
  rm(rankings)

  aucThreshold <- round(aucMaxRank)            # 5%: same as .ini, but Stein said 3%
  maxAUC <- aucThreshold * nrow(gSetRanks)     # database.gene_count  ->  IS THIS CORRECT?

  # Apply by columns (i.e. to each ranking)
  auc <- sapply(gSetRanks, .calcAUC, aucThreshold, maxAUC)

  c(auc, missing=missing, nGenes=nGenes)
}


# oneRanking <- gSetRanks[,3, with=F]
.calcAUC <- function(oneRanking, aucThreshold, maxAUC)
{
  x <- unlist(oneRanking)

  x <- sort(x[x<aucThreshold])
  y <- 1:length(x)
  sum(diff(c(x, aucThreshold)) * y)/maxAUC
}
