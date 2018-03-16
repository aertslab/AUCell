## AUCell_calcAUC()
##
# Input:
# Output:

test_AUCell_calcAUC <- function()
{
  library(AUCell)
  ##################################################
  # Previous steps
  set.seed(123)
  exprMatrix <- matrix(data=sample(c(rep(0, 5000), sample(1:3, 5000, replace=TRUE))),
                       nrow=20)
  rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
  colnames(exprMatrix) <- paste("Cell", 1:500, sep="")

  cells_rankings <- suppressWarnings(AUCell_buildRankings(exprMatrix, plotStats=FALSE, verbose=FALSE))
  ##################################################

  # a) Character vector (i.e. only one gene-set)
  fewGenes <- sample(rownames(exprMatrix), 10)
  cellsAUC <- suppressWarnings(AUCell_calcAUC(fewGenes, cells_rankings, aucMaxRank=5))

  testthat::expect_equal(nrow(cellsAUC), 1)
  testthat::expect_equal(ncol(cellsAUC), 500)

  # b) List
  otherGenes <- sample(rownames(exprMatrix), 5)
  geneSets <- list(geneSet1=fewGenes,
                   geneSet2=otherGenes)
  cellsAUC <- suppressWarnings(AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=5))

  testthat::expect_equal(nrow(cellsAUC), 2)
  testthat::expect_equal(ncol(cellsAUC), 500)

  # c) GeneSet object (from GSEABase)
  geneSetOne <- GSEABase::GeneSet(fewGenes, setName="geneSetOne")
  cellsAUC <- suppressWarnings(AUCell_calcAUC(geneSetOne, cells_rankings, aucMaxRank=5))

  testthat::expect_equal(nrow(cellsAUC), 1)
  testthat::expect_equal(ncol(cellsAUC), 500)

  # d) GeneSetCollection object (from GSEABase)
  geneSetTwo <- GSEABase::GeneSet(otherGenes, setName="geneSetTwo")
  geneSets <- GSEABase::GeneSetCollection(geneSetOne, geneSetTwo)
  cellsAUC <- suppressWarnings(AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=5))

  testthat::expect_equal(nrow(cellsAUC), 2)
  testthat::expect_equal(ncol(cellsAUC), 500)

  testthat::expect_equal(class(cellsAUC)[1], "aucellResults")
  testthat::expect_equal(SummarizedExperiment::assayNames(cellsAUC)[1], "AUC")

  ### Multicore
  # cellsAUC_multicore <- suppressWarnings(AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=5, nCores=2))
  # testthat::expect_equal(getAUC(cellsAUC), getAUC(cellsAUC_multicore))
}

test_that("AUCell_calcAUC tests", test_AUCell_calcAUC())
