## AUCell_plot()
##
# Input:
# Output:

test_AUCell_plotHist <- function()
{
  # TO DO: add tests for _plotTsne()s
  library(AUCell)
  ##################################################
  # Previous steps
  set.seed(123)
  exprMatrix <- matrix(data=sample(c(rep(0, 5000), sample(1:3, 5000, replace=TRUE))),
                       nrow=20)
  rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
  colnames(exprMatrix) <- paste("Cell", 1:500, sep="")

  cells_rankings <- suppressWarnings(AUCell_buildRankings(exprMatrix, plotStats=FALSE, verbose=FALSE))
  fewGenes <- sample(rownames(exprMatrix), 10)
  cellsAUC <- suppressWarnings(AUCell_calcAUC(fewGenes, cells_rankings, aucMaxRank=5))
  ##################################################

  testthat::expect_equal(length(AUCell_plotHist(cellsAUC)), 1)
  testthat::expect_equal(length(AUCell_plotHist(cellsAUC[1,])), 1)
}

test_that("AUCell_plotHist tests", test_AUCell_plotHist())
