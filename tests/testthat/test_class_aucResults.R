## Class: aucellResults
##
# Input:
# Output:

test_aucellResults <- function()
{
  ##################################################
  # Create objects
  set.seed(123)
  exprMatrix <- matrix(data=sample(c(rep(0, 5000), sample(1:3, 5000, replace=TRUE))),
                       nrow=20)
  rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
  colnames(exprMatrix) <- paste("Cell", 1:500, sep="")

  cells_rankings <- suppressWarnings(AUCell_buildRankings(exprMatrix, plotStats=FALSE, verbose=FALSE))
  fewGenes <- sample(rownames(exprMatrix), 10)
  cellsAUC <- suppressWarnings(AUCell_calcAUC(fewGenes, cells_rankings, aucMaxRank=5))
  ##################################################

  # Ranking
  testthat::expect_equal(class(cells_rankings)[1], "aucellResults")
  testthat::expect_equal(dim(cells_rankings), c(20,500))

  testthat::expect_equal(dim(getRanking(cells_rankings)), c(20,500))
  testthat::expect_error(getAUC(cells_rankings))

  # AUC
  testthat::expect_equal(class(cellsAUC)[1], "aucellResults")
  testthat::expect_equal(dim(cellsAUC), c(1,500))

  testthat::expect_equal(dim(getAUC(cellsAUC)), c(1,500))
  testthat::expect_error(getRanking(cellsAUC))
}
