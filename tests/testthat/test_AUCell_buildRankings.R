## AUCell_buildRankings()
##
# Input:
# Output:

test_AUCell_buildRankings <- function()
{
  ##################################################
  ### Fake dataset
  set.seed(123)
  exprMatrix <- matrix(data=sample(c(rep(0, 5000), sample(1:3, 5000, replace=TRUE))),
                       nrow=20)
  rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
  colnames(exprMatrix) <- paste("Cell", 1:500, sep="")
  ##################################################

  cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=T, verbose=FALSE)

  .check_AUCell_buildRankings_output(cells_rankings)

  ### Multicore
  set.seed(123)
  cells_rankings_multicore_1 <- AUCell_buildRankings(exprMatrix, plotStats=T, verbose=FALSE, nCores=4)
  set.seed(123)
  cells_rankings_multicore_2 <- AUCell_buildRankings(exprMatrix, plotStats=T, verbose=FALSE, nCores=4)

  testthat::expect_equal(getRanking(cells_rankings_multicore_1), getRanking(cells_rankings_multicore_2))
  .check_AUCell_buildRankings_output(cells_rankings_multicore_1)
}

.check_AUCell_buildRankings_output <- function(cells_rankings)
{
  testthat::expect_equal(ncol(cells_rankings), 500)
  testthat::expect_equal(nrow(cells_rankings), 20)

  testthat::expect_equal(class(cells_rankings)[1], "aucellResults")
  testthat::expect_equal(SummarizedExperiment::assayNames(cells_rankings)[1], "ranking")
}
