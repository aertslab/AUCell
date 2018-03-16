## AUCell_plot()
##
# Input:
# Output:

test_AUCell_exploreThresholds <- function()
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
  geneSets <- GSEABase::GeneSetCollection(list(
    GSEABase::GeneSet(sample(rownames(exprMatrix), 10), setName="gs1"),
    GSEABase::GeneSet(sample(rownames(exprMatrix), 10), setName="gs2")))
  cellsAUC <- suppressWarnings(AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=5))
  ##################################################

  thresholds <- AUCell_exploreThresholds(cellsAUC, plotHist=FALSE)

  .check_AUCell_exploreThresholds(thresholds, geneSets)

  ### Multicore (fails check)
  # set.seed(123)
  # thresholds_multicore_1 <- AUCell_exploreThresholds(cellsAUC, plotHist=FALSE, nCores=2)
  # set.seed(123)
  # thresholds_multicore_2 <- AUCell_exploreThresholds(cellsAUC, plotHist=FALSE, nCores=2)
  #
  # testthat::expect_equal(thresholds_multicore_1, thresholds_multicore_2)
  .check_AUCell_exploreThresholds(thresholds, geneSets)
}

.check_AUCell_exploreThresholds <- function(thresholds, geneSets)
{
  testthat::expect_equal(length(thresholds), 2)
  testthat::expect_equal(names(thresholds), names(geneSets))
  testthat::expect_equal(names(thresholds[[1]]$aucThr), c("selected", "thresholds", "comment"))
}

test_that("AUCell_exploreThresholds tests", test_AUCell_exploreThresholds())
