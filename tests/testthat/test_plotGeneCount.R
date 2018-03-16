## AUCell_plot()
##
# Input:
# Output:

test_plotGeneCount <- function()
{
  library(AUCell)
  ##################################################
  # Previous steps
  set.seed(123)
  exprMatrix <- matrix(data=sample(c(rep(0, 5000), sample(1:3, 5000, replace=TRUE))),
                       nrow=20)
  rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
  colnames(exprMatrix) <- paste("Cell", 1:500, sep="")

  testOutput <- plotGeneCount(exprMatrix, verbose=FALSE)
  ##################################################

  testthat::expect_equal(class(testOutput), "numeric")
  testthat::expect_equal(length(testOutput), 6)
}

test_that("plotGeneCount tests", test_plotGeneCount())
