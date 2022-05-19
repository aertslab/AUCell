## AUCell_buildRankings()
##
# Input:
# Output:

test_AUCell_buildRankings <- function()
{
  library(AUCell)
  ##################################################
  ### Fake dataset
  set.seed(123)
  exprMatrix <- matrix(data=sample(c(rep(0, 5000), sample(1:3, 5000, replace=TRUE))),
                       nrow=20)
  rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
  colnames(exprMatrix) <- paste("Cell", 1:500, sep="")
  ##################################################

  cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE, verbose=FALSE)

  .check_AUCell_buildRankings_output(cells_rankings)

  ### Other input classes:
  eset <- Biobase::ExpressionSet(assayData=exprMatrix)
  rEset <- AUCell_buildRankings(eset, plotStats=FALSE)
  testthat::expect_equal(class(rEset)[1], "aucellResults")

  sexp <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=exprMatrix))
  rSexp <- AUCell_buildRankings(sexp, plotStats=FALSE)
  testthat::expect_equal(class(rSexp)[1], "aucellResults")

  sparseMat <- as(exprMatrix, "dgCMatrix")
  rSparse <- AUCell_buildRankings(sparseMat, plotStats=FALSE, splitByBlocks=TRUE)
  testthat::expect_equal(class(rSparse)[1], "aucellResults")
}

.check_AUCell_buildRankings_output <- function(cells_rankings)
{
  testthat::expect_equal(ncol(cells_rankings), 500)
  testthat::expect_equal(nrow(cells_rankings), 20)
  
  testthat::expect_equal(getRanking(cells_rankings)['Gene13','Cell2'], 1)
  

  testthat::expect_equal(class(cells_rankings)[1], "aucellResults")
  testthat::expect_equal(SummarizedExperiment::assayNames(cells_rankings)[1], "ranking")
}


test_that("AUCell_buildRankings tests", test_AUCell_buildRankings())
