## AUCell_run()
##
# Input:
# Output:

test_AUCell_run <- function()
{
  library(AUCell)
  ##################################################
  # Data
  set.seed(123)
  exprMatrix <- matrix(data=sample(c(rep(0, 5000), sample(1:3, 5000, replace=TRUE))),
                       nrow=20)
  rownames(exprMatrix) <- paste("Gene", 1:20, sep="")
  colnames(exprMatrix) <- paste("Cell", 1:500, sep="")
  ##################################################
  
  ### Test input: Geneset
  # a) Character vector (i.e. only one gene-set)
  fewGenes <- sample(rownames(exprMatrix), 10)
  cellsAUC <- suppressWarnings(AUCell_run(exprMatrix, fewGenes, aucMaxRank=5))

  testthat::expect_equal(nrow(cellsAUC), 1)
  testthat::expect_equal(ncol(cellsAUC), 500)

  # b) List
  otherGenes <- sample(rownames(exprMatrix), 5)
  geneSets <- list(geneSet1=fewGenes,
                   geneSet2=otherGenes)
  cellsAUC <- suppressWarnings(AUCell_run(exprMatrix, geneSets, aucMaxRank=5))

  testthat::expect_equal(nrow(cellsAUC), 2)
  testthat::expect_equal(ncol(cellsAUC), 500)

  # c) GeneSet object (from GSEABase)
  geneSetOne <- GSEABase::GeneSet(fewGenes, setName="geneSetOne")
  cellsAUC <- suppressWarnings(AUCell_run(exprMatrix, geneSetOne, aucMaxRank=5))

  testthat::expect_equal(nrow(cellsAUC), 1)
  testthat::expect_equal(ncol(cellsAUC), 500)

  # d) GeneSetCollection object (from GSEABase)
  geneSetTwo <- GSEABase::GeneSet(otherGenes, setName="geneSetTwo")
  geneSets <- GSEABase::GeneSetCollection(geneSetOne, geneSetTwo)
  set.seed(123)
  cellsAUC <- suppressWarnings(AUCell_run(exprMatrix, geneSets, aucMaxRank=5))

  testthat::expect_equal(nrow(cellsAUC), 2)
  testthat::expect_equal(ncol(cellsAUC), 500)

  testthat::expect_equal(class(cellsAUC)[1], "aucellResults")
  testthat::expect_equal(SummarizedExperiment::assayNames(cellsAUC)[1], "AUC")
  
  ### Test input: ExprMatrix
  ### Other input classes:
  set.seed(123)
  sparseMat <- as(exprMatrix, "dgCMatrix")
  rSparse <- AUCell_run(sparseMat, geneSets, aucMaxRank=5)
  testthat::expect_equal(class(rSparse)[1], "aucellResults")
  
  set.seed(123)
  sexp <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=exprMatrix))
  rSexp <- AUCell_run(sexp, geneSets, aucMaxRank=5)
  testthat::expect_equal(class(rSexp)[1], "aucellResults")
}

test_that("AUCell_run tests", test_AUCell_run())
