## Methods created to manipulate GeneSets and GeneSetCollection
##
# Input:
# Output:

test_GeneSet_methods <- function()
{
  library(AUCell)
  ##################################################
  library(GSEABase)

  genes_1 <- GeneSet(paste("Gene", 1:20, sep=""), setName="geneSet1")
  genes_2 <- GeneSet(paste("Gene", 18:22, sep=""), setName="geneSet2")

  geneSets <- GeneSetCollection(genes_1, genes_2)

  ##################################################

  testthat::expect_equal(nGenes(genes_1), 20)
  testthat::expect_equal(nGenes(geneSets), setNames(c(20, 5), names(geneSets)))

  geneSets_subset <- subsetGeneSets(geneSets, paste("Gene", 15:20, sep=""))
  testthat::expect_equal(nGenes(geneSets_subset),
                         setNames(c(6, 3), names(geneSets)))
  testthat::expect_equal(length(unique(unlist(geneIds(geneSets_subset)))), 6)

  geneSets_newNames <- setGeneSetNames(geneSets, c("one", "two"))
  testthat::expect_equal(names(geneSets_newNames), c("one", "two"))
}

test_that("GeneSet_methods tests", test_GeneSet_methods())
