
#' @title GeneSet methods
#' @description Functions to manipulate GeneSet and GeneSetCollection objects (from package GSEABase)
#' @import GSEABase
#' @param geneSet One gene-set (\code{\link[GSEABase]{GeneSet}})
#' @param geneSets Gene-set collection (\code{\link[GSEABase]{GeneSetCollection}})
#' @param geneNames Gene names (for subset)
#' @param newNames New names (to assign to the gene sets)
#' @return
#' - **nGenes()**: provides the number of genes in the gene-set, or each of the gene-sets in a collection
#'
#' - **subsetGeneSets()**: Subsets each of the gene-sets in a collection to contain only the genes inthe given list. Equivalent to intersect(), but keeping the original gene-set name.
#'
#' - **setGeneSetNames()**: Modifies the name of each gene-set in a collection
#' @docType methods
#' @rdname GeneSet-methods
#'
#' @example
#'  genes_1 <- GeneSet(paste("Gene", 1:20, sep=""), setName="geneSet1")
#'  genes_2 <- GeneSet(paste("Gene", 18:22, sep=""), setName="geneSet2")
#'  geneSets <- GeneSetCollection(genes_1, genes_2)
#'
#'  nGenes(genes_1)
#'  nGenes(geneSets)
#'
#'  subsetGeneSets(geneSets, paste("Gene", 15:20, sep=""))
#'
#'  geneSets_newNames <- setGeneSetNames(geneSets, c("one", "two"))
#'  names(geneSets_newNames)
#'
#' @export
setGeneric("nGenes", signature="geneSet",
           function(geneSet)
           {
             standardGeneric("nGenes")
           })

#' @rdname GeneSet-methods
#' @aliases nGenes,GeneSet-method
setMethod("nGenes", "GeneSet",
  function(geneSet)
  {
    length(geneIds(geneSet))
  })

#' @rdname GeneSet-methods
#' @aliases nGenes,GeneSetCollection-method
setMethod("nGenes", "GeneSetCollection",
  function(geneSet)
  {
    setNames(sapply(geneSet, nGenes), names(geneSet))
  })

# equivalent to intersect, but keeping original name
#' @rdname GeneSet-methods
#' @export
setGeneric("subsetGeneSets",
           function(geneSets, geneNames)
           {
             standardGeneric("subsetGeneSets")
           })

#' @rdname GeneSet-methods
#' @aliases subsetGeneSets,GeneSetCollection-method
setMethod("subsetGeneSets", signature="GeneSetCollection",
          function(geneSets, geneNames)
          {
            GeneSetCollection(lapply(geneSets, function(x) {
              y <- x & geneNames
              setName(y) <- setName(x)
              y
            }))
          })

#' @rdname GeneSet-methods
#' @export
setGeneric("setGeneSetNames",
           function(geneSets, newNames)
           {
             standardGeneric("setGeneSetNames")
           })


#' @rdname GeneSet-methods
#' @aliases setGeneSetNames,GeneSetCollection-method
setMethod("setGeneSetNames", signature="GeneSetCollection",
          function(geneSets, newNames)
          {
            if(length(geneSets) != length(newNames))
              stop("The number of gene sets and new names do not match.")

            geneSets <- unlist(geneSets)
            for(i in seq_len(length(newNames)))
            {
              setName(geneSets[[i]]) <- newNames[i]
            }
            GeneSetCollection(geneSets)
          })
