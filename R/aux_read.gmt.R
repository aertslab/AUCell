# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @title Functions to manipulate GMT files (gene-set format)
#' @description read.gmt: Reads the gene-sets in a .gmt file.
#'
#' write.gmt: Write several gene-sets into a single file (input: list).
#'
#' write.gmt.single: Writes a gene set into a file (input: character vector).
#' @param geneSets [write.gmt] List containing the gene-sets. The gene-sets should be provided as a 'named list' in which each element is a gene-set (i.e. \code{list(geneSet1=c("gene1", "gene2"))})
#' @param fileName File name to read or write into
#' @param gmtDir [write.gmt] Directory, if diferent from \code{getwd()}
#' @param genes [write.gmt.single] Character vector containing the genes to write into the file
#' @param setName [write.gmt.single] Name of the gene set
#' @return \code{read.gmt} returns a list of gene sets (the gene-set names as \code{names()}, and the elements of the gene-set as character vector).
#' @details File format: One gene-set per line, tab separated elements (i.e. genes). The first two elements of the line are the gene-set name (The second element can be replaced by a description of the gene set).
#'
#' Example:
#'
#' <geneSetName> <tab> <geneSetName> <tab> <tab-separated genes>
#' @example inst/examples/example_gmt.R


#' @export
read.gmt <- function(fileName) {
  tmp <- unname(sapply(readLines(fileName), function(x) strsplit(as.character(x),"\t")))
  tmp <- tmp[which(lengths(tmp)>0)]
  names(tmp) <- sapply(tmp, function(x) x[1])
  lapply(tmp, function(x) x[3:length(x)])
}

#' @rdname read.gmt
#' @export
write.gmt <- function(geneSets, fileName=NULL, gmtDir=NULL)
{
  geneSetList <- geneSets
  currDir <- getwd()
  if(!is.null(gmtDir)) setwd(gmtDir)

  if(is.null(fileName))
  {
    lapply(names(geneSetList), function(gsName) write.gmt(geneSetList[[gsName]], gsName))
  }else{

    gmtLines <- sapply(names(geneSetList), function(gsName){
      paste(c(rep(gsName, 2), geneSetList[[gsName]]), collapse="\t")
    })
    writeLines(gmtLines, con=fileName)
  }

  setwd(currDir)
}

#' @rdname read.gmt
#' @export
write.gmt.single <- function(genes, setName, fileName=NULL)
{
  gmtLine <- paste(c(rep(setName, 2), genes), collapse="\t")

  if(is.null(fileName)) fileName <- paste(setName, ".gmt", sep="")
  writeLines(gmtLine, con=fileName)
}


# merge.gmt <- function(signaturesDir, gmtFileName="mergedGmts.gmt"){
#   gmtFiles <- list.files(signaturesDir)[grep(".gmt", list.files(signaturesDir))]
#   gmtFiles <- sapply(gmtFiles, function(gtmFile) readLines(paste(signaturesDir, gtmFile, sep="/")))
#   writeLines(gmtFiles, con=gmtFileName)
#   return(gmtFileName)
# }
