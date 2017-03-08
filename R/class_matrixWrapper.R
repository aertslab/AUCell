
#' @title Wrapper to the matrix that stores the AUC or the cell rankings.
#' @description
#' This class is a wrapper for the AUC matrix and cell rankings.
#'
#' The matrices are stored in the slot @matrix slot with the matrix, but they can be accessed through the regular methods (i.e. nrow, rownames... )
#'
#' Types:
#'  - "AUC": The matrix contains the AUC for the gene-sets (or region-sets) & cells.
#'  - "ranking": The @matrix contains the gene rankings for each cell.
#'
#' Methods: See examples.
#' @import BiocGenerics
#' @example inst/examples/example_class_matrixWrapper.R
#' @export matrixWrapper
#' @exportClass matrixWrapper
matrixWrapper <- setClass(
  # Set the name for the class
  Class="matrixWrapper",

  # Define the slots
  slots = c(
    matrix = "matrix",
    colType = "character", # cell
    rowType = "character", # gene/region (ranking) or gene-set/region-set (AUC)
    matrixType= "character"   # AUC or Ranking
  )
)

#' @export
setMethod("show",
          signature="matrixWrapper",
          definition = function(object) {
            message <- paste(object@matrixType, " for ", nrow(object@matrix)," ", object@rowType, "s (rows) and ", ncol(object@matrix)," ", object@colType, "s (columns).\n", sep="")

            if(object@matrixType == "AUC") message <- paste(message, "\nTop-left corner of the AUC matrix:\n", sep="")
            if(object@matrixType == "Ranking") message <- paste(message, "\nTop-left corner of the ranking:\n", sep="")

            cat(message)
            show(head(object@matrix[,0:min(5, ncol(object@matrix)),drop=FALSE]))
          }
)

#' @export
setMethod('[', signature(x="matrixWrapper"),
          definition=function(x, i, j, drop){
            x@matrix[i,j, drop=drop]
          })

#' @export
setMethod("subset",
          signature="matrixWrapper",
          definition = function(x, elements, select=x@rowType) {

            if(length(select) > 1) stop()

            if(grepl("cell", tolower(select))) {
              x@matrix <- x@matrix[,elements, drop=FALSE]
            }else{
              x@matrix <- x@matrix[elements,, drop=FALSE]
            }
            x
          }
)

#' @export
setMethod("dim",
          signature="matrixWrapper",
          definition = function(x) {
            dim(x@matrix)
          }
)

#' @export
setMethod("rownames",
          signature="matrixWrapper",
          definition = function(x) {
            rownames(x@matrix)
          }
)

#' @export
setMethod("colnames",
          signature="matrixWrapper",
          definition = function(x) {
            colnames(x@matrix)
          }
)

#' @export
setMethod("nrow",
          signature="matrixWrapper",
          definition = function(x) {
            nrow(x@matrix)
          }
)

#' @export
setMethod("ncol",
          signature="matrixWrapper",
          definition = function(x) {
            ncol(x@matrix)
          }
)


##### Access the matrix:
#' @export
setGeneric(name="getAuc", def=function(object) standardGeneric("getAuc"))
setMethod("getAuc",
          signature="matrixWrapper",
          definition = function(object) {
            if(grepl("auc", tolower(object@matrixType))) {
              object@matrix
            }else{
              stop("This object does not contain an AUC matrix.")
            }
          }
)

#' @export
setGeneric(name="getRanking", def=function(object) standardGeneric("getRanking"))
setMethod("getRanking",
          signature="matrixWrapper",
          definition = function(object) {
            if(grepl("ranking", tolower(object@matrixType))) {
              object@matrix
            }else{
              stop("This object does not contain a ranking.")
            }
          }
)

