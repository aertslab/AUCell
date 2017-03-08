
#' @title Class to store AUC
#' @description Class with only one slot: @AUC, a matrix which contains the AUC for the gene-sets & cells.
#' @examples
#' fakeAUC <- cellsAUC(AUC=matrix(1/1:30, nrow=3)) # gene-sets by cells.
#' fakeAUC
#'
#' nGeneSets(fakeAUC)
#' nCells(fakeAUC)
#' cellNames(fakeAUC)
#' geneSetNames(fakeAUC)
#'
#' class(fakeAUC)
#' class(fakeAUC@AUC)
#' dim(fakeAUC@AUC)
#'
#' # To subset & access the AUC slot:
#' fakeAUC[1:2,]
#' fakeAUC[,3:4]
#' auc(fakeAUC)[,1:5] # fakeAUC@AUC[,1:5]
#'
#' # To subset the AUC slot, but keeping the class:
#' subset(fakeAUC, 1:2, "geneSets")
#' subset(fakeAUC, 3:4, "cells")
#' @export

cellsAUC <- setClass(
  # Set the name for the class
  Class="cellsAUC",

  # Define the slots
  slots = c(
    AUC = "matrix"
  )
)

setMethod("show",
          signature="cellsAUC",
          definition = function(object) {
            message(paste("AUC for",
                          nrow(object@AUC), "gene-sets and", ncol(object@AUC), "cells"))
            if(!is.null(rownames(object@AUC))) message(paste("First gene-sets:", paste(head(rownames(object@AUC), 3), collapse=", ")))
            if(!is.null(colnames(object@AUC))) message(paste("First cells:", paste(head(colnames(object@AUC), 3), collapse=", ")))
          }
)

setMethod('[', signature(x="cellsAUC"),
          definition=function(x, i, j, drop){
              x@AUC[i,j, drop=drop]
          })

setMethod("subset",
          signature="cellsAUC",
          definition = function(x, elements, select=c("geneSets", "cells")[1]) {
              if(length(select) > 1) stop()

              if(grepl("geneset", tolower(select))) x@AUC <- x@AUC[elements,, drop=FALSE]
              if(grepl("cell", tolower(select))) x@AUC <- x@AUC[,elements, drop=FALSE]

              x
          }
)

#' @export
setGeneric(name="auc", def=function(object) standardGeneric("auc"))
setMethod("auc",
          signature="cellsAUC",
          definition = function(object) {
            object@AUC
          }
)

#' @export
setGeneric(name="geneSetNames", def=function(object) standardGeneric("geneSetNames"))
setMethod("geneSetNames",
          signature="cellsAUC",
          definition = function(object) {
            rownames(object@AUC)
            }
)

#' @export
setGeneric(name="cellNames", def=function(object) standardGeneric("cellNames"))
setMethod("cellNames",
          signature="cellsAUC",
          definition = function(object) {
            colnames(object@AUC)
          }
)

#' @export
setGeneric(name="nGeneSets", def=function(object) standardGeneric("nGeneSets"))
setMethod("nGeneSets",
          signature="cellsAUC",
          definition = function(object) {
              nrow(object@AUC)
          }
)

#' @export
setGeneric(name="nCells", def=function(object) standardGeneric("nCells"))
setMethod("nCells",
          signature="cellsAUC",
          definition = function(object) {
              ncol(object@AUC)
          }
)


