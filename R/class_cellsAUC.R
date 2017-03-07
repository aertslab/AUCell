
#' @title Class to store AUC
#' @description Class with only one slot: @AUC, a matrix which contains the AUC for the gene-sets & cells.
#' @example
#' fakeAUC <- motifsAUC(AUC=matrix(1/1:30, nrow=3)) # gene-sets by cells.
#' fakeAUC
#' dim(fakeAUC@AUC)
#' class(fakeAUC@AUC)
#' fakeAUC@AUC[,1:5]
#' auc(fakeAUC)[,1:5]
#' cellNames(fakeAUC)
#' geneSetNames(fakeAUC)
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

setGeneric(name="auc", def=function(object) standardGeneric("auc"))
setMethod("auc",
          signature="cellsAUC",
          definition = function(object) {
            object@AUC
          }
)

setGeneric(name="geneSetNames", def=function(object) standardGeneric("geneSetNames"))
setMethod("geneSetNames",
          signature="cellsAUC",
          definition = function(object) {
            rownames(object@AUC)
            }
)

setGeneric(name="cellNames", def=function(object) standardGeneric("cellNames"))
setMethod("cellNames",
          signature="cellsAUC",
          definition = function(object) {
            colnames(object@AUC)
          }
)
