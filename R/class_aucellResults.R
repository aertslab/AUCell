#' @title Wrapper to the matrix that stores the AUC or the cell rankings.
#' @aliases getAUC getRanking show
#' @description
#' This class extends SummarizedExperiment to contain the AUC matrix and cell
#' rankings (as 'assays').
#'
#' The results are stored in the assays slot, but they can be accessed through
#' the regular methods (i.e. nrow, rownames... )
#'
#' Types:
#'
#'  - "AUC": The assays contains the AUC for the gene-sets (or region-sets)
#'  & cells.
#'
#'  - "ranking": The assays contains the gene rankings for each cell.
#'
#' @param object Results from \code{AUCell_buildRanking}
#' or \code{AUCell_calcAUC}.
#' @return
#' \itemize{
#' \item show: Prints a summary of the object
#' \item getAUC: Returns the matrix containing the AUC
#' \item getRanking: Returns the matrix containing the rankings
#' \item cbind: Combines objects by columns (cbind on \code{assays}); other other slots are conserved from the first object provided as argument.
#' Both, ranking and AUC are calculated by column (cell or sample). Therefore, it is fine to merge objects as long as they come from equivalent datasets (and keep same genes/genesets, etc...)
#' }
#' @method show aucellResults
#' @method getAUC aucellResults
#' @method getRanking aucellResults
#' @example inst/examples/example_aucellResults_class.R
#' @rawNamespace import(data.table, except = shift)
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importMethodsFrom SummarizedExperiment assay assays assayNames
#' @importFrom utils head
#' @rdname aucellResults-class
#' @export aucellResults
aucellResults <- setClass("aucellResults",
         contains="SummarizedExperiment",
         representation=representation(
           nGenesDetected = "numeric"
         )
)
#' @importFrom R.utils capitalize
#' @rdname aucellResults-class
#' @aliases show,aucellResults-method
setMethod("show",
  signature="aucellResults",
  definition = function(object) {
    message <- paste(R.utils::capitalize(assayNames(object)), " for ",
           nrow(object)," ", names(dimnames(assay(object)))[1],
           " (rows) and ",
           ncol(object)," ", names(dimnames(assay(object)))[2],
           " (columns).\n", sep="")

    if(assayNames(object) == "AUC")
      message <- paste(message,
                       "\nTop-left corner of the AUC matrix:\n", sep="")
    if(assayNames(object) == "ranking")
      message <- paste(message,
                       "\nTop-left corner of the ranking:\n", sep="")

    cat(message)
    show(head(assay(object)[,seq_len(min(5, ncol(object))),drop=FALSE]))

    if(is.numeric(object@nGenesDetected) &&
       (length(object@nGenesDetected)>0)) {
      cat("\nQuantiles for the number of genes detected by cell:\n")
      show(object@nGenesDetected)
    }
  }
)
##### Access the matrix:
#' @name getAUC
#' @rdname aucellResults-class
# Export generic, for RcisTarget
#' @export getAUC 
setGeneric(name="getAUC",
           def=function(object) standardGeneric("getAUC"))

#' @rdname aucellResults-class
#' @aliases getAUC,aucellResults-method
#' @exportMethod getAUC
setMethod("getAUC",
          signature="aucellResults",
          definition = function(object) {
            if("AUC" %in% assayNames(object)) {
              SummarizedExperiment::assays(object)[["AUC"]]
            }else{
              stop("This object does not contain an AUC matrix.")
            }
          }
)

##### Access the rankings:
# setGeneric
# @method test data.frame
#' @name getRanking
#' @rdname aucellResults-class
# Export generic, for RcisTarget
#' @export getRanking 
setGeneric(name="getRanking",
           def=function(object) standardGeneric("getRanking"))

#' @rdname aucellResults-class
#' @aliases getRanking,aucellResults-method
#' @exportMethod getRanking
setMethod("getRanking",
          signature="aucellResults",
          definition = function(object) {
            if("ranking" %in% assayNames(object)) {
              SummarizedExperiment::assays(object)[["ranking"]]
            }else{
              stop("This object does not contain a ranking.")
            }
          }
)

setGeneric("cbind", signature="...")
##### Combine objects (by colums):
#' @name cbind
#' @rdname aucellResults-class
#' @export
setMethod("cbind", "aucellResults", function(..., deparse.level=1) {
  args <- list(...)
  
  objectType <- unique(sapply(args, function(x) names(assays(x))))
  if(length(objectType)>1) stop("All objects should be of the same type (e.g. ranking OR AUC).")
  dimNames <- apply(sapply(args, function(x) names(dimnames(assay(x)))), 1, function(x) unique(x))
  if(length(dimNames)!=2)  stop("Dimnames do not match.")
  
  allAssays <- lapply(args, assay, withDimnames=TRUE)
  if(length(unique(sapply(allAssays, nrow)))>1)  stop("Number of rows (",dimNames[1],") do not match.")
  if(!all(apply(rbind(sapply(allAssays, rownames)), 1, function(x) length(unique(x)))==1))  stop("Rownames (",dimNames[1],") do not match.")
  if(any(table(unlist(sapply(allAssays, colnames)))>1))  stop("Some column IDs (",dimNames[2],") are duplicated.")
  
  allAssays <- do.call(cbind, allAssays)
  names(dimnames(allAssays)) <- dimNames
  
  old.validity <- S4Vectors:::disableValidity()
  S4Vectors:::disableValidity(TRUE)
  on.exit(S4Vectors:::disableValidity(old.validity))
  
  #### Slots used:
  # new("aucellResults",
  #     SummarizedExperiment::SummarizedExperiment(assays=list(ranking=exprMat)),
  #     nGenesDetected=nGenesDetected)
  # new("aucellResults",
  #     SummarizedExperiment::SummarizedExperiment(assays=list(AUC=aucMatrix)))
  ####
  
  out <- callNextMethod()
  BiocGenerics:::replaceSlots(out, 
                              assays=setNames(list(allAssays),  objectType),
                              check=FALSE)
})
