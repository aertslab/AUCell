# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @importFrom stats setNames
#' @importFrom methods new
#'
#' @title Update AUCell results
#' @description Updates the AUC scores provided by AUCell from a previous version.
#' @param oldAucObject Object to update
#' @param objectType Either "AUC" or "ranking" indicating the object type
#' @return Updated version of the object as \code{\link{aucellResults}}.
#' @examples
#' oldAuc <- matrix(data=1:2000, nrow=50, ncol=40)
#' updateAucellResults(oldAuc)
#' @export
#'
updateAucellResults <- function(oldAucObject, objectType="AUC")
{
  objectType <- tolower(objectType)
  if(!objectType %in% c("auc", "ranking")) stop("'objectType' should be either 'AUC' or 'ranking'.")

  newAucObject <- NULL
  if(methods::is(oldAucObject,"matrixWrapper")) 
  {
    aucMatrix <- oldAucObject@matrix
    #library(SummarizedExperiment)
    if(objectType== "auc") newAucObject <- new("aucellResults", SummarizedExperiment::SummarizedExperiment(assays=list(AUC=aucMatrix)))
    if(objectType== "ranking") newAucObject <- new("aucellResults", SummarizedExperiment::SummarizedExperiment(assays=list(ranking=aucMatrix)))
  }

  if(methods::is(oldAucObject,"matrix"))
  {
    aucMatrix <- oldAucObject
    # library(SummarizedExperiment)
    if(objectType== "auc") newAucObject <- new("aucellResults", SummarizedExperiment::SummarizedExperiment(assays=list(AUC=aucMatrix)))
    if(objectType== "ranking") newAucObject <- new("aucellResults", SummarizedExperiment::SummarizedExperiment(assays=list(ranking=aucMatrix)))
  }

  return(newAucObject)
}
