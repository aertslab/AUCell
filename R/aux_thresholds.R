#' @rdname AUCell_exploreThresholds
#' @aliases getThresholdSelected
#' @export
getThresholdSelected <- function(aucellThresholds){
  sapply(aucellThresholds, function(x) {
    ret <- x # already only one threshold
    if(length(x) > 1) ret <- unname(x$aucThr$selected)
    return(ret)
  })
}

#' @rdname AUCell_exploreThresholds
#' @aliases getThresholdSelected
#' @export
getAssignments <- function(aucellThresholds){
  lapply(aucellThresholds, function(x) x$assignment)
}
