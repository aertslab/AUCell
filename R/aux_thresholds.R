#' @rdname AUCell_exploreThresholds
#' @aliases getThresholdSelected
#' @export
getThresholdSelected <- function(aucellThresholds){
  sapply(aucellThresholds, function(x) {
    if(length(x) > 1) return(unname(x$aucThr$selected))
    x
  })
}

#' @rdname AUCell_exploreThresholds
#' @aliases getThresholdSelected
#' @export
getAssignments <- function(aucellThresholds){
  lapply(aucellThresholds, function(x) x$assignment)
}
