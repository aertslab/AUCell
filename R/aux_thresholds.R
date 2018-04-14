#' @rdname AUCell_exploreThresholds
#' @aliases getThresholdSelected
#' @export
getThresholdSelected <- function(aucellThresholds){
  sapply(aucellThresholds, function(x) unname(x$aucThr$selected))
}

#' @rdname AUCell_exploreThresholds
#' @aliases getThresholdSelected
#' @export
getAssignments <- function(aucellThresholds){
  lapply(aucellThresholds, function(x) x$assignment)
}
