#' @rdname AUCell_exploreThresholds
#' @aliases getThresholdSelected
#' @export
getThresholdSelected <- function(aucellThresholds){
  sapply(aucellThresholds, function(x) {
    ret <- x # already only one threshold
    if(length(x) > 1) {
      #otherwise do not modify, it might already be the selected
      if(("aucThr" %in% names(x)) & ("selected" %in% names(x$aucThr)))  
        ret <- unname(x$aucThr$selected)
    }
    return(ret)
  })
}

#' @rdname AUCell_exploreThresholds
#' @aliases getThresholdSelected
#' @export
getAssignments <- function(aucellThresholds){
  lapply(aucellThresholds, function(x) x$assignment)
}
