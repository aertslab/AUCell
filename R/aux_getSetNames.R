#' @export
getSetNames <- function(aucMat, patterns, startChr="^", endChr=" |_")
{
  setNames(sapply(patterns, function(x) grep(paste0(startChr,x,endChr), rownames(aucMat), value=T)),
           names(patterns))
}
