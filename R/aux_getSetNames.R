
#' @title Get gene set name (excluding number of genes/regions after space)
#' @description Returns the gene set name (i.e. selects the given pattern)
#' @param aucMat AUC matrix
#' @param patterns patterns
#' @param startChr  Character to indicate the start (typically "^")
#' @param endChr Character at the end of the gene set name, i.e. space or "_" after a transcription factor name
#' @return Returns the gene set name (i.e. selects the given pattern)

#' @export
getSetNames <- function(aucMat, patterns, startChr="^", endChr=" |_")
{
  setNames(sapply(patterns, function(x) grep(paste0(startChr,x,endChr), rownames(aucMat), value=T)),
           names(patterns))
}
