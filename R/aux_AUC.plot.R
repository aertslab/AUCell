# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @title Plot AUC histogram
#' @description Plots the distribution of AUC across the cells (for a given gene-set) as an histogram.
#' @param auc A row from the AUC matrix returned by \code{\link{AUCell.calcAUC}}
#' @param gSetName Title for the plot (recommended: Gene set name or description)
#' @param aucThr AUC value planned to use as threshold (to make sure the X axis includes it), if any. Otherwise, the X axis extends to cover all the AUC values plotted.
#' @param nBreaks Number of 'bars' to plot (breaks argument for hist function).
#' @param returnInfo Whether to return the histogram object (TRUE/FALSE).
#' @param ... Other arguments to pass to \code{\link{hist}} function.
#' @return Histogram object (if requested with returnInfo=TRUE).
#' @seealso See the package vignette for examples and more details: \code{vignette("AUCell")}
#' @example inst/examples/example_AUCplot.R
#' @export
AUC.plot <- function(auc, gSetName="AUC distribution for this geneset", aucThr=max(auc), nBreaks=100, returnInfo=FALSE, ...)
{
  ret <- hist(auc, breaks=nBreaks, col="#6666aa80", border="#5588bb", main=gSetName, xlab="AUC histogram", xlim=c(0,max(c(auc, aucThr))), ...)
  if(returnInfo) return(ret)
}
