# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @title Plot AUC histogram
#' @description Plots the distribution of AUC across the cells (for each gene-set) as an histogram.
#' @param cellsAUC Subset of the object returned by \code{\link{AUCell_calcAUC}} (i.e. including only the gene-sets to plot)
#' @param aucThr AUC value planned to use as threshold (to make sure the X axis includes it), if any. Otherwise, the X axis extends to cover only the AUC values plotted.
#' @param nBreaks Number of 'bars' to plot (breaks argument for hist function).
#' @param ... Other arguments to pass to \code{\link{hist}} function.
#' @return List of histogram objects (invisible).
#' @seealso See the package vignette for examples and more details: \code{vignette("AUCell")}
#' @example inst/examples/example_AUCell_plot.R
#' @export
AUCell_plot <- function(cellsAUC, aucThr=max(cellsAUC), nBreaks=100, ...)
{
  if(class(cellsAUC)[1]=="aucellResults") {
    cellsAUC <- getAUC(cellsAUC)
  }
  if(!is.matrix(cellsAUC)) stop("The first argument should be the output from 'AUCell_calcAUC()'")

  ret <- lapply(seq_len(nrow(cellsAUC)), function(gsn) {
    .auc_plot(auc=cellsAUC[gsn,], gSetName=rownames(cellsAUC)[gsn], aucThr=aucThr, nBreaks=nBreaks, ...)
    })
  names(ret) <- rownames(cellsAUC)

  invisible(ret)
}

# only for one gene-set
.auc_plot <- function(auc, gSetName, aucThr, nBreaks, ...)
{
  ret <- hist(auc, breaks=nBreaks, col="#6666aa80", border="#5588bb",
          main=gSetName, xlab="AUC histogram", xlim=c(0,max(c(auc, aucThr))), ...)
  return(ret)
}
