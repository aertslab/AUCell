# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @title Plot AUC histogram
#' @description Plots the distribution of AUC across the cells
#' (for each gene-set) as an histogram.
#' @param cellsAUC Subset of the object returned by \code{\link{AUCell_calcAUC}}
#' (i.e. including only the gene-sets to plot)
#' @param aucThr AUC value planned to use as threshold
#' (to make sure the X axis includes it), if any.
#' Otherwise, the X axis extends to cover only the AUC values plotted.
#' @param nBreaks Number of 'bars' to plot (breaks argument for hist function).
#' @param onColor Color for the bars that pass the AUC threshold
#' @param offColor Color for the bars that do not pass the AUC threshold
#' @param ... Other arguments to pass to \code{\link{hist}} function.
#' @return List of histogram objects (invisible).
#' @seealso See the package vignette for examples and more details:
#' \code{vignette("AUCell")}
#' @example inst/examples/example_AUCell_plotHist.R
#' @export
AUCell_plotHist <- function(cellsAUC, aucThr=max(cellsAUC), nBreaks=100, 
                            onColor="dodgerblue4", offColor="slategray2", ...)
{
  if(methods::is(cellsAUC,"aucellResults")) {
    cellsAUC <- getAUC(cellsAUC)
  }
  if(!is.matrix(cellsAUC))
    stop("The first argument should be the output from 'AUCell_calcAUC()'")

  ret <- lapply(seq_len(nrow(cellsAUC)), function(gsn) {
    .auc_plotHist(auc=cellsAUC[gsn,], gSetName=rownames(cellsAUC)[gsn],
              aucThr=aucThr, nBreaks=nBreaks, 
              onColor=onColor, offColor=offColor, ...)
    })
  names(ret) <- rownames(cellsAUC)

  invisible(ret)
}
# @export
#AUCell_plot <- AUCell_plotHist # Alias (backwards compatibility) TO DO: Replace by conditional?

# only for one gene-set
.auc_plotHist <- function(auc, gSetName, aucThr, 
                        nBreaks=100, 
                        onColor="dodgerblue4", offColor="slategray2", 
                        ...)
{
  maxLim <- max(c(auc, aucThr), na.rm=T)
  # col <- rep("lightskyblue2", nBreaks)
  # col[passThr] <- "dodgerblue2"
  col <- rep(offColor, nBreaks)
  borderCol <- rep("lightsteelblue1", nBreaks)
  passThr <- hist(auc, breaks=nBreaks, plot=FALSE)$breaks >= aucThr
  col[passThr] <- onColor
  borderCol[passThr] <- "cornflowerblue"
  
  ret <- hist(auc, breaks=nBreaks, xlim=c(0,maxLim),
    col=col, border=borderCol,
    main=gSetName, xlab="AUC histogram", ...)
  return(ret)
}

# Default plots for binary and AUC are available in priv_plots.R
