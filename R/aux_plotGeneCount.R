# Used by buildRankings
# Input: matrix or data.table
# If data.table, do not include the "rn" column (rownames)


# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
#' @import data.table
#'
#' @title to do
#' @description to do
#' @param to do
#' @return to do
#' @examples  # to do
#' @export
plotGeneCount <- function(exprMat, verbose=TRUE)
{
  if((is.data.table(exprMat)) && (colnames(exprMat)[1] == "rn")) stop('The data.table contains a column with rownames (to skip, i.e. dt[,-"rn", with=FALSE]')
  countByCell <- colSums(exprMat>0)
  sampleType <- "cell"
  plot.new()
  par(fig=c(0,1,0,.8), new=TRUE)
  na <- hist(countByCell, main="", col="skyblue", xlab="# of genes", ylab="# of cells", cex.axis=.75, xlim=c(0,max(countByCell)))
  par(fig=c(0,1,0.55,1), new=TRUE)
  na <- boxplot(countByCell, ylim=c(0,max(countByCell)), range=0,
                col="skyblue", cex.axis=.75, horizontal=TRUE, axes=FALSE)
  mtext("# of genes detected by cell", side=3, outer=TRUE, line=-3, cex=2)

  if(verbose) message("Quantiles for the number of genes detected by cell: \n(Non-detected genes are shuffled at the end of the ranking. Keep in mind when choosing the threshold for calculating the AUC).")
  ret <- c(min=min(countByCell),quantile(countByCell, c(.01,.05, .10, .50, 1)))
  print(ret)
  invisible(ret)
}
