# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
#' @title to do
#' @description to do
#' @param to do
#' @return to do
#' @examples #to do
#' @export
AUC.plot <- function(auc, gSetName="AUC distribution for this geneset", aucThr=max(auc), returnInfo=FALSE)
{
  ret <- hist(auc, breaks=100, col="#6666aa80", border="#5588bb", main=gSetName, xlab="AUC histogram", xlim=c(0,max(c(auc, aucThr))))
  if(returnInfo) return(ret)
}
