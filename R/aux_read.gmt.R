# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)
#' @title to do
#' @description to do
#' @param to do
#' @return to do
#' @examples #to do
#' @export
read.gmt <- function(fileName) {
  tmp <- sapply(readLines(fileName), function(x) strsplit(x,"\t"))
  names(tmp) <- sapply(tmp, function(x) x[1])
  lapply(tmp, function(x) x[3:length(x)])
}
