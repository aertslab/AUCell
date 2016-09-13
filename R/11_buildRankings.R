# Version 2: Using data.table



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
AUCell.buildRankings <- function(exprMat, plotStats=TRUE, nCores=1, verbose=TRUE)
{
  if(!is.data.table(exprMat)) exprMat <- data.table(exprMat, keep.rownames=TRUE)
  setkey(exprMat, "rn") # (reorders rows)

  if(plotStats)
  {
    msg <- tryCatch(plotGeneCount(exprMat[,-"rn", with=FALSE], verbose=verbose),
                                                    error = function(e) {
                                                      return(e)
                                                    })
    if("error" %in% class(msg)) warning(paste("There has been an error in plotGeneCount() [Message: ", msg$message, "]. Proceeding to calculate the rankings...", sep=""))
  }

  colsNam <- colnames(exprMat)[-1] # 1=row names
  if(nCores==1)
  {
    # The rankings are saved in exprMat (i.e. By reference)
    exprMat[, (colsNam):=lapply(-.SD, frank, ties.method="random"), .SDcols=colsNam]

  }else
  {
    suppressMessages(require("doMC", quietly=TRUE)); require("doRNG", quietly=TRUE)
    registerDoMC(nCores)
    if(verbose) message(paste("Using", getDoParWorkers(), "cores."))

    suppressWarnings(colsNamsGroups <- split(colsNam, (1:length(colsNam)) %% nCores)) # Expected warning: Not multiple
    rowNams <- exprMat$rn
    exprMat <- foreach(colsGr=colsNamsGroups, .combine="cbind") %dorng%
    {
      # exprMat[, (colsGr):=lapply(-.SD, frank, ties.method="random"), .SDcols=colsGr]  # Edits by reference: how to make it work in paralell...?
      subMat <- exprMat[,colsGr, with=FALSE]
      subMat[, (colsGr):=lapply(-.SD, frank, ties.method="random"), .SDcols=colsGr]
    }
    exprMat <- data.table(rn=rowNams, exprMat[,colsNam, with=FALSE]) # Keep initial order & recover rownames
    setkey(exprMat, "rn")
  }
  return(exprMat)
}

