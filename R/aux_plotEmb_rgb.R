#' @title plotEmb_rgb
#' @description Colors the embeddings (t-SNE/Umap) based on the activity of 3 (groups of) geneSets
#' @param aucMat AUC matrix (continuous or binary)
#' @param embedding AUC matrix (continuous or binary)
#' @param geneSetsByCol Gene sets to plot
#' @param aucType "AUC" or "Binary"
#' @param aucMaxContrast To increase the AUC contrast decrease the value.
#' @param offColor Color por the cells completelly off. To deactivate (color as black), set to NULL.
#' @param showPlot Whether to plot the coloured embeddings.
#' @param showLegend Whether to plot add a legend to the plot.
#' @param ... Other arguments to pass to the \code{plot} function.
#' @return The cell colors (invisible)
#' @examples 
# par(mfrow=c(1,2))
# 
# setNames <- c("Dlx1", "Sox8")
# cellCol <- plotEmb_rgb(mat2col, emb, setNames, aucType="AUC", aucMaxContrast=0.6)
# text(-5,-23, attr(cellCol,"red"), col="red", cex=.7)
# text(-10,-18, attr(cellCol,"green"), col="green", cex=.7)
# 
# setNames <- list(red=c("Dlx1","Dlx5"),
#                      green=c("Neurod1"),
#                      blue=c("Sox8"))
# cellCol <- plotEmb_rgb(mat2col, emb, setNames, aucType="Binary")


#' @export
plotEmb_rgb <- function(aucMat, embedding, geneSetsByCol, 
                        aucType="AUC", 
                        aucMaxContrast=0.8, offColor="#c0c0c030", 
                        showPlot=TRUE, showLegend=TRUE, ...)
{
  # Check format
  if(length(geneSetsByCol)>3) stop("To plot more than three regulons, group them by color.")
  if(any(!names(geneSetsByCol) %in% c("red","green", "blue"))) 
    stop('If a list, the names of geneSetsByCol should be "red","green", and/or"blue".')
  if(aucMaxContrast<=0 | aucMaxContrast>1) stop("aucMaxContrast should be between 0 and 1")
  if(!tolower(aucType) %in% c("binary","bin", "auc","continuous","cont")) stop("aucType should be binary or auc/continuous")
  
  if(is.null(names(geneSetsByCol)) && length(geneSetsByCol)<=3) {
    names(geneSetsByCol) <- c("red","green", "blue")[seq_along(geneSetsByCol)]
    geneSetsByCol <- as.list(geneSetsByCol)
  }
  
  geneSetsByCol <- geneSetsByCol[lengths(geneSetsByCol)>0]
  
  # Average of binary...
  if(tolower(aucType) %in% c("binary","bin"))
    cellColChan <- sapply(geneSetsByCol, function(modsCol) apply(getAUC(aucMat)[modsCol,, drop=FALSE], 2, mean))
  # Color if all modules are "on"
  # cellColChan <- sapply(geneSetsByCol, function(modsCol) apply(aucMat[modsCol,, drop=FALSE], 2, function(x) as.numeric(sum(x)==length(x))))
  
  # AUC
  if(tolower(aucType) %in% c("auc","continuous","cont")) {
    cellColChan <- sapply(geneSetsByCol, function(regCol) {
      aucAvg <- base::colMeans(getAUC(aucMat[regCol,, drop=FALSE]))
      setNames(sapply(as.numeric(aucAvg/(max(aucAvg)*aucMaxContrast)), min, 1), names(aucAvg))
    })
  }
  
  # Apply color
  missingCol <- setdiff(c("red","green", "blue"), colnames(cellColChan))
  if(length(missingCol)>0)
    cellColChan <- cbind(cellColChan, matrix(rep(0, nrow(cellColChan)*length(missingCol)),ncol=length(missingCol), dimnames=list(NULL,missingCol)))
  cellCol <- apply(cellColChan, 1, function(x) rgb(x["red"], x["green"], x["blue"], alpha=.8))
  if(!is.null(offColor)) cellCol[which(cellCol=="#000000CC")] <- offColor # mostly for binary
  names(cellCol) <- colnames(aucMat)
  
  if(aucType=="binary") attr(cellCol,"Description") <- "Color: average of BINARY regulon activity"
  if(aucType=="auc") attr(cellCol,"Description") <- "Color: average of regulon activity (AUC)"
  attr(cellCol, "red") <- paste0(geneSetsByCol$red, collapse=", ")
  attr(cellCol, "green") <- paste0(geneSetsByCol$green, collapse=", ")
  attr(cellCol, "blue") <- paste0(geneSetsByCol$blue, collapse=", ")
  
  if(showPlot) 
  {
    plot(embedding, col=cellCol[rownames(embedding)], pch=16, 
         sub=attr(cellCol, "Description"), axes=FALSE, ...)
    
    if(showLegend)
    {
      cellColChan[which(cellColChan < (aucMaxContrast/2), arr.ind = T)] <- 0
      cellColChan <- cellColChan[which(apply(cellColChan, 1, function(x) any(x>0))),]
      cellLabels <- setNames(colnames(cellColChan)[apply(cellColChan, 1, function(x) which.max(x))], rownames(cellColChan))
      labsCoords <- t(sapply(split(data.frame(embedding), as.character(cellLabels[rownames(embedding)])), colMeans));
      
      geneSetsByCol[rownames(labsCoords)] <- sapply(geneSetsByCol[rownames(labsCoords)], paste, collapse=", ")
      for(i in rownames(labsCoords)) text(mean(labsCoords[i,1]), mean(labsCoords[i,2]), geneSetsByCol[i], cex=1, col="black")
    }
  }
  invisible(cellCol)
}

