# .auc_plotGradientTsne
# .auc_plotBinaryTsne
# .cellProps_plotTsne  # For APP, only one prop
#     plotTsne_cellProps: For SCENIC, iterator: plots all properties

########################################################################
# txt <- paste0("Cells with AUC > ", signif(gsTrheshold, 2))
.auc_plotBinaryTsne <- function(tSNE, selectedCells,  
                                title="", txt="Coloured cells pass the threshold",
                                cex=1,
                                alphaOn=1, alphaOff=0.2,
                                onColor="dodgerblue4", offColor="lightgray", 
                                borderColor=adjustcolor("lightgray", alpha.f=.1),
                                ...)  
{
  pointBg <-  setNames(rep(adjustcolor(offColor, alpha.f=alphaOff), nrow(tSNE)), rownames(tSNE))
  pointBg[selectedCells] <- adjustcolor(onColor , alpha.f=alphaOn) # "royalblue4"
  pointBorder <- setNames(rep(borderColor, nrow(tSNE)), rownames(tSNE))
  
  ### Plot
  plot(tSNE, pch=21, cex=cex,
       col=pointBorder, bg=pointBg,
       main=title, xlab=txt,
       axes=FALSE, ylab="", ...)
  box(which = "plot", col="grey")
}

.auc_plotGradientTsne <- function(tSNE, cellProp,
                                  colorsForPal=c("mistyrose", "red"),
                                  title="", txt="Regulon activity (AUC)", 
                                  alphaOn=1, alphaOff=0.2,
                                  offColor="lightgray", 
                                  borderColor=adjustcolor("darkgrey", alpha.f=.3),
                                  cex=1,
                                  ...)
{
  # Setup palette
  nBreaks <- 5
  colorPal <- grDevices::colorRampPalette(colorsForPal)(nBreaks)
  colorPal <- adjustcolor(colorPal, alpha.f=alphaOn)
  
  # Assign cell color
  cellColor <- setNames(rep("black", nrow(tSNE)), rownames(tSNE))
  pointBorder <- setNames(rep(borderColor, nrow(tSNE)), rownames(tSNE))
  if (sum(cellProp) > 0)
  {
    # On
    cellColor[names(cellProp)] <- colorPal[as.numeric(cut(cellProp, breaks=nBreaks, right=FALSE,include.lowest=TRUE))]
    # Off
    lowLim <- 0
    tryCatch({
      lowLim <- as.numeric(gsub("]", "", strsplit(levels(cut(cellProp, breaks=100))[1], ",")[[1]][2]))
    },error = function(e) {
      save(cellProp, file="error_cellProp.RData")
      message("There was an error trying to plot the t-SNE. Please report it in https://github.com/aertslab/AUCell (file for debugging: error_cellProp.RData). ",
              "\nError message: ", e$message)
    })
  } else
  {
    cellColor[names(cellProp)] <- c(adjustcolor(offColor, alpha.f=alphaOff), length(cellProp))
  }
  
  # Plot
  pointBg <- cellColor
  plot(tSNE, pch=21, cex=cex,
       col=pointBorder, bg=pointBg,
       main=title, xlab=txt,
       axes=FALSE, ylab="", ...)
  box(which = "plot", col="grey")
}


########################################################################
####                    CELL properties                            #####
#' @title plotTsne_cellProps
#' @description Plots the t-SNE coloured based on the known the cell properties
#' @param tSNE t-SNE coordinates (e.g. \code{tSNE$Y})
#' @param cellInfo Dataframe with cell phenodata
#' @param colVars Colors for the cell properties (optional)
#' @param cex Scaling factor for the dots in the scatterplot
#' @param sub Subtitle (e.g. tSNE type)
#' @param gradientCols Gradient colors for numerical variables
#' @param showLegend Whether to show the legend
#' @return Plots the t-SNE
#' @example inst/examples/example_AUCell_plotTSNE.R
#' @export
plotTsne_cellProps <- function(tSNE, cellInfo, colVars=NULL,
                               cex=1, sub="", 
                               gradientCols=c("yellow", "orange","red"),
                               showLegend=TRUE)
{
  if(!is.matrix(tSNE)) stop("tSNE should be a matrix.")
  
  # If it is null, plot all cells black
  if(is.null(cellInfo)) {
    cellInfo <- data.frame(row.names=rownames(tSNE), tSNE=rep("xx",nrow(tSNE)))
    showLegend <- FALSE
  }
  
  cellInfo <- data.frame(cellInfo)[rownames(tSNE),,drop=F]
  for(varName in colnames(cellInfo))
  {
    .cellProps_plotTsne(tSNE, cellInfo, varName, colVars=colVars, cex=cex,
                        sub=sub, gradientCols=gradientCols, showLegend=showLegend)
  }
}


.cellProps_plotTsne <- function(tSNE, cellInfo, varName=NULL, colVars=NULL, cex=1,
                                sub="", gradientCols=c("yellow", "orange","red"),
                                showLegend=TRUE)
{
  if(!is.null(varName))
  {
    if(!methods::is(cellInfo[,varName],"numeric"))
    {
      if(is.null(colVars[[varName]])) {
        varLevels <- as.character(unique(cellInfo[,varName]))
        colVars[varName] <- list(setNames(rep("black", length(varLevels)), varLevels))
        if(length(colVars[[varName]])>1) colVars[[varName]] <- setNames(rainbow(length(unique(cellInfo[,varName]))), unique(cellInfo[,varName]))
      } 
      cellColor <- setNames(rep("#30303010", nrow(tSNE)), rownames(tSNE))
      thisProp <- setNames(as.character(cellInfo[, varName]), rownames(cellInfo))
      thisProp <- thisProp[which(!is.na(thisProp))]
      propLevels <- levels(factor(thisProp))
      cellColor[names(thisProp)] <- colVars[[varName]][thisProp]
    }else
    {
      # if(is.null(colVars[[varName]])) {
      #   colVars[[varName]] <- setNames(rainbow(length(unique(cellInfo[,varName]))), unique(cellInfo[,varName]))
      # } 
      colorPal <- grDevices::colorRampPalette(gradientCols)
      cellColor <- setNames(adjustcolor(colorPal(10), alpha.f=.8)[as.numeric(cut(cellInfo[,varName],breaks=10, right=FALSE,include.lowest=TRUE))], rownames(cellInfo))
    }
    
    colsLegend <- colVars[[varName]]
    if(methods::is(cellInfo[,varName],"numeric"))
      colsLegend <- setNames(gradientCols[c(length(gradientCols),1)], signif(c(max(cellInfo[,varName]), min(cellInfo[,varName])),2))
    
    plot(tSNE, pch=16, cex=cex,
         col=cellColor[rownames(tSNE)],
         main=varName, 
         sub=sub,
         axes=FALSE, xlab="", ylab="")
    if(showLegend) legend(min(tSNE[,1]), max(tSNE[,2]), names(colsLegend), col=colsLegend,
                          bg=NULL,border=NULL, box.lwd="none", bty = "n", cex=.8, pch=16)
    box(which = "plot", col="grey")
  }
}
