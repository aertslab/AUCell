# txt <- paste0("Cells with AUC > ", signif(gsTrheshold, 2))
.auc_plotBinaryTsne <- function(tSNE, selectedCells,  
                                title="", txt="Coloured cells pass the threshold",
                                cex=1,
                                alphaOn=1, alphaOff=0.2,
                                onColor="dodgerblue4", offColor="lightgray", 
                                borderColor=adjustcolor("lightgray", alpha=.1),
                                ...)  
{
  pointBg <-  setNames(rep(adjustcolor(offColor, alpha=alphaOff), nrow(tSNE)), rownames(tSNE))
  pointBg[selectedCells] <- adjustcolor(onColor , alpha=alphaOn) # "royalblue4"
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
                                  borderColor=adjustcolor("darkgrey", alpha=.3),
                                  cex=1,
                                  ...)
{
  # Setup palette
  nBreaks <- 5
  colorPal <- grDevices::colorRampPalette(colorsForPal)(nBreaks)
  colorPal <- adjustcolor(colorPal, alpha=alphaOn)
  
  # Assign cell color
  cellColor <- setNames(rep("black", nrow(tSNE)), rownames(tSNE))
  pointBorder <- setNames(rep(borderColor, nrow(tSNE)), rownames(tSNE))
  if (sum(cellProp) > 0)
  {
    # On
    cellColor[names(cellProp)] <- colorPal[as.numeric(cut(cellProp, breaks=nBreaks, right=F,include.lowest=T))]
    # Off
    lowLim <- as.numeric(gsub("]", "", strsplit(levels(cut(cellProp, breaks=100))[1], ",")[[1]][2]))
    cellColor[which(cellProp <= max(0, lowLim))] <- adjustcolor(offColor, alpha=alphaOff)
  } else
  {
    cellColor[names(cellProp)] <- c(sadjustcolor(offColor, alpha=alphaOff), length(cellProp))
  }
  
  # Plot
  pointBg <- cellColor
  plot(tSNE, pch=21, cex=cex,
       col=pointBorder, bg=pointBg,
       main=title, xlab=txt,
       axes=FALSE, ylab="", ...)
  box(which = "plot", col="grey")
}
