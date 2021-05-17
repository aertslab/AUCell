# Aux function for AUCell.assignCells
#' @importFrom stats sd
#' @importFrom stats ks.test
#' @importFrom stats rnorm
#' @importFrom stats qnorm
#' @importFrom stats dnorm
#' @importFrom stats density
#' @importFrom utils installed.packages
#' @importFrom utils capture.output
#' @importFrom graphics lines
#' @importFrom graphics points
#' @importFrom graphics hist
#' @importFrom graphics rect
#' @importFrom graphics abline
#' @importFrom graphics text

#### This version:
####    - AUC calculation, same as always,
####    - Threshold: Testing diferent thresholds
####       (3 versions are calculated, and the maximum is chosen automatically)
####    TO DO: decide what to do with global distr (i.e. for HK...)

# V6: smallestPopPercent renamed
# auc <- cells_AUC[,5] # grep("regulation of stress", colnames(AUCell.auc))
.auc_assignmnetThreshold_v6 <- function(aucRow, plotHist=TRUE,
          smallestPopPercent=.25, densAdjust=2, thrP=0.01, nBreaks=100)
{
  gSetName <- rownames(aucRow)[1]
  auc <- aucRow[gSetName,]


  nCells <- length(auc)
  skipGlobal <- TRUE
  skipRed <- FALSE
  skipSmallDens <- FALSE
  commentMsg <- ""
  aucThrs <- c()

  notPopPercent <- 1 - smallestPopPercent
  if(sum(auc==0) > (nCells*notPopPercent))
  {
    skipGlobal <- FALSE
    commentMsg <- paste(commentMsg,
                 round((sum(auc==0)/nCells)*100),
                 "% (more than ", notPopPercent,"%) of AUC are zero. ", sep="")
  }

  meanAUC <- mean(auc)
  sdAUC <- sd(auc)
  maybeNormalDistr <- !suppressWarnings(
    ks.test(auc, rnorm(max(100,length(auc)),mean=meanAUC, sd = sdAUC),
            alternative = "less")$p.value < .01)
  if(maybeNormalDistr){
    commentMsg <- paste0(commentMsg,
            "The AUC might follow a normal distribution (random gene-set?). ")
    skipGlobal <- FALSE

    # aucThrs["outlierOfGlobal"] <- meanAUC + 2*sdAUC
    aucThrs["outlierOfGlobal"] <- qnorm(1-thrP, mean=meanAUC, sd=sdAUC)
  }

  #V6
  histogram <- hist(c(0, auc/max(auc)), breaks=100, plot=FALSE)$count
  if((sum(histogram[1:5]) / sum(histogram)) >= notPopPercent*.75) {
    skipGlobal <- FALSE
    skipRed <- TRUE
    skipSmallDens <- TRUE
  }
  if((sum(histogram[1:10]) / sum(histogram)) >= notPopPercent*.50) {
    skipSmallDens <- TRUE
    skipGlobal <- FALSE
    # skipRed <- TRUE ?
    aucThrs["tenPercentOfMax"] <- max(auc)*.10
  }
  # print(skipRed)

  densCurve <- density(auc, adjust=densAdjust, cut=0)
  maximumsDens <- NULL
  inflPoints <- diff(sign(diff(densCurve$y)))
  maximumsDens <- which(inflPoints==-2)
  globalMax <- maximumsDens[which.max(densCurve$y[maximumsDens])]
  minimumDens <- which(inflPoints==2)
  smallMin <- NULL
  if(!skipSmallDens)
    smallMin <- data.table::last(minimumDens[which(minimumDens < globalMax)]) #1prev to max
  minimumDens <- c(smallMin,
        minimumDens[which(minimumDens > globalMax)]) # all after maximum

  # Density-based threshold (V4):
  # First minimum after the biggest maximum   (adjust=2)
  densTrh <- NULL
  if(length(minimumDens)>0) # && (!skipMinimumDens))
  {
    densTrh <- densCurve$x[min(minimumDens)]
    # Commented on V6
    # Only keep if it is a real inflextion point
    # (i.e. next max at least 5% of the global max)
    if(length(maximumsDens)>0)
    {
      nextMaxs <- maximumsDens[which(densCurve$x[maximumsDens] > densTrh)]
      if((max(densCurve$y[nextMaxs])/max(densCurve$y))<.05)
      {
        densTrh <- NULL
        # print(gSetName)
      }
    }
  }

  ## TO DO: Check special cases with many zeroes
  auc <- sort(auc)
  distrs <- list()
  distrs[["Global_k1"]] <- list(mu=c(meanAUC, NA), sigma=c(sdAUC, NA), x=auc)


  if("mixtools" %in% rownames(installed.packages()))
  {
    na <- capture.output(distrs[["k2"]] <-
        tryCatch(mixtools::normalmixEM(auc, fast=FALSE, k=2, verb=FALSE),
        # With fast, if there are many zeroes, it fails quite often
        error = function(e) {
          return(NULL)
        }))

    na <- capture.output(distrs[["k3"]] <-
        tryCatch(mixtools::normalmixEM(auc, fast=FALSE, k=3, verb=FALSE),
        error = function(e) {
          return(NULL)
        }))

    if(is.null(distrs[["k2"]]) && is.null(distrs[["k3"]]))
    {
      if(sum(auc==0)<(nCells*notPopPercent*.5))
        skipGlobal <- FALSE    # only if not too many zeroes??
      # qpois(1-(thrP/nCells), 1, log = FALSE)
      # plot(sort(rpois(auc, lambda=var(auc)), decreasing=TRUE))

      # commented V6 why was it here??
      # qPop <- quantile(auc, 1-smallestPopPercent)
      # if(sum(auc<qPop) >0)
      #   distrs[["k2"]] <- list(mu=c(mean(auc[auc<qPop]), NA),
      #                          sigma=c(sd(auc[auc<qPop]), NA),
      #                          lambda=c(1,NA), x=auc)
    }
    # if(!skipGlobal) print(gSetName) Warning?

    if(!is.null(distrs[["k2"]]))
    {
      compL <- which.min(distrs[["k2"]][["mu"]])
      compR <- which.max(distrs[["k2"]][["mu"]])
      ### Check distributions
      # Second distribution is "taller" than first one
      height1 <- .4/distrs[["k2"]][["sigma"]][compL]*
        distrs[["k2"]][["lambda"]][compL]
      height2 <- .4/distrs[["k2"]][["sigma"]][compR]*
        distrs[["k2"]][["lambda"]][compR]
      taller <- height1 < height2
      # Use global distr:
      # Mean of the global distr is included within the SD of the first
      # & Both means are included within the mean+SD of the Global distribution
      globalInclInFirst <-
        (distrs[["Global_k1"]]$mu[1] <
        (distrs[["k2"]][["mu"]][compL]+(1.5*distrs[["k2"]][["sigma"]][compL])))
      includedInGlobal <-
        ((distrs[["k2"]][["mu"]][compL] >
        (distrs[["Global_k1"]]$mu[1]-distrs[["Global_k1"]]$sigma[1])) &&
          (distrs[["k2"]][["mu"]][compR] <
          (distrs[["Global_k1"]]$mu[1]+distrs[["Global_k1"]]$sigma[1])))
      if(taller || (globalInclInFirst && includedInGlobal))
      {
        skipGlobal <- FALSE

        if(globalInclInFirst && includedInGlobal)
          commentMsg <- paste(commentMsg,
                "The global distribution overlaps the partial distributions. ")
        if(taller && !includedInGlobal)
          commentMsg <- paste(commentMsg, "The right distribution is taller. ")
      }
    }
  }else{
    warning("Package 'mixtools' is not available to calculate the sub-distributions.")
  }

  glProb <- 1-(thrP/nCells + smallestPopPercent)   ## CORRECT?!?!
  aucThrs["Global_k1"] <- qnorm(glProb,# qnorm(1-(thrP/nCells),
                                mean=distrs[["Global_k1"]][["mu"]][1],
                                sd=distrs[["Global_k1"]][["sigma"]][1])
  if(!is.null(distrs[["k2"]]))
  {
    k2_L <- which.min(distrs[["k2"]][["mu"]]) # (sometimes the indexes are shifted)
    aucThrs["L_k2"] <- qnorm(1-(thrP/nCells),
                             mean=distrs[["k2"]][["mu"]][k2_L],
                             sd=distrs[["k2"]][["sigma"]][k2_L])
  }

  if(!is.null(distrs[["k3"]]))
  {
    k3_R <- which.max(distrs[["k3"]][["mu"]]) # R: right distribution
    k3_R_threshold <- qnorm(thrP,
                            mean=distrs[["k3"]][["mu"]][k3_R],
                            sd=distrs[["k3"]][["sigma"]][k3_R])
    if(k3_R_threshold > 0) aucThrs["R_k3"] <- k3_R_threshold
  }

  if(!is.null(densTrh))
  {
    aucThrs["minimumDens"] <- densTrh
  }

  aucThr <- aucThrs
  if(skipGlobal)
    aucThr <- aucThrs[which(!names(aucThrs) %in% "Global_k1")]
    # TO DO: Decide when to merge with GLOBAL

  if(skipRed)
    aucThr <- aucThrs[which(!names(aucThrs) %in% "L_k2")]
    # TO DO: Decide when to merge with GLOBAL

  aucThr <- aucThr[which.max(aucThr)] # to keep name
  if((length(aucThr)>0) && (names(aucThr) == "minimumDens"))
  {
    maximumsDens <- maximumsDens[which(densCurve$y[maximumsDens]>1)]
    if(length(maximumsDens) > 2)
    {
      tmp <- cbind(minimumDens[seq_len(length(maximumsDens)-1)],
                   maximumsDens[-1])
      FCs <- densCurve$y[tmp[,2]]/densCurve$y[tmp[,1]]
      if(any(FCs > 1.5))
        warning(gSetName,
          ":\tCheck the AUC histogram. ",
          "'minimumDens' was selected as the best threshold, ",
          "but there might be several distributions in the AUC.")
    }
  }

  if("minimumDens" %in% names(aucThrs))
    aucThr <- aucThrs["minimumDens"]
  if(length(aucThr)==0)
    aucThr <- aucThrs[which.max(aucThrs)]
  if(length(aucThr)==0) #should not happen
    aucThr <- 1
  if(length(aucThr)>1) #should not happen
    aucThr <- unlist(aucThr[which.max(aucThr)])
    

  if(plotHist)
  {
    histInfo <- AUCell_plotHist(aucRow,
                         aucThr=aucThr,
                         nBreaks=nBreaks)
    histMax <- max(histInfo[[gSetName]]$counts)

    # Plot density
    densCurve$y <- densCurve$y*(histMax/max(densCurve$y))
    thisLwd <- ifelse(
      (aucThrs["minimumDens"]==aucThr) &&
        (!is.null(aucThr) && !is.null(aucThrs["minimumDens"])),
      3,
      1)
    lines(densCurve, lty=1, lwd=thisLwd, col="blue")
    if(!is.null(minimumDens))
      points(densCurve$x[minimumDens], densCurve$y[minimumDens],
             pch=16, col="darkblue")

    ### Plot distributions
    scalFact <- 1
    # if(!skipGlobal)
    # {
    aucDistr <- dnorm(distrs[["Global_k1"]][["x"]],
                      mean=distrs[["Global_k1"]][["mu"]][1],
                      sd=distrs[["Global_k1"]][["sigma"]][1])
    scalFact <- (histMax/max(aucDistr))*.95

    thisLwd <- ifelse(aucThrs["Global_k1"]==aucThr, 3, 1)
    lines(distrs[["Global_k1"]][["x"]],
          scalFact * aucDistr,
          col="darkgrey", lwd=thisLwd, lty=2)

    if(!is.null(distrs[["k2"]]))
    {
      aucDistr <- dnorm(distrs[["k2"]][["x"]],
                        mean=distrs[["k2"]][["mu"]][k2_L],
                        sd=distrs[["k2"]][["sigma"]][k2_L])
      scalFact <- (histMax/max(aucDistr))*.95


      thisLwd <- ifelse(aucThrs["k2"]==aucThr, 3, 1)
      lines(distrs[["k2"]][["x"]],
            scalFact * aucDistr,
            col="red", lwd=thisLwd, lty=2)

      rect(distrs[["k2"]][["mu"]][k2_L]-distrs[["k2"]][["sigma"]][k2_L],
           histMax-(histMax*.02),
           distrs[["k2"]][["mu"]][k2_L]+distrs[["k2"]][["sigma"]][k2_L],
           histMax, col="#70000030", border="#00009000")
    }

    # print(aucThrs)
    if((!is.null(distrs[["k3"]])) && ("R_k3" %in% names(aucThrs)))
    {
      k3_L <- which.min(distrs[["k3"]][["mu"]]) # (index position not constant)

      aucDistr2 <- dnorm(distrs[["k3"]][["x"]],
                         mean=distrs[["k3"]][["mu"]][k3_R],
                         sd=distrs[["k3"]][["sigma"]][k3_R])
      scalFact2 <- scalFact *
        (distrs[["k3"]][["lambda"]][k3_R]/distrs[["k3"]][["lambda"]][k3_L])

      thisLwd <- ifelse(aucThrs["k3"]==aucThr, 3, 1)
      lines(distrs[["k3"]][["x"]],
            scalFact2*aucDistr2,
            col="magenta", lwd=thisLwd, lty=2)

      rect(distrs[["k3"]][["mu"]][k3_R]-distrs[["k3"]][["sigma"]][k3_R],
           histMax-(histMax*.02),
           distrs[["k3"]][["mu"]][k3_R]+distrs[["k3"]][["sigma"]][k3_R],
           histMax, col="#80808030", border="#80808030")
    }

    ## Add threshold lines
    aucThrs <- aucThrs[!is.na(aucThrs)]
    if(length(aucThrs)>0)
    {
      pars <- list()
      pars[["Global_k1"]] <- c(col1="#909090", col2="black", pos=.9)
      pars[["L_k2"]] <- c(col1="red", col2="darkred", pos=.8)
      # pars[["Max"]] <- c(col1="grey", col2="black", pos=.4)
      pars[["R_k3"]] <- c(col1="magenta", col2="magenta", pos=.6)
      pars[["minimumDens"]] <- c(col1="blue", col2="darkblue", pos=.4)
      pars[["tenPercentOfMax"]] <- c(col1="darkgreen", col2="darkgreen", pos=.9)
      pars[["outlierOfGlobal"]] <- c(col1="darkgreen", col2="darkgreen", pos=.9)

      for(thr in names(aucThrs))
      {
        thisLwd <- ifelse(aucThrs[thr]==aucThr, 5, 2)
        thisLty <- ifelse(aucThrs[thr]==aucThr, 1, 3)

        abline(v=aucThrs[thr], col=pars[[thr]][1], lwd=thisLwd, lty=thisLty)
        xPos <- aucThrs[thr]*1.01
        if(aucThrs[thr] > (max(auc)*.8))
          xPos <- 0
        if(aucThrs[thr]==aucThr)
          text(xPos, histMax*as.numeric(pars[[thr]][3]),
               pos=4, col=pars[[thr]][2], cex=.8,
               paste("AUC > ", signif(aucThrs[thr],2),
                     "\n(",sum(auc>aucThrs[thr])," cells)", sep=""))
      }
    }
  }
  return(list(selected=aucThr,
       thresholds=cbind(threshold=aucThrs,
                        nCells=sapply(aucThrs, function(x) sum(auc>x))),
       comment=commentMsg))
}
