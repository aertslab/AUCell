#' @importFrom DelayedArray DelayedArray blockApply colAutoGrid
blocked_AUCell <- function(exprMat, geneSets, keepZeroesAsNA=FALSE, normAUC=TRUE, aucMaxRank=ceiling(0.05*nrow(exprMat)), BPPARAM=NULL) {
    collected <- blockApply(DelayedArray(exprMat), FUN=.blocked_AUC_internal, 
        geneSets=geneSets,
        keepZeroesAsNA=keepZeroesAsNA,
        normAUC=normAUC,
        aucMaxRank=aucMaxRank,
        BPPARAM=BPPARAM,
        grid=colAutoGrid(exprMat)
    )
    do.call(cbind, collected)
}

.blocked_AUC_internal <- function(block, geneSets, keepZeroesAsNA, normAUC, aucMaxRank) {
    ranked <- .AUCell_buildRankings(block, plotStats=FALSE, keepZeroesAsNA=keepZeroesAsNA, verbose=FALSE)
    .AUCell_calcAUC(geneSets, ranked, normAUC=normAUC, aucMaxRank=aucMaxRank, verbose=FALSE)
}
