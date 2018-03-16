\dontrun{
library(shiny); library(rbokeh)

# Create app 
aucellApp <- AUCell_createViewerApp(auc=cells_AUC, thresholds=selectedThresholds, 
                                    tSNE=cellsTsne, exprMat=exprMat)

# Run (the exact commands depend on the R settings, see Shiny's doc for help)
options(shiny.host="0.0.0.0") 
savedSelections <- runApp(aucellApp) 
}
