# Help files will be automatically generated from the coments starting with #'
# (https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html)

#' @import shiny
#'
#' @title Create AUCell viewer app
#' @description Creates a Shiny app to explore AUCell results
#' @param auc AUC object returned by \code{\link{AUCell_calcAUC}}
#' @param thresholds Thresholds corresponding to each gene set (optional)
#' @param tSNE t-SNE coordinates for the cells (optional). 
#' The row names should correspond to the cell ID. 
#' The column names should be "tsne1" and "tsne2".
#' @return Thresholds and cells selected within the app (as list).
#' @examples
#' \dontrun{
#' # Create the Shiny app: 
#' aucellApp <- AUCell_createViewerApp(auc=cells_AUC, 
#'                thresholds=selectedThresholds,
#'                tSNE=cellsTsne)
#'                
#' # This object contains the $ui and $server required to lauch the app, e.g.:
#' savedSelections <- shinyApp(ui=aucellApp$ui, server=aucellApp$server)
#' 
#' # (How to launch the app might depend on the local settings)
#' options(shiny.host="0.0.0.0") 
#' savedSelections <- runApp(aucellApp) 
#' }
#' @export
AUCell_createViewerApp <- function(auc, thresholds=NULL, tSNE=NULL)
{
  library(shiny)
  if(class(auc)!="aucellResults") stop("Please provide an aucellResults object.")
  if(is.null(thresholds)) thresholds <- setNames(rep(0, nrow(auc)), rownames(auc))
  
  tSNE.df <- data.frame(tSNE, cell=rownames(tSNE), t(getAUC(auc)[,rownames(tSNE)]))
  #colnames(tSNE.df)[which(!colnames(tSNE.df) %in% c("tsne1", "tsne2", "cell", rownames(auc)))] # to add other props?

  app <- list()
  app$thresholds <- thresholds
  app$cells <- list()
  
  ################################################
  # UI

  # Choose according to whether the t-SNE is provided
  if(!is.null(tSNE))
  {
    library(rbokeh)
    app$ui <- fluidPage(
        titlePanel("AUCell"),
        tabsetPanel(
          tabPanel("Threshold selection", 
                   sidebarPanel(
                     selectInput(inputId = "geneSet",
                                 label = "Gene set:",
                                 choices=rownames(auc)
                     ),
                     uiOutput("threshold_slider"),
                     actionButton("saveThr", "Save threshold"),
                     plotOutput(outputId = "histPlot")),
                   
                   mainPanel(plotOutput(outputId = "tsnePlot2"),
                             sliderInput(inputId = "size",
                                         label = "Point size:",
                                         min = 0.01,
                                         max = 3,
                                         value = 0.5)
                   ) 
                   
          ),
          tabPanel("Cell selection", 
                   column(6,
                          selectInput(inputId = "geneSetBokeh",
                                      label = "Gene set:",
                                      choices=rownames(auc)),
                          rbokehOutput("tsne_rbokeh"),
                          sliderInput(inputId = "size_bokeh",
                                       label = "Point size:",
                                       min = 0.01,
                                       max = 10,
                                       value = 1)
                          ),
                   column(6,
                          wellPanel(
                          fixedRow(textInput(inputId = "cellGroupName", 
                                             label="Group name", value = "group1"),
                          actionButton("saveCells", "Save cells")),
                          textOutput("cellSelectedText")))
                   #column(6, DT::dataTableOutput("cellSelectedTable"))
          )
        ),
      title="AUCell"
    )
  }else{
    # If no t-SNE: Histogram on the main panel
    app$ui <- fluidPage(
      titlePanel("AUCell"),
      sidebarPanel(
        selectInput(inputId = "geneSet",
                    label = "Gene set:",
                    choices=rownames(auc)
        ),
        # Input: Slider for the threshold ----
        uiOutput("threshold_slider"),
        actionButton("saveThr", "Save threshold")
        ),
      
      mainPanel(
        plotOutput(outputId = "histPlot")
      ),
      title="AUCell"
    )
  }
  
  ################################################
  # Server
  app$server <- function(input, output, session) {
    
    # Reactive inputs: 
    output$threshold_slider <- renderUI ({
      sliderInput(inputId = "threshold",
                  label = "Threshold:",
                  min = signif(max(min(getAUC(auc)[input$geneSet,])-0.01,0),2),
                  max = signif(max(getAUC(auc)[input$geneSet,])+0.01,2),
                  value = app$thresholds[input$geneSet])
    })
    
    # Reactive plots:
    output$histPlot <- renderPlot({
      AUCell_plot(auc[input$geneSet,], aucThr=max(getAUC(auc)[input$geneSet,])+0.01) 
      abline(v=input$threshold)
    })
    
    output$tsnePlot2 <- renderPlot({
      par(mfrow=c(1,2))
      
      # Binary
      plot(tSNE, main=paste(input$geneSet),
           sub="Blue cells pass the threshold",
           col=c("#e0e0e020","blue")[as.numeric(getAUC(auc)[input$geneSet,rownames(tSNE)]>input$threshold)+1],
           pch=16, cex=input$size)

      # Continuous
      nBreaks <- 5 # Number of levels in the color palettes
      colorPal <- grDevices::colorRampPalette(c("lightgray", "pink", "red"))(nBreaks)
      
      # Assign cell color
      cellColor <- setNames(rep("black", nrow(tSNE)), rownames(tSNE))
      cellColor[colnames(auc)] <- setNames(colorPal[cut(getAUC(auc)[input$geneSet,], breaks=nBreaks)], colnames(auc))
      
      # Plot
      plot(tSNE, main=input$geneSet,
           sub="AUC value",
           col=cellColor[rownames(tSNE)], pch=16, cex=input$size) 
      
    })
    
    output$tsne_rbokeh <- renderRbokeh({
      figure() %>%
        ly_points(tsne1, tsne2, data=tSNE.df, hover=cell, size=input$size_bokeh,
            color = getAUC(auc)[input$geneSetBokeh,rownames(tSNE.df)], legend=FALSE, lname = "cells") %>%
      set_palette(continuous_color = pal_gradient(c("lightgrey", "pink", "red"))) %>%
        tool_lasso_select(callback = shiny_callback(id="cellsSelected"), "cells")
    })
    
    output$cellSelectedTable <- DT::renderDataTable({
      data.frame("Cells selected"= rownames(tSNE)[input$cellsSelected])
    })
    
    output$cellSelectedText <- renderText({
      paste("Cells selected:\n", paste(rownames(tSNE)[input$cellsSelected], collapse="\n"), sep="")
    })
    
    # Save thresholds
    observeEvent(input$saveThr, {
      app$thresholds[input$geneSet] <<- input$threshold
      message(input$geneSet, " threshold replaced by ", app$thresholds[input$geneSet] )
    })
    observeEvent(input$saveCells, {
      app$cells[[input$cellGroupName]] <<- list(rownames(tSNE)[input$cellsSelected])
      message("Selected cells (", input$cellGroupName,"): ", app$cells[[input$cellGroupName]])
      if(grepl("group([[:digit:]]+)", input$cellGroupName)) {
        updateTextInput(session,
          inputId="cellGroupName", 
          value = paste0("group", as.numeric(gsub("group", "", input$cellGroupName))+1))
      }
    })
    
    # Return selections
    onStop(function() {
      message("App stopped. Returning the thresholds & cells selected.")
      stopApp(returnValue = app[c("thresholds","cells")]) 
    })
  } 
  
  return(app)
} 



