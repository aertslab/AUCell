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
#' @param exprMat Expression matrix (optional)
#' @param cellInfo Phenodata (optional)
#' @param colVars Color for the phenodata variables (as list, optional)
#' @param includeCellSelectionTab Whether to include the cell selection tab
#' @note 
#' With lasso: "To make a multiple selection, press the SHIFT key. To clear the selection, press the ESC key."
#' @return Thresholds and cells selected within the app (as list).
#' @example inst/examples/example_AUCell_createViewerApp.R
#' @export
AUCell_createViewerApp <- function(auc, thresholds=NULL, tSNE=NULL, 
                                   exprMat=NULL, cellInfo=NULL, colVars=NULL,
                                   includeCellSelectionTab=TRUE) # includeThresholdSelectionTab=TRUE, 
{
  if(!methods::is(auc,"aucellResults")) stop("Please provide an aucellResults object.")
  if(is.null(thresholds)) thresholds <- setNames(rep(0, nrow(auc)), rownames(auc))
  
  commonCells <- as.character(intersect(colnames(auc), rownames(tSNE)))
  tSNE.df <- data.frame(tSNE[commonCells,,drop=FALSE], cell=commonCells, t(getAUC(auc)[,commonCells, drop=FALSE]), stringsAsFactors=FALSE)
  colnames(tSNE.df)[1:2] <- c("tsne1","tsne2")
  #colnames(tSNE.df)[which(!colnames(tSNE.df) %in% c("tsne1", "tsne2", "cell", rownames(auc)))] # to add other props?
  
  app <- list()
  app$thresholds <- getThresholdSelected(thresholds)
  app$cells <- list()
  
  ################################################
  # UI
  
  # Choose according to whether the t-SNE is provided
  if(!is.null(tSNE))
  {
    if(!all(c("rbokeh") %in% rownames(installed.packages()))) # "shiny"
    {
      if(includeCellSelectionTab) warning("Package rbokeh is not available.")
    }else{
      #requireNamespace(rbokeh)
      #requireNamespace(shiny)
      
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
                     br(),
                     plotOutput(outputId = "histPlot")),
                   
                   mainPanel(plotOutput(outputId = "tsnePlot"),
                             # Extra properties (e.g. expression or cell info)
                             conditionalPanel(c("false","true")[as.numeric(!is.null(cellInfo) | !is.null(exprMat))+1],
                                              fluidRow(
                                                conditionalPanel(c("false","true")[as.numeric(!is.null(exprMat))+1],
                                                                 column(6,
                                                                        uiOutput("gene_selection")
                                                                 )),
                                                conditionalPanel(c("false","true")[as.numeric(!is.null(cellInfo))+1],
                                                                 column(6,
                                                                        selectInput(inputId = "phenodata_selection",
                                                                                    label = "Cell info:",
                                                                                    choices=colnames(cellInfo),
                                                                                    selected=colnames(cellInfo)[1])
                                                                 )),
                                                plotOutput(outputId = "tsnePlot_expression_cellInfo")
                                              )),
                             
                             
                             
                             sliderInput(inputId = "size",
                                         label = "Point size:",
                                         min = 0.01,
                                         max = 3,
                                         value = 0.5)
                   ) 
                   
          ),
          if(includeCellSelectionTab) 
            tabPanel("Cell selection", 
               column(6,
                      selectInput(inputId = "geneSetBokeh",
                                  label = "Gene set:",
                                  choices=rownames(auc)),
                      rbokeh::rbokehOutput("tsne_rbokeh"),
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
    }
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
    
    output$gene_selection <- renderUI ({
      possibleGene <- gsub( "\\s\\(\\d+g)", "",input$geneSet)
      possibleGene <- gsub( "_extended", "", possibleGene)
      
      gene <- ""
      if(possibleGene %in% rownames(exprMat)) gene <- possibleGene
      
      # selectInput(inputId = "geneExpression",
      #             label = "Gene expression:",
      #             choices=rownames(exprMat),
      #             selected=possibleGene)
      textInput(inputId = "geneExpression",
                label = "Gene expression:",
                value = gene)
    })
    
    # Reactive plots:
    output$histPlot <- renderPlot({
      AUCell_plotHist(auc[input$geneSet,], aucThr=input$threshold)
      abline(v=input$threshold)
    })
    
    output$tsnePlot <- renderPlot({
      # nRow <- 1
      # if(TRUE) nRow <- 2
      par(mfrow=c(1,2))
      # Binary
      passThreshold <- which(getAUC(auc)[input$geneSet,rownames(tSNE)]>input$threshold)
      .auc_plotBinaryTsne(tSNE, cex=input$size,
                          selectedCells=passThreshold,  
                          title=paste(input$geneSet), 
                          txt="Blue cells pass the threshold")  
      
      # Continuous
      .auc_plotGradientTsne(tSNE, cellProp=getAUC(auc)[input$geneSet,],
                            title=input$geneSet, txt="AUC value",
                            cex=input$size)
    })
    
    
    output$tsnePlot_expression_cellInfo <- renderPlot({
      if(is.null(exprMat) & is.null(cellInfo)){
        return(NULL)
      }else{
        
        par(mfrow=c(1,2))
        
        if(is.null(input$geneExpression) || input$geneExpression=="")
        {
          plot.new()
        }else{
          if(input$geneExpression %in% rownames(exprMat))
          {
            .auc_plotGradientTsne(tSNE, 
                                  cellProp=setNames(exprMat[input$geneExpression, rownames(tSNE)], rownames(tSNE)),
                                  title=paste0(input$geneExpression, " expression"), txt="",
                                  colorsForPal = c("goldenrod1", "darkorange", "brown"),
                                  cex=input$size)
            
            # Add legend (not included)
            legend(min(tSNE[,1]), max(tSNE[,2]), 
                   c("0", "", "", signif(max(exprMat[input$geneExpression, rownames(tSNE)]),2)), 
                   border="lightgrey",
                   fill=c("white", "goldenrod1", "darkorange", "brown"), 
                   box.lwd="none", bty = "n", cex=input$size*.8)
          }else{
            plot.new()
          }
        } 
        
        if(is.null(cellInfo))
        {
          plot.new()
        }else{
          .cellProps_plotTsne(tSNE, cellInfo, 
                              varName=input$phenodata_selection, cex=input$size,
                              colVars=colVars, sub="")
        }
      }
    })
    
    if("rbokeh" %in% rownames(installed.packages()) & includeCellSelectionTab)
    {
      output$tsne_rbokeh <- rbokeh::renderRbokeh({
        rbokeh::figure(logo=NULL) %>%
          rbokeh::ly_points(tsne1, tsne2, data=tSNE.df, hover=cell, size=input$size_bokeh,
                            color = getAUC(auc)[input$geneSetBokeh,rownames(tSNE.df)], legend=FALSE, lname = "cells") %>%
          rbokeh::set_palette(continuous_color = rbokeh::pal_gradient(c("lightgrey", "pink", "red"))) %>%
          rbokeh::tool_lasso_select(callback = rbokeh::shiny_callback(id="cellsSelected"), "cells")
      })
    }else{
      output$tsne_rbokeh <- NULL
    }
    
    # output$cellSelectedTable <- DT::renderDataTable({
    #   data.frame("Cells selected"= rownames(tSNE)[input$cellsSelected])
    # })
    
    output$cellSelectedText <- renderText({
      paste("Cells selected:\n", paste(rownames(tSNE)[input$cellsSelected+1], collapse="\n"), sep="")
    })
    
    # Save thresholds
    observeEvent(input$saveThr, {
      app$thresholds[input$geneSet] <<- input$threshold
      message(input$geneSet, " threshold replaced by ", app$thresholds[input$geneSet] )
    })
    observeEvent(input$saveCells, {
      app$cells[[input$cellGroupName]] <<- list(rownames(tSNE)[input$cellsSelected+1])
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
