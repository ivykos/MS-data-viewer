library(shiny)
library(shinythemes)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(egg)
library(RColorBrewer)
library(shinyWidgets)
library(tidyverse)
library(tools)
library(rlang)
library(dplyr)
#library(ggtips)


# Getting the file names
##### Okay, what we need to do is download all the data and
##### separate it into CTRL and MS groups
rdsfiles <- list.files(pattern = "\\.Rds$")
csvfiles <- list.files(pattern = "\\.csv$" )
types <- c("Astrocyte","Endothelia","Ependymal","GABAergic",
           "GLUTamatergic","Macrophage","Microglia","OPC","Oligodendrocyte")

# --- Front End ---
ui <- shinyUI(fluidPage(theme = shinytheme("spacelab"), pageWithSidebar(
  
  # Title
  headerPanel("Data Viewer"),
  
  # Sidebar to select a dataset
  sidebarPanel(
    selectInput("MS_obj", "MS dataset:", #data_set is a seurat object
                choices = rdsfiles),
    selectInput("CTRL_obj", "Control dataset",
                choices = rdsfiles),
    selectInput("MS_tissue_csv", 
              "MS Tissue Coordinates File:",
              choices = csvfiles),
    selectInput("CTRL_tissue_csv", "Control Tissue Coordinates File",
                choices = csvfiles),
    selectInput("celltype", "Cell Type for Cell2Location", choices = types),
    textInput("feature", label = "Gene:"),
    
    
  ),
  
  # Different analyses available
  mainPanel(
    tabsetPanel(
      tabPanel('UMAP', plotOutput("umap")),
      tabPanel('Tissue', plotOutput("tissue")),
      tabPanel('Gene Expression', plotOutput("genex")),
      tabPanel('Violin Plots', plotOutput("vln")),
      tabPanel('Cell2Location Prediction', plotOutput("c2l"))
      
      
    ))
  
)))

# --- Back end ---
server <- shinyServer(function(input, output, session) {
  
  # Return the requested datasets
  MSdatasetInput <- reactive({
    MSdf <- readRDS(input$MS_obj)
    return(MSdf)
  })
  
  CTRLdatasetInput <- reactive({
    CTRLdf <- readRDS(input$CTRL_obj)
    return(CTRLdf)
  })
 
  # Return the tissue coordinates
  MStissueInput <- reactive({
    MSinFile <- req(input$MS_tissue_csv)
    return(MSinFile)
  })
  
  CTRLtissueInput <- reactive({
    CTRLinFile <- req(input$CTRL_tissue_csv)
    return(CTRLinFile)
  })
  # Return requested gene or feature
  featInput <- reactive({
    text <- req(input$feature)
    return(text)
    
  })
  
  # Return the cell type
  cellInput <- reactive({
    cell <- req(input$celltype)
    return(cell)
  })
  
  #Return sample name
  sampleInput <- reactive({
    samp <- tools::file_path_sans_ext(as.character(tissueInput()))
    return(samp)
  })
  
  MSsampleInput <- reactive({
    samp <- tools::file_path_sans_ext(as.character(MStissueInput()))
    return(samp)
  })
  
  CTRLsampleInput <- reactive({
    samp <- tools::file_path_sans_ext(as.character(CTRLtissueInput()))
    return(samp)
  })

  #Generate the tissue plot
  output$tissue <- renderPlot({
    obj <- MSdatasetInput()
    tiss <- MStissueInput()
    obj2 <- CTRLdatasetInput()
    tiss2 <- CTRLtissueInput()
    tissue_patch(obj, tiss, obj2, tiss2)
    
 })
  # Retrieve the UMAP projection
  output$umap <- renderPlot({
    obj <- MSdatasetInput()
    obj2 <- CTRLdatasetInput()
    umap_patch(obj, obj2)
    
    
  })
  # Plot the expression of a gene
  output$genex <- renderPlot({
    obj <- MSdatasetInput()
    tiss <- MStissueInput()
    feat <- featInput()
    obj2 <- CTRLdatasetInput()
    tiss2 <- CTRLtissueInput()
    gene_patch(obj, obj2, tiss, tiss2, feat)
  })
  
  #Generate violin plots
  output$vln <- renderPlot({
    obj <- MSdatasetInput()
    obj2 <- CTRLdatasetInput()
    feat <- featInput()
    vln_patch(obj, obj2, feat)
  })
  
  output$c2l <- renderPlot({
    sample <- MSsampleInput()
    sample2 <- CTRLsampleInput()
    csv <- MStissueInput()
    csv2 <- CTRLtissueInput()
    pred <- "cell2loc_broad_preds_norm.csv"
    celltype <- cellInput()
    #cell2loc(sample,pred,csv,celltype)
    c2l_patch(sample, sample2, csv, csv2, celltype, pred)
  })
})


shinyApp(ui, server)