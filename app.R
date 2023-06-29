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
library(dplyr)
library(data.table)

# Getting the file names
rdsfiles <- list.files(pattern = "\\.Rds")
csvfiles <- list.files(pattern = "\\.csv")
types <- c("Astrocyte","Endothelia","Ependymal","GABAergic",
           "GLUTamatergic","Macrophage","Microglia","OPC","Oligodendrocyte")
bulk_types <- c("Astro","Endo","Epen","GABA",
           "GLUT","Macro","Micro","OPC","Oligo")
regions <- c("DLPFC","WM","Pulvinar")
# --- Front End ---
ui <- shinyUI(fluidPage(theme = shinytheme("spacelab"), pageWithSidebar(
  
  # Title
  headerPanel("Data Viewer"),
  
  # Sidebar to select a dataset
  sidebarPanel(
    selectInput("obj", "Choose a dataset:", #data_set is a seurat object
                choices = rdsfiles),
    selectInput("tissue_csv", 
              "Load tissue positions .csv file:",
              choices = csvfiles),
    selectInput("celltype", "Cell Type for Cell2Location", choices = types),
    textInput("feature", label = "sn Gene:")
    ),
  
  # Different analyses available
  mainPanel(
    tabsetPanel(
      tabPanel('UMAP', plotOutput("umap")),
      tabPanel('Tissue', plotOutput("tissue")),
      tabPanel('sn Gene Expression', plotOutput("genex")),
      tabPanel('Violin Plots', plotOutput("vln")),
      tabPanel('Cell2Location Prediction', plotOutput("c2l"))
      
      
    ),
    tabsetPanel(
      tabPanel("Pseudobulk RNA-seq", plotOutput("bulk"))
    ),
    tabsetPanel(
      selectInput("bulk_cell","Cell Type", choices = bulk_types),
      selectInput("bulk_region","Brain Region", choices = regions),
      textInput("genes_list", "Gene(s)")
    )
    )
  
)))

# --- Back end ---
server <- shinyServer(function(input, output, session) {
  
  # Return the requested datasets
  datasetInput <- reactive({
    df <- readRDS(input$obj, input$obj)
    return(df)
  })
  
  # Return the tissue coordinate
  tissueInput <- reactive({
    inFile <- req(input$tissue_csv)
    return(inFile)
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
  
  bulkcellInput <- reactive({
    ct <- req(input$bulk_cell)
    return(ct)
  })
  
  bulkRegionInput <- reactive({
    br <- (input$bulk_region)
    return(br)
  })
  
  bulkGenes <- reactive({
    bg <- req(input$genes_list)
  })

  #Generate the tissue plot
  output$tissue <- renderPlot({
    obj <- datasetInput()
    tiss <- tissueInput()
    transfer_clusters(obj, tiss)
    
 })
  # Retrieve the UMAP projection
  output$umap <- renderPlot({
    obj <- datasetInput()
    DimPlot(obj, reduction = "umap")
    
    
  })
  # Plot the expression of a gene
  output$genex <- renderPlot({
    obj <- datasetInput()
    tiss <- tissueInput()
    feat <- featInput()
    get_expression(obj, feat, tiss)
    
  })
  #Generate violin plots
  output$vln <- renderPlot({
    obj <- datasetInput()
    feat <- featInput()
    VlnPlot(obj, features = feat, group.by = "seurat_clusters")
  })
  
  output$c2l <- renderPlot({
    sample <- sampleInput()
    csv <- tissueInput()
    pred <- "cell2loc_broad_preds_norm.csv"
    celltype <- cellInput()
    cell2loc(sample,pred,csv,celltype)
    
    
    
    
  })
  
  output$bulk <- renderPlot({
    cell <- bulkcellInput()
    region <- bulkRegionInput()
    gene_list <- bulkGenes()
    bulk_plot(cell, region, gene_list)
  })

  
})


shinyApp(ui, server)