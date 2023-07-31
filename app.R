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
      tabPanel('Cell2Location Spatially Predicted Proportions', plotOutput("c2l"))
      
      
    ),
    tabsetPanel(
      tabPanel("Pseudobulk RNA-seq", plotOutput("bulk"))
    ),
    tabsetPanel(
      selectInput("bulk_cell","Cell Type", choices = bulk_types),
      selectInput("bulk_region","Brain Region", choices = regions),
      textInput("genes_list", "Gene(s)")
    ),
    tabsetPanel(
      tabPanel("Cell2Location Predicted Proportions", imageOutput("image"))
    ),
    tabsetPanel(
      selectInput("celltypes", "Cell Type", choices = types),
      selectInput("region", "Region", choices = regions)
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
  
  #Return cell type for bulk RNA-seq visualization
  bulkcellInput <- reactive({
    ct <- req(input$bulk_cell)
    return(ct)
  })
  
  #Return specified region for bulk RNA-seq visualization
  bulkRegionInput <- reactive({
    br <- (input$bulk_region)
    return(br)
  })
  
  #Return desired gene for bulk RNA-seq visualization
  bulkGenes <- reactive({
    bg <- req(input$genes_list)
  })
  
  #Return specified celltype for cell2loc region predictions
  image_celltype <- reactive({
    ct2 <- input$celltypes
    return(ct2)
  })
  
  #Return specified region for cell2loc region predictions
  image_region <- reactive({
    reg <- input$region
    return(reg)
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
  #output$c2lp <- renderPlot({
  #  obj <- datasetInput()
    #feat <- featInput()
    #VlnPlot(obj, features = feat, group.by = "seurat_clusters")
  #})
  
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
  
  output$image <- renderImage({
    cell <- image_celltype()
    region <- image_region()
    path_to_image <- as.character(paste(cell, region, ".png", sep = "_"))
    list(
      src = file.path(path_to_image)
    )
  }, deleteFile = FALSE)
  
  

  
})


shinyApp(ui, server)