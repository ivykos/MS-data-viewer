library(shiny)
library(shinythemes)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(egg)
library(RColorBrewer)
library(shinyWidgets)
library(tidyverse)


# Getting the file names
rdsfiles <- list.files(pattern = "\\.Rds$")
csvfiles <- list.files(pattern = "\\.csv$" )
# --- Front End ---
ui <- shinyUI(fluidPage(theme = shinytheme("cerulean"), pageWithSidebar(
  
  # Title
  headerPanel("Seurat + Visium Data Viewer"),
  
  # Sidebar to select a dataset
  sidebarPanel(
    selectInput("obj", "Choose a dataset:", #data_set is a seurat object
                choices = rdsfiles),
    selectInput("tissue_csv", 
              "Load tissue positions .csv file",
              choices = csvfiles),
    textInput("feature", label = "Gene"),
    
    
  ),
  
  # Different analyses available
  mainPanel(
    tabsetPanel(
      tabPanel('UMAP', plotOutput("umap")),
      tabPanel('Tissue', plotOutput("tissue")),
      tabPanel('Gene Expression', plotOutput("genex")),
      #tabPanel("Test", tableOutput("contents"))
    ))
  
)))

# --- Back end ---
server <- shinyServer(function(input, output) {
  
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
  
})


shinyApp(ui, server)