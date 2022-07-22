library(shiny)
library(shinythemes)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(egg)
library(RColorBrewer)
library(shinyWidgets)


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
  tissueInput <- reactive({
    inFile <- req(input$tissue_csv)
    return(inFile)
  })
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
  
  
})


shinyApp(ui, server)