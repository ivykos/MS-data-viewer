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

# --- Front End ---
ui <- shinyUI(fluidPage(theme = shinytheme("cerulean"), pageWithSidebar(
  
  # Title
  headerPanel("Seuratplusvisium Data Viewer"),
  
  # Sidebar to select a dataset
  sidebarPanel(
    selectInput("dataset", "Choose a dataset:", 
                choices = rdsfiles),
    textInput("feature", label = "Gene")
  ),
  
  # 
  mainPanel(
    tabsetPanel(
      tabPanel('Tissue', plotOutput("tissue")),
      tabPanel('UMAP', plotOutput("umap")),
      tabPanel('Genes', plotOutput("genex"))
      
    ))
  
)))

# --- Back end ---
server <- shinyServer(function(input, output) {
  
  # Return the requested dataset
  datasetInput <- reactive({
    df <- readRDS(input$dataset, input$dataset)
    return(df)
  })
  
  # Retrieve the UMAP projection
  output$umap <- renderPlot({
    dataset <- datasetInput()
    DimPlot(dataset, reduction = "umap")
    
  })
  
  # 
  output$gene <- renderPlot({
    dataset <- datasetInput()
    FeaturePlot(dataset, reduction = "umap")
    
  })
  
})


shinyApp(ui, server)