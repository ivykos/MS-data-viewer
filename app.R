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
  headerPanel("Seuraplusvisium Data Viewer"),
  
  # Sidebar to select a dataset
  sidebarPanel(
    selectInput("dataset", "Choose a dataset:", 
                choices = rdsfiles),
  ),
  
  # 
  mainPanel(
    tabsetPanel(
      tabPanel('UMAP', plotOutput("umap")),
      tabPanel('Gene', plotOutput("gene"))
      
    ))
  
)))

# --- Back end ---
server <- shinyServer(function(input, output) {
  
  # Return the requested dataset
  datasetInput <- reactive({
    df <- readRDS(paste0("/temp2/data/", input$dataset))
    return(df)
  })
  
  # Generate a UMAP of the dataset
  output$umap <- renderPlot({
    dataset <- datasetInput()
    plot(input$dataset, reduction = "umap")
    
  })
  
  # Generate a Feature of the dataset
  output$gene <- renderPlot({
    dataset <- datasetInput()
    FeaturePlot(dataset, reduction = "umap")
    
  })
  
})


shinyApp(ui, server)