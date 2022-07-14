library(shiny)
options(shiny.maxRequestSize=3000*1024^2) #3 GB

#Function to retrieve available data sets in the suerat object
list_datasets <- function(obj){
  sets <- unique(obj@meta.data$dataset)
}

# ---Frontend---
ui <- fluidPage(
  # App title
  titlePanel("SeuratPlusVisium Data Viewer"),
  
  # File Upload
  fileInput(inputId = "upload", 
            label = "Upload Seurat Object in .Rds format"),
  )

# ---Backend---
require(Seurat)
require(seuratplusvisium)
require(ggplot2)

server <- function(input, output, session) {
  
  observe({
    inFile <- input$upload
    if (is.null(inFile))
      return(NULL)
    seurat <- readRDS(inFile)
    sets <- list_datasets(inFile)
  })
}
shinyApp(ui, server)