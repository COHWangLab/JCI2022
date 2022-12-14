library(shiny)
library(shinythemes)
source("helpers.R")
#Refdata <- readRDS("data/MCL_cluster_reference.rds")




ui <- fluidPage(
   theme = "cerulean", 
   titlePanel("MCL Cluster"),
  
  sidebarLayout(
    
    position = 'left',
    
    sidebarPanel(
      helpText("Project your data to clusters"),
      fileInput("file", label = h3("File input"),placeholder = "No file selected"),
      hr("Please prepare table as data/example_input.txt !", n = 80),
      fluidRow(column(6, verbatimTextOutput("value"))),
      sliderInput("PlotWidth", "PlotWidth:", min = 2, max = 20, value = 5,
                  step = 1),
      sliderInput("PlotHeight", "PlotHight:", min = 2, max = 20, value = 5,
                  step = 1),
      sliderInput("RefDotSize", "RefDotSize:", min = .5, max = 10, value = 2,
                  step = .5),
      sliderInput("InDotSize", "InDotSize:", min = .5, max = 10, value = 4,
                  step = .5),
      submitButton(text = "Submit"),
      
    ),
    
    mainPanel(
      
      tabsetPanel(
        
        tabPanel("Plot", plotOutput("plot"),
                 downloadButton('downloadPlot', 'Download Plot')
                 ),
        
        tabPanel("Table", tableOutput("table"),
                 downloadButton('downloadTable', 'Download Data')
                 )
        
      )
    
    )
  
))




server <- function(input, output) {
  
  dataInput <- reactive({
    
    inFile <- input$file
    
    if (is.null(inFile)){
      
      return(Refdata$ref)
      
    }else{
      
    dataInput <- read_input(inFile)
    return(dataInput)
  }
    })
  
  output$value <- renderPrint({
    n <- ncol( dataInput() ) - 150
    print( paste0(n," samples are loaded."))
    sample_name <- colnames( dataInput()[,-c(1:150)] )
    print( paste0("Samples are : ",paste(sample_name, collapse = " ")))
  })
  

  UmapNew <- reactive({ 
    UmapNew <- project_sample( dataInput() ) 
    return(UmapNew) 
    
    })
  
  plotInput <- reactive({
    p <- plot_umap( umap_new = UmapNew(), 
                    RefDot = input$RefDotSize,
                    InDot = input$InDotSize )
    return(p)
    })
  
  output$table <- renderTable({ UmapNew() })
  output$downloadTable <- downloadHandler(
    filename = function() { paste(input$PlotFileName, '.tsv', sep='') },
    content = function(file) { write.table( UmapNew() , file ,sep='\t',quote=F,row.names = F) }
  )
  
  output$plot <- renderPlot({  plotInput() })
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste('UMAPPlot_', Sys.Date(), '.pdf', sep='')
    },
    content = function(file) {
      ggsave(file, plotInput(), width = input$PlotWidth, height = input$PlotHeight, dpi = 600, units = "in")
    }
  ) 
  
  
}

shinyApp(ui = ui, server = server)
