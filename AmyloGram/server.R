library(shiny)
library(markdown)
library(DT)
source("functions.R")
load("AmyloGram.RData")

options(shiny.maxRequestSize=10*1024^2)

shinyServer(function(input, output) {
  
  prediction <- reactive({
    
    if (!is.null(input[["seq_file"]]))
      res <- read_txt(input[["seq_file"]][["datapath"]])
    input[["use_area"]]
    isolate({
      if (!is.null(input[["text_area"]]))
        if(input[["text_area"]] != "")
          res <- read_txt(textConnection(input[["text_area"]]))
    })
    
    if(exists("res")) {
      if(length(res) > 300) {
        #dummy error, just to stop further processing
        stop("Too many sequences.")
      } else {
        predict_AmyloGram(AmyloGram_model, res)
      }
    } else {
      NULL
    }
  })
  
  output$dynamic_ui <- renderUI({
    if(!is.null(prediction())) {
      div(tags$h3("Download results"),
          tags$p(""),
          downloadButton("download_short", "Download output"),
          tags$p("Refresh page (press F5) to start a new query with signalHsmm."))
    }
  })
  
  
  output$pred_table <- renderTable({
    prediction()
  })
  
  output$dynamic_tabset <- renderUI({
    if(is.null(prediction())) {
      
      tabPanel(title = "Sequence input",
               h3("Paste sequences (FASTA format required) into the field below:"), 
               tags$style(type="text/css", "textarea {width:100%}"),
               tags$textarea(id = "text_area", rows = 22, cols = 60, ""),
               p(""),
               actionButton("use_area", "Submit data from field above"),
               p(""),
               fileInput('seq_file', 'Submit .fasta or .txt file:'))
      
      
    } else {
      tabPanel("Short output", tableOutput("pred_table"))
    }
  })
  
  #name for downloads
  file_name <- reactive({
    if(is.null(input[["seq_file"]][["name"]])) {
      part_name <- "AmyloGram_results"
    } else {
      part_name <- strsplit(input[["seq_file"]][["name"]], ".", fixed = TRUE)[[1]][1]
    }
    part_name
  })
  
  
  output$download_short <- downloadHandler(
    filename  = function() { 
      paste0(file_name(), "_pred.csv") 
    },
    content <- function(file) {
      write.csv(prediction(), file)
    }
  )
  
  
})