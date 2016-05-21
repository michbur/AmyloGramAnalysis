library(shiny)
library(seqinr)
library(biogram)
library(randomForest)
library(dplyr)
library(signalHsmm)
library(shinyAce)
library(markdown)
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
      predict_classifier(res)
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
      
      tabPanel("Paste sequences here:", aceEditor("text_area", value="", height = 150),
               actionButton("use_area", "Submit data from field above"),
               p(""),
               fileInput('seq_file', 'Submit .fasta or .txt file:'))
      
      
    } else {
      tabsetPanel(
        tabPanel("Summary", tableOutput("pred_table"))
      )
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