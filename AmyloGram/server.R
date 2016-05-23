library(shiny)
library(markdown)
library(DT)
source("functions.R")
load("AmyloGram.RData")

options(shiny.maxRequestSize=10*1024^2)

options(DT.options = list(dom = "Brtip",
                          buttons = c("copy", "csv", "excel", "print")
))

my_DT <- function(x)
  datatable(x, escape = FALSE, extensions = 'Buttons', 
            filter = "top", rownames = FALSE)


shinyServer(function(input, output) {
  
  prediction <- reactive({
    
    if (!is.null(input[["seq_file"]]))
      input_sequences <- read_txt(input[["seq_file"]][["datapath"]])
    input[["use_area"]]
    isolate({
      if (!is.null(input[["text_area"]]))
        if(input[["text_area"]] != "")
          input_sequences <- read_txt(textConnection(input[["text_area"]]))
    })
    
    if(exists("input_sequences")) {
      if(length(input_sequences) > 100) {
        #dummy error, just to stop further processing
        stop("Too many sequences.")
      } else {
        predict_AmyloGram(AmyloGram_model, input_sequences)
      }
    } else {
      NULL
    }
  })
  
  decision <- reactive({
    if(!is.null(prediction()))
      make_decision(prediction(), 0.5)
  })
  
  
  output$dynamic_ui <- renderUI({
    if(!is.null(prediction())) {
          tags$p("Refresh page (press F5) to start a new query with signalHsmm.")
    }
  })
  
  output$pred_table <- DT::renderDataTable({
    formatRound(my_DT(decision()), 2, 4)
  })
  
  output$dynamic_tabset <- renderUI({
    if(is.null(prediction())) {
      
      tabPanel(title = "Sequence input",
               tags$textarea(id = "text_area", style = "width:100%",
                             placeholder="Paste sequences (FASTA format required) here...", rows = 22, cols = 60, ""),
               p(""),
               actionButton("use_area", "Submit data from field above"),
               p(""),
               fileInput('seq_file', 'Submit .fasta or .txt file:'))
      
      
    } else {
      tabPanel("Short output", DT::dataTableOutput("pred_table"))
    }
  })
  

})