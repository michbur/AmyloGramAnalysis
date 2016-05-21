library(shiny)

shinyUI(fluidPage(
  
  headerPanel("AmyloGram"),
  
  sidebarLayout(
    sidebarPanel(
      includeMarkdown("readme.md"),
      pre(includeText("prots.txt")),
      uiOutput("dynamic_ui")
    ),
    
    mainPanel(
      uiOutput("dynamic_tabset")    
    )
  )))
