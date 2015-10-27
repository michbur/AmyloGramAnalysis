library(shiny)
library(ggplot2)
library(DT)
load("properties.RData")

props_order <- c(1, 6, 10, 13, 16, 17, 18, 19, 20, #hydroph
                 3, 9, #polarity
                 2, 5, 8, 14, 21:26, #contact
                 4, 7, 11, #beta prob 
                 12, #stability
                 15) #size


my_theme <- theme(plot.background=element_rect(fill = "transparent",
                                               colour = "transparent"),
                  panel.grid.major = element_line(colour="grey", linetype = "dashed", size = 0.5),
                  panel.grid.major = element_line(colour="lightgrey", linetype = "dashed", size = 0.5),
                  panel.background = element_rect(fill = "transparent",colour = "black"),
                  legend.background = element_rect(fill = "NA"),
                  legend.position = "bottom",
                  strip.background = element_rect(fill = "NA", colour = "NA"))

shinyServer(function(input, output, session) {
  
  observe({
    if(input[["selectall"]] != 0) {
      if (input[["selectall"]]%%2 == 0) {
        updateCheckboxGroupInput(session, "checkProp", label = h3("Choose 1 or more properties"), 
                                 choices = as.list(plot_id)[props_order],
                                 selected = as.character(1))
      } else {
        updateCheckboxGroupInput(session, "checkProp", label = h3("Choose 1 or more properties"), 
                                 choices = as.list(plot_id)[props_order],
                                 selected = as.character(1L:26))
      }
    }
  })
  
  output$value <- DT::renderDataTable({ 
    tab <- t(plot_values[, as.numeric(input$checkProp)])
    rownames(tab) <- names(plot_id)[as.numeric(input$checkProp)]
    datatable(tab, escape = FALSE, extensions = 'TableTools', 
              filter = "top", options = list(
                dom = 'T<"clear">lfrtip',
                tableTools = list(sSwfPath = copySWF('www'))))
  })
  
  
  output$plot <- renderPlot({
    ggplot(mplot_values[mplot_values[["id"]] %in% input$checkProp, ], aes(x = id, fill = id, y = value)) +
      geom_bar(stat="identity", position = "dodge") + 
      facet_wrap(~ aa, ncol = 5) +
      scale_x_discrete("Amino acid") + 
      scale_y_continuous("Normalized value") +
      scale_fill_discrete(name="Property name", 
                          labels = names(plot_id)[as.numeric(input$checkProp)]) +
      guides(fill=guide_legend(ncol=2)) + my_theme
    
  })
  
  corr_df <- reactive({
    melt(cor(plot_values[, as.numeric(input$checkProp)] - 0.5))
  })
  
  
  output$corplot <- renderPlot({
    nms <- names(plot_id)[as.numeric(input$checkProp)]
    corm <- corr_df()
    
    br_nms <- sapply(nms, function(single_nm) {
      len <- nchar(single_nm)
      space_pos <-  which(strsplit(single_nm, "")[[1]] == " ")
      break_pos <- space_pos[which.min(abs(len/2 - space_pos))]
      paste0(substr(single_nm, 1, break_pos - 1), "\n", substr(single_nm, break_pos + 1, len))
    })
    
    corm[["Var1"]] <- factor(corm[["Var1"]], labels = br_nms)
    corm[["Var2"]] <- factor(corm[["Var2"]], labels = br_nms)
    corm[["value_fac"]] <- round(corm[["value"]], 4)
    
    ggplot(data = corm, aes(x=Var1, y=Var2, fill=value, label = value_fac)) + 
      geom_tile() +
      geom_text() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(-1,1), name="Correlation\ncoefficient") +
      my_theme
  })
  
  thr_list <- reactive({
    corm <- corr_df()
    nms <- names(plot_id)[as.numeric(input$checkProp)]
    
    corm[["Var1"]] <- factor(corm[["Var1"]], labels = nms)
    corm[["Var2"]] <- factor(corm[["Var2"]], labels = nms)

    lapply(levels(corm[["Var1"]]), function(i)
      corm[corm[["Var1"]] == i & corm[["Var2"]] != i & abs(corm[["value"]]) > input[["thresh"]], ])
  })
  
  
  output$cortab <- DT::renderDataTable({ 
    nms <- names(plot_id)[as.numeric(input$checkProp)]
    
    res <- data.frame(prop_name = nms, n_similar = sapply(thr_list(), nrow, USE.NAMES = FALSE))
    datatable(res, escape = FALSE, extensions = 'TableTools', 
              filter = "top", options = list(
                dom = 'T<"clear">lfrtip',
                tableTools = list(sSwfPath = copySWF('www')), pageLength = 30, server = TRUE))
  })
  
  output$cortext <- renderPrint({ 
    nms <- names(plot_id)[as.numeric(input$checkProp)]
    res <- sapply(thr_list(), function(i) paste0(as.character(i[["Var2"]]), collapse = ",\n"))
    names(res) <- nms
    chosen_props <- as.numeric(input[["cortab_rows_selected"]])
    for(i in chosen_props) {
      cat("Property:", names(res)[i], "\n\n")
      cat(res[i], "\n\n")
    }
  })
  
  output$pcaplot <- renderPlot({
    pca_res <- princomp(plot_values[, as.numeric(input$checkProp)])
    ggdat <- cbind(data.frame(pca_res[["scores"]]), aa = rownames(pca_res[["scores"]]))
    ggplot(ggdat, aes(x = Comp.1, y = Comp.2, colour = Comp.3, label = aa)) +
      geom_point(size = 5) +
      scale_colour_gradient2(low = "red", high = "yellow", mid = "orange") +
      my_theme +
      geom_text(hjust=1.8, vjust=1.8)
  })
  
})
