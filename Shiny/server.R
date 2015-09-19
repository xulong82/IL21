library(shiny)
library(igraph)
library(ggplot2)

load("gene_expression_and_summary_10922.rdt")
load("igraph.rdt")

shinyServer(function(input, output) {
  
  output$igraph <- renderPlot({
    if (input$cell == 1) igraph.dt <- igraphList$NN
    if (input$cell == 2) igraph.dt <- igraphList$NP
    if (input$cell == 3) igraph.dt <- igraphList$PP
    
    if (input$layout == 1) igraph.dt$layout <- layout.sphere
    if (input$layout == 2) igraph.dt$layout <- layout.circle
    if (input$layout == 3) igraph.dt$layout <- layout.fruchterman.reingold 
    
    plot.igraph(igraph.dt, vertex.frame.color = "white", edge.arrow.size = 0.3)
  }, height = 800)
  
  output$plot <- renderPlot({
    geneId <- input$gene
    
    if (geneId %in% rownames(myTpm))
      expr1 <- t(myTpm[geneId, ])
    else
      expr1 <- rep(0, 6)
    
    cell <- c("Naive", "Naive", "Act", "Act", "Act-Il21", "Act-Il21")
    graph.dt <- data.frame(value = expr1, 
      cell = factor(cell, levels = c("Naive", "Act", "Act-Il21")))
    rownames(graph.dt) <- NULL
    colnames(graph.dt) <- c("value", "cell")
    
    ggplot(graph.dt, aes(x = cell, y = value)) + 
      geom_boxplot(aes(color = cell)) +
      geom_point(colour = "grey50", size = 5) +
      geom_point(aes(color = cell), size = 3) +
      xlab("") + ylab("") + ggtitle(input$gene) +
      scale_color_manual(values = c("grey30", "dodgerblue3", "firebrick1")) +
      theme(axis.title.y = element_text(vjust = 1, color = "grey30"),
            plot.title = element_text(size=20, face="bold", vjust=2))
  })

  output$summary <- renderTable({
    geneId <- input$gene
    table[table$query == geneId, ]
  }, include.rownames = FALSE)

  output$downloadData <- downloadHandler(
    filename = "TPM.xlsx", content = function (file) file.copy("./TPM.xlsx", file)
  )

  output$downloadDataM1 <- downloadHandler(
    filename = "Supp_Table1.xlsx", content = function (file) file.copy("./Supp_Table1.xlsx", file)
  )

  output$downloadDataM2 <- downloadHandler(
    filename = "Supp_Table2.xlsx", content = function (file) file.copy("./Supp_Table2.xlsx", file)
  )

})
