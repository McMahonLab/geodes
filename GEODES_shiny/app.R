#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(reshape2)
library(cowplot)

# Read and process Mendota data
sig.mendota.key <- read.csv("data/WGCNA_mendota_results.csv", header = T)
mendota.eigenvectors <- read.csv("data/WGCNA_mendota_eigenvectors.csv", header = T, row.names = 1)
clusters <- c(2, 3, 6, 7, 8, 13, 21, 25, 26, 30, 36, 37, 39, 47, 53, 55)
plot.sig.mendota.key <- sig.mendota.key[which(sig.mendota.key$Cluster %in% clusters),]
mendota.eigenvectors$Timepoint <- rownames(mendota.eigenvectors)
mendota.eigenvectors <- melt(mendota.eigenvectors)
plot.colors <- NA
plot.colors[which(mendota.eigenvectors$value > 0)] <- "green"
plot.colors[which(mendota.eigenvectors$value < 0)] <- "red"
mendota.eigenvectors$Sign <- plot.colors
mendota.eigenvectors$Timepoint <- factor(mendota.eigenvectors$Timepoint, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
modules <- paste("ME", clusters, sep = "")
plot.mendota.eigenvectors <- mendota.eigenvectors[which(mendota.eigenvectors$variable %in% modules), ]
ME1 <- ggplot(data = plot.mendota.eigenvectors, aes(x = Timepoint, y = variable, fill = value)) + geom_tile() + labs(y = "Cluster", x = "Time") + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow", midpoint = 0) + scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels=c("5AM", "9AM", "1PM", "5PM", "9PM", "1AM", "5AM", "9AM", "1PM", "5PM", "9PM", "1AM"))
phylum_totals <- table(plot.sig.mendota.key$Phylum)
plot.sig.mendota.key$Cluster <- paste("ME", plot.sig.mendota.key$Cluster, sep = "")
keep.phyla <- names(phylum_totals)[which(phylum_totals > 75)]
plot.sig.mendota.key <- plot.sig.mendota.key[which(plot.sig.mendota.key$Phylum %in% keep.phyla), ]
plot.sig.mendota.key$Cluster <- factor(plot.sig.mendota.key$Cluster, levels = rev(c("ME53", "ME36", "ME21", "ME3", "ME39", "ME6", "ME37", "ME7", "ME26", "ME8", "ME2", "ME55", "ME30", "ME47", "ME25", "ME13")))
ME2 <- ggplot(data = plot.sig.mendota.key, aes(y = log(Totals), x = Cluster, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = NULL, y = "Log of Total Reads") + coord_flip() + scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#33a02c", "#fb9a99", "#b2df8a", "#1f78b4", "grey", "#a6cee3"))

plot.mendota <- plot_grid(ME1, ME2, labels = c("A", "B"))


# Read and process Sparkling data
sig.spark.key <- read.csv("data/WGCNA_sparkling_results.csv", header = T)
spark.eigenvectors <- read.csv("data/WGCNA_sparkling_eigenvectors.csv", header = T, row.names = 1)
clusters <- c(4, 6, 10, 13, 14, 15, 19, 27, 31)
plot.sig.spark.key <- sig.spark.key[which(sig.spark.key$Cluster %in% clusters),]
spark.eigenvectors$Timepoint <- rownames(spark.eigenvectors)
spark.eigenvectors <- melt(spark.eigenvectors)
plot.colors <- NA
plot.colors[which(spark.eigenvectors$value > 0)] <- "green"
plot.colors[which(spark.eigenvectors$value < 0)] <- "red"
spark.eigenvectors$Sign <- plot.colors
spark.eigenvectors$Timepoint <- factor(spark.eigenvectors$Timepoint, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
modules <- paste("ME", clusters, sep = "")
plot.spark.eigenvectors <- spark.eigenvectors[which(spark.eigenvectors$variable %in% modules), ]
SP1 <- ggplot(data = plot.spark.eigenvectors, aes(x = Timepoint, y = variable, fill = value)) + geom_tile() + labs(y = "Cluster", x = "Time") + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow", midpoint = 0) + scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels=c("5AM", "9AM", "1PM", "5PM", "9PM", "1AM", "5AM", "9AM", "1PM", "5PM", "9PM", "1AM"))
phylum_totals <- table(plot.sig.spark.key$Phylum)
plot.sig.spark.key$Cluster <- paste("ME", plot.sig.spark.key$Cluster, sep = "")
keep.phyla <- names(phylum_totals)[which(phylum_totals > 50)]
plot.sig.spark.key <- plot.sig.spark.key[which(plot.sig.spark.key$Phylum %in% keep.phyla), ]
plot.sig.spark.key$Cluster <- factor(plot.sig.spark.key$Cluster, levels = rev(c("ME19", "ME13", "ME10", "ME15", "ME4", "ME31", "ME6", "ME27", "ME14")))
SP2 <- ggplot(data = plot.sig.spark.key, aes(y = log(Totals), x = Cluster, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = NULL, y = "Log of Total Reads") + coord_flip() + scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#33a02c", "#1f78b4", "grey"))

plot.spark <- plot_grid(SP1, SP2, labels = c("A", "B"))

# Read and process Mendota data
sig.trout.key <- read.csv("data/WGCNA_trout_results.csv", header = T)
trout.eigenvectors <- read.csv("data/WGCNA_trout_eigenvectors.csv", header = T, row.names = 1)
clusters <- c(1, 3, 4, 11, 16, 18)
plot.sig.trout.key <- sig.trout.key[which(sig.trout.key$Cluster %in% clusters),]
trout.eigenvectors$Timepoint <- rownames(trout.eigenvectors)
trout.eigenvectors <- melt(trout.eigenvectors)
plot.colors <- NA
plot.colors[which(trout.eigenvectors$value > 0)] <- "green"
plot.colors[which(trout.eigenvectors$value < 0)] <- "red"
trout.eigenvectors$Sign <- plot.colors
trout.eigenvectors$Timepoint <- factor(trout.eigenvectors$Timepoint, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
modules <- paste("ME", clusters, sep = "")
plot.trout.eigenvectors <- trout.eigenvectors[which(trout.eigenvectors$variable %in% modules), ]
TB1 <- ggplot(data = plot.trout.eigenvectors, aes(x = Timepoint, y = variable, fill = value)) + geom_tile() + labs(y = "Cluster", x = "Time") + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow", midpoint = 0) + scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels=c("5AM", "9AM", "1PM", "5PM", "9PM", "1AM", "5AM", "9AM", "1PM", "5PM", "9PM", "1AM"))
phylum_totals <- table(plot.sig.trout.key$Phylum)
plot.sig.trout.key$Cluster <- paste("ME", plot.sig.trout.key$Cluster, sep = "")
keep.phyla <- names(phylum_totals)[which(phylum_totals > 50)]
plot.sig.trout.key <- plot.sig.trout.key[which(plot.sig.trout.key$Phylum %in% keep.phyla), ]
plot.sig.trout.key$Cluster <- factor(plot.sig.trout.key$Cluster, levels = rev(c("ME1", "ME4", "ME3", "ME11", "ME18", "ME16")))
TB2 <- ggplot(data = plot.sig.trout.key, aes(y = log(Totals), x = Cluster, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = NULL, y = "Log of Total Reads") + coord_flip() + scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#1f78b4", "grey", "#a6cee3"))

plot.trout <- plot_grid(TB1, TB2, labels = c("A", "B"))


# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Gene Expression in Oligotrophic, Dystrophic, and Eutrophic Systems"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         radioButtons("radio", h3("Choose a lake:"),
                       choices = list("Lake Mendota" = 1, "Sparkling Lake" = 2,
                                      "Trout Bog" = 3), selected = 1),
         numericInput("num", label = h3("Choose a cluster:"), value = 1), 
         hr(),
         fluidRow(column(3, verbatimTextOutput("value"))),
         helpText("Choose a lake to see an overview of diel trends. Each row represents a cluster of genes with similar expression profiles. The heatmap shows the eigenvector of each cluster - this is a unitless measure. Yellow is higher expression and blue is lower expression. The 2nd panel shows the total number of reads in each cluster and their phylum-level assignment. Type in a cluster number to look more closely at their composition. The Eigenvector tab shows overall trend for each cluster - note that negative values are not negative expression, just lower expression. The Phyla tab shows phylum-level assignments unstacked, without unclassified reads. The Gene tab is a searchable table of the genes and their products in each cluster.")
         ),

      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(
          tabPanel("Heatmap", plotOutput("heatmap")),
          tabPanel("Eigenvector", plotOutput("eigenvector")),
          tabPanel("Phyla", plotOutput("phyla")),
          tabPanel("Genes", dataTableOutput("genes"))
        )
      )
   )
)
# Define server logic required to draw a histogram
server <- function(input, output) {
   output$heatmap <- renderPlot({
     image <- switch(input$radio, 
                     "1" = plot.mendota,
                     "2" = plot.spark,
                     "3" = plot.trout)
     
     image
   })
   output$eigenvector <- renderPlot({
     long_eigenvectors <- switch(input$radio, 
                                 "1" = mendota.eigenvectors,
                                 "2" = spark.eigenvectors,
                                 "3" = trout.eigenvectors)
     cluster <- paste("ME", input$num, sep = "")
     ggplot(data = long_eigenvectors[which(long_eigenvectors$variable == cluster), ], aes(x = Timepoint, y = value, fill = Sign)) + geom_bar(stat = "identity") + labs(title = cluster) + scale_fill_manual(values = c("yellow", "dodgerblue")) + theme(legend.position = "none")
   })
   output$phyla <- renderPlot({
     key <- switch(input$radio, 
                                 "1" = plot.sig.mendota.key,
                                 "2" = plot.sig.spark.key,
                                 "3" = plot.sig.trout.key)
     
     phylum_totals <- table(key$Phylum)
     key$Cluster <- paste("ME", input$num, sep = "")
     keep.phyla <- names(phylum_totals)[which(phylum_totals > 50)]
     key <- key[which(key$Phylum %in% keep.phyla), ]
     ggplot(data = key[which(key$Phylum != "Unclassified"),], aes(y = Totals, x = Phylum, fill = Phylum)) + geom_bar(stat = "identity") + labs(title = paste("Cluster", input$num)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20)) + theme(legend.position = "none")
     
   })
   output$genes <- renderDataTable({
     key <- switch(input$radio, 
                   "1" = sig.mendota.key,
                   "2" = sig.spark.key,
                   "3" = sig.trout.key)
     key <- key[which(key$Cluster == input$num),]
     key
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

