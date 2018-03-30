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
sig.mendota.key <- read.csv("data/WGCNA_mendota_results_nophoto_2018-03-30.csv", header = T)
mendota.eigenvectors <- read.csv("data/WGCNA_mendota_eigenvectors_nophoto_2018-03-30.csv", header = T, row.names = 1)

# Fix taxonomy
sig.mendota.key$Taxonomy <- gsub("Bacteria;", "", sig.mendota.key$Taxonomy)
sig.mendota.key$Taxonomy <- gsub("Eukaryota;", "", sig.mendota.key$Taxonomy)
sig.mendota.key$Phylum <- sapply(strsplit(as.character(sig.mendota.key$Taxonomy),";"), `[`, 1)
sig.mendota.key$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("None", "Unclassified", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("unclassified unclassified", "Unclassified", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("unclassified", "Unclassified", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("Unclassified ", "Unclassified", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", sig.mendota.key$Phylum)
sig.mendota.key$Phylum <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", sig.mendota.key$Phylum)


clusters <- c(2, 4, 5, 6, 9, 11, 14, 15, 16, 18, 27, 28, 29, 31, 33)
plot.sig.mendota.key <- sig.mendota.key[which(sig.mendota.key$Cluster %in% clusters),]
mendota.eigenvectors$Timepoint <- rownames(mendota.eigenvectors)
mendota.eigenvectors <- melt(mendota.eigenvectors)
plot.colors <- NA
plot.colors[which(mendota.eigenvectors$value > 0)] <- "dodgerblue"
plot.colors[which(mendota.eigenvectors$value < 0)] <- "yellow"
mendota.eigenvectors$Sign <- plot.colors
mendota.eigenvectors$Timepoint <- factor(mendota.eigenvectors$Timepoint, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
modules <- paste("ME", clusters, sep = "")
plot.mendota.eigenvectors <- mendota.eigenvectors[which(mendota.eigenvectors$variable %in% modules), ]
plot.mendota.eigenvectors$variable <- gsub("ME", "Cluster", plot.mendota.eigenvectors$variable)
plot.mendota.eigenvectors$variable <- factor(plot.mendota.eigenvectors$variable, levels = c("Cluster15", "Cluster28", "Cluster9", "Cluster11", "Cluster4", "Cluster14", "Cluster18", "Cluster5", "Cluster6", "Cluster16", "Cluster2", "Cluster31", "Cluster27", "Cluster29", "Cluster33"))
ME1 <- ggplot(data = plot.mendota.eigenvectors, aes(x = Timepoint, y = variable, fill = value)) + geom_tile() + labs(y = "Cluster", x = "Time") + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow", midpoint = 0) + scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels=c("5AM", "9AM", "1PM", "5PM", "9PM", "1AM", "5AM", "9AM", "1PM", "5PM", "9PM", "1AM"))
phylum_totals <- table(plot.sig.mendota.key$Phylum)
plot.sig.mendota.key$Cluster <- paste("Cluster", plot.sig.mendota.key$Cluster, sep = "")
keep.phyla <- names(phylum_totals)
plot.sig.mendota.key <- plot.sig.mendota.key[which(plot.sig.mendota.key$Phylum %in% keep.phyla), ]
plot.sig.mendota.key$Cluster <- factor(plot.sig.mendota.key$Cluster, levels = c("Cluster15", "Cluster28", "Cluster9", "Cluster11", "Cluster4", "Cluster14", "Cluster18", "Cluster5", "Cluster6", "Cluster16", "Cluster2", "Cluster31", "Cluster27", "Cluster29", "Cluster33"))
ME2 <- ggplot(data = plot.sig.mendota.key, aes(y = log(Totals), x = Cluster, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = NULL, y = "Log of Total Reads") + coord_flip()
  
#scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#33a02c", "#fb9a99", "#b2df8a", "#1f78b4", "grey", "#a6cee3"))

plot.mendota <- plot_grid(ME1, ME2, labels = c("A", "B"))


# Read and process Sparkling data
sig.spark.key <- read.csv("data/WGCNA_sparkling_results_nophoto_2018-03-30.csv", header = T)
spark.eigenvectors <- read.csv("data/WGCNA_sparkling_eigenvectors_nophoto_2018-03-30.csv", header = T, row.names = 1)

sig.spark.key$Taxonomy <- gsub("Bacteria;", "", sig.spark.key$Taxonomy)
sig.spark.key$Taxonomy <- gsub("Eukaryota;", "", sig.spark.key$Taxonomy)
sig.spark.key$Phylum <- sapply(strsplit(as.character(sig.spark.key$Taxonomy),";"), `[`, 1)
sig.spark.key$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("None", "Unclassified", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("unclassified unclassified", "Unclassified", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("unclassified", "Unclassified", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("Unclassified ", "Unclassified", sig.spark.key$Phylum)
sig.spark.key$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", sig.spark.key$Phylum)


clusters <- c(4, 6, 7, 9, 10, 11, 12, 16)
plot.sig.spark.key <- sig.spark.key[which(sig.spark.key$Cluster %in% clusters),]
spark.eigenvectors$Timepoint <- rownames(spark.eigenvectors)
spark.eigenvectors <- melt(spark.eigenvectors)
plot.colors <- NA
plot.colors[which(spark.eigenvectors$value > 0)] <- "dodgerblue"
plot.colors[which(spark.eigenvectors$value < 0)] <- "yellow"
spark.eigenvectors$Sign <- plot.colors
spark.eigenvectors$Timepoint <- factor(spark.eigenvectors$Timepoint, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
modules <- paste("ME", clusters, sep = "")
plot.spark.eigenvectors <- spark.eigenvectors[which(spark.eigenvectors$variable %in% modules), ]
plot.spark.eigenvectors$variable <- gsub("ME", "Cluster", plot.spark.eigenvectors$variable)
SP1 <- ggplot(data = plot.spark.eigenvectors, aes(x = Timepoint, y = variable, fill = value)) + geom_tile() + labs(y = "Cluster", x = "Time") + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow", midpoint = 0) + scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels=c("5AM", "9AM", "1PM", "5PM", "9PM", "1AM", "5AM", "9AM", "1PM", "5PM", "9PM", "1AM"))
phylum_totals <- table(plot.sig.spark.key$Phylum)
keep.phyla <- names(phylum_totals)
plot.sig.spark.key$Cluster <- paste("Cluster", plot.sig.spark.key$Cluster, sep = "")
plot.sig.spark.key <- plot.sig.spark.key[which(plot.sig.spark.key$Phylum %in% keep.phyla), ]
plot.sig.spark.key$Cluster <- factor(plot.sig.spark.key$Cluster, levels = rev(c("Cluster9", "Cluster7", "Cluster6", "Cluster4", "Cluster16", "Cluster12", "Cluster11", "Cluster10")))
SP2 <- ggplot(data = plot.sig.spark.key, aes(y = log(Totals), x = Cluster, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = NULL, y = "Log of Total Reads") + coord_flip()

#scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#33a02c", "#1f78b4", "grey"))

plot.spark <- plot_grid(SP1, SP2, labels = c("A", "B"))

# Read and process Trout Bog data
sig.trout.key <- read.csv("data/WGCNA_trout_results_nophoto_2018-03-30.csv", header = T)
trout.eigenvectors <- read.csv("data/WGCNA_trout_eigenvectors_nophoto_2018-03-30.csv", header = T, row.names = 1)

sig.trout.key$Taxonomy <- gsub("Bacteria;", "", sig.trout.key$Taxonomy)
sig.trout.key$Taxonomy <- gsub("Eukaryota;", "", sig.trout.key$Taxonomy)
sig.trout.key$Phylum <- sapply(strsplit(as.character(sig.trout.key$Taxonomy),";"), `[`, 1)
sig.trout.key$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("None", "Unclassified", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("unclassified unclassified", "Unclassified", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("unclassified", "Unclassified", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("Unclassified ", "Unclassified", sig.trout.key$Phylum)
sig.trout.key$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", sig.trout.key$Phylum)

clusters <- c(0, 2, 3, 10, 15)
plot.sig.trout.key <- sig.trout.key[which(sig.trout.key$Cluster %in% clusters),]
trout.eigenvectors$Timepoint <- rownames(trout.eigenvectors)
trout.eigenvectors <- melt(trout.eigenvectors)
plot.colors <- NA
plot.colors[which(trout.eigenvectors$value > 0)] <- "dodgerblue"
plot.colors[which(trout.eigenvectors$value < 0)] <- "yellow"
trout.eigenvectors$Sign <- plot.colors
trout.eigenvectors$Timepoint <- factor(trout.eigenvectors$Timepoint, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
modules <- paste("ME", clusters, sep = "")
plot.trout.eigenvectors <- trout.eigenvectors[which(trout.eigenvectors$variable %in% modules), ]
plot.trout.eigenvectors$variable <- gsub("ME", "Cluster", plot.trout.eigenvectors$variable)
TB1 <- ggplot(data = plot.trout.eigenvectors, aes(x = Timepoint, y = variable, fill = value)) + geom_tile() + labs(y = "Cluster", x = "Time") + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow", midpoint = 0) + scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels=c("5AM", "9AM", "1PM", "5PM", "9PM", "1AM", "5AM", "9AM", "1PM", "5PM", "9PM", "1AM"))
phylum_totals <- table(plot.sig.trout.key$Phylum)
plot.sig.trout.key$Cluster <- paste("Cluster", plot.sig.trout.key$Cluster, sep = "")
keep.phyla <- names(phylum_totals)
plot.sig.trout.key <- plot.sig.trout.key[which(plot.sig.trout.key$Phylum %in% keep.phyla), ]
plot.sig.trout.key$Cluster <- factor(plot.sig.trout.key$Cluster, levels = rev(c("Cluster3", "Cluster2", "Cluster15", "Cluster10", "Cluster0")))
TB2 <- ggplot(data = plot.sig.trout.key, aes(y = log(Totals), x = Cluster, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = NULL, y = "Log of Total Reads") + coord_flip() 
#scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#1f78b4", "grey", "#a6cee3"))

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

