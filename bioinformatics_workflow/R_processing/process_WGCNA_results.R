# Process WGCNA results

# Packages
library(ggplot2)
library(cowplot)
library(reshape2)

##### Mendota
# Read data
sig.mendota.key <- read.csv("D:/geodes_data_tables/WGCNA_mendota_results.csv", header = T)
eigenvectors <- read.csv("D:/geodes_data_tables/WGCNA_mendota_eigenvectors.csv", header = T, row.names = 1)

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

# Which clusters to plot? Visualize and choose
eigenvectors$Timepoint <- rownames(eigenvectors)
long_eigenvectors <- melt(eigenvectors)
plot.colors <- NA
plot.colors[which(long_eigenvectors$value > 0)] <- "green"
plot.colors[which(long_eigenvectors$value < 0)] <- "red"
long_eigenvectors$Sign <- plot.colors
long_eigenvectors$Timepoint <- factor(long_eigenvectors$Timepoint, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
for(i in 0:57){
  cluster = paste("ME", i, sep = "")
  p <- ggplot(data = long_eigenvectors[which(long_eigenvectors$variable == cluster), ], aes(x = Timepoint, y = value, fill = Sign)) + geom_bar(stat = "identity") + labs(title = cluster) + scale_fill_manual(values = c("green", "red")) + theme(legend.position = "none")
  print(p)
}

clusters <- c(2, 3, 6, 7, 8, 13, 21, 25, 26, 30, 36, 37, 39, 47, 53, 55)
plot.sig.mendota.key <- sig.mendota.key[which(sig.mendota.key$Cluster %in% clusters),]

# Panel A
modules <- paste("ME", clusters, sep = "")
plot.long_eigenvectors <- long_eigenvectors[which(long_eigenvectors$variable %in% modules), ]
ME1 <- ggplot(data = plot.long_eigenvectors, aes(x = Timepoint, y = variable, fill = value)) + geom_tile() + labs(y = "Cluster", x = "Time") + scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0) + scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels=c("5AM", "9AM", "1PM", "5PM", "9PM", "1AM", "5AM", "9AM", "1PM", "5PM", "9PM", "1AM")) + theme(legend.position = "none")

# Panel B: Phylum barcharts
phylum_totals <- table(plot.sig.mendota.key$Phylum)
plot.sig.mendota.key$Cluster <- paste("ME", plot.sig.mendota.key$Cluster, sep = "")
keep.phyla <- names(phylum_totals)[which(phylum_totals > 75)]
plot.sig.mendota.key <- plot.sig.mendota.key[which(plot.sig.mendota.key$Phylum %in% keep.phyla), ]
plot.sig.mendota.key$Cluster <- factor(plot.sig.mendota.key$Cluster, levels = rev(c("ME53", "ME36", "ME21", "ME3", "ME39", "ME6", "ME37", "ME7", "ME26", "ME8", "ME2", "ME55", "ME30", "ME47", "ME25", "ME13")))
ME2 <- ggplot(data = plot.sig.mendota.key, aes(y = log(Totals), x = Cluster, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = NULL, y = "Log of Total Reads") + coord_flip() + scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#33a02c", "#fb9a99", "#b2df8a", "#1f78b4", "grey", "#a6cee3"))

plot_grid(ME1, ME2, labels = c("A", "B"))


# Get genes from clusters

x <- sig.mendota.key[which(sig.mendota.key$Cluster == 13),]
x <- x[order(x$Totals),]
x[(dim(x)[1] - 50): dim(x)[1],]

# Plot phyla, not including unclassified
phyla_breakdown <- aggregate(Totals ~ Phylum, x, sum)
ggplot(phyla_breakdown[grep("Unclassified", phyla_breakdown$Phylum, invert = T), ], aes(x = Phylum, y = Totals)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#### Sparkling
# Read data
sig.spark.key <- read.csv("D:/geodes_data_tables/WGCNA_sparkling_results.csv", header = T)
eigenvectors <- read.csv("D:/geodes_data_tables/WGCNA_sparkling_eigenvectors.csv", header = T, row.names = 1)

# Fix taxonomy
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

# Which clusters to plot? Visualize and choose
eigenvectors$Timepoint <- rownames(eigenvectors)
long_eigenvectors <- melt(eigenvectors)
plot.colors <- NA
plot.colors[which(long_eigenvectors$value > 0)] <- "green"
plot.colors[which(long_eigenvectors$value < 0)] <- "red"
long_eigenvectors$Sign <- plot.colors
long_eigenvectors$Timepoint <- factor(long_eigenvectors$Timepoint, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
for(i in 0:46){
  cluster = paste("ME", i, sep = "")
  p <- ggplot(data = long_eigenvectors[which(long_eigenvectors$variable == cluster), ], aes(x = Timepoint, y = value, fill = Sign)) + geom_bar(stat = "identity") + labs(title = cluster) + scale_fill_manual(values = c("green", "red")) + theme(legend.position = "none")
  print(p)
}

clusters <- c(4, 6, 10, 13, 14, 15, 19, 27, 31)
plot.sig.spark.key <- sig.spark.key[which(sig.spark.key$Cluster %in% clusters),]

# Panel A
modules <- paste("ME", clusters, sep = "")
plot.long_eigenvectors <- long_eigenvectors[which(long_eigenvectors$variable %in% modules), ]
SP1 <- ggplot(data = plot.long_eigenvectors, aes(x = Timepoint, y = variable, fill = value)) + geom_tile() + labs(y = "Cluster", x = "Time") + scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0) + scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), labels=c("5AM", "9AM", "1PM", "5PM", "9PM", "1AM", "5AM", "9AM", "1PM", "5PM", "9PM", "1AM")) + theme(legend.position = "none")

# Panel B: Phylum barcharts
phylum_totals <- table(plot.sig.spark.key$Phylum)
plot.sig.spark.key$Cluster <- paste("ME", plot.sig.spark.key$Cluster, sep = "")
keep.phyla <- names(phylum_totals)[which(phylum_totals > 50)]
plot.sig.spark.key <- plot.sig.spark.key[which(plot.sig.spark.key$Phylum %in% keep.phyla), ]
plot.sig.spark.key$Cluster <- factor(plot.sig.spark.key$Cluster, levels = rev(c("ME19", "ME13", "ME10", "ME15", "ME4", "ME31", "ME6", "ME27", "ME14")))
SP2 <- ggplot(data = plot.sig.spark.key, aes(y = log(Totals), x = Cluster, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = NULL, y = "Log of Total Reads") + coord_flip() + scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#33a02c", "#1f78b4", "grey"))

plot_grid(SP1, SP2, labels = c("A", "B"))


# Get genes from clusters

x <- sig.spark.key[which(sig.spark.key$Cluster == 14),]
x <- x[order(x$Totals),]
x[(dim(x)[1] - 50): dim(x)[1],]

# Plot phyla, not including unclassified
phyla_breakdown <- aggregate(Totals ~ Phylum, x, sum)
ggplot(phyla_breakdown[grep("Unclassified", phyla_breakdown$Phylum, invert = T), ], aes(x = Phylum, y = Totals)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#### Trout
# Read data
sig.trout.key <- read.csv("D:/geodes_data_tables/WGCNA_trout_results.csv", header = T)
eigenvectors <- read.csv("D:/geodes_data_tables/WGCNA_trout_eigenvectors.csv", header = T, row.names = 1)

# Fix taxonomy
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

# Which clusters to plot? Visualize and choose
eigenvectors$Timepoint <- rownames(eigenvectors)
long_eigenvectors <- melt(eigenvectors)
plot.colors <- NA
plot.colors[which(long_eigenvectors$value > 0)] <- "green"
plot.colors[which(long_eigenvectors$value < 0)] <- "red"
long_eigenvectors$Sign <- plot.colors
long_eigenvectors$Timepoint <- factor(long_eigenvectors$Timepoint, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
for(i in 0:24){
  cluster = paste("ME", i, sep = "")
  p <- ggplot(data = long_eigenvectors[which(long_eigenvectors$variable == cluster), ], aes(x = Timepoint, y = value, fill = Sign)) + geom_bar(stat = "identity") + labs(title = cluster) + scale_fill_manual(values = c("green", "red")) + theme(legend.position = "none")
  print(p)
}

clusters <- c(1, 3, 4, 11, 16, 18)
plot.sig.trout.key <- sig.trout.key[which(sig.trout.key$Cluster %in% clusters),]

# Panel A
modules <- paste("ME", clusters, sep = "")
plot.long_eigenvectors <- long_eigenvectors[which(long_eigenvectors$variable %in% modules), ]
TB1 <- ggplot(data = plot.long_eigenvectors, aes(x = Timepoint, y = variable, fill = value)) + geom_tile() + labs(y = "Cluster", x = "Time") + scale_fill_gradient2(low = "red", mid = "white", high = "green", midpoint = 0) + scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7"), labels=c("5AM", "9AM", "1PM", "5PM", "9PM", "1AM", "5AM")) + theme(legend.position = "none")

# Panel B: Phylum barcharts
phylum_totals <- table(plot.sig.trout.key$Phylum)
plot.sig.trout.key$Cluster <- paste("ME", plot.sig.trout.key$Cluster, sep = "")
keep.phyla <- names(phylum_totals)[which(phylum_totals > 50)]
plot.sig.trout.key <- plot.sig.trout.key[which(plot.sig.trout.key$Phylum %in% keep.phyla), ]
plot.sig.trout.key$Cluster <- factor(plot.sig.trout.key$Cluster, levels = rev(c("ME1", "ME4", "ME3", "ME11", "ME18", "ME16")))
TB2 <- ggplot(data = plot.sig.trout.key, aes(y = log(Totals), x = Cluster, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = NULL, y = "Log of Total Reads") + coord_flip() + scale_fill_manual(values = c("#fdbf6f", "#e31a1c", "#1f78b4", "grey", "#a6cee3"))

plot_grid(TB1, TB2, labels = c("A", "B"))


# Get genes from clusters

x <- sig.trout.key[which(sig.trout.key$Cluster == 16),]
x <- x[order(x$Totals),]
x[(dim(x)[1] - 50): dim(x)[1],]

# Plot phyla, not including unclassified
phyla_breakdown <- aggregate(Totals ~ Phylum, x, sum)
ggplot(phyla_breakdown[grep("Unclassified", phyla_breakdown$Phylum, invert = T), ], aes(x = Phylum, y = Totals)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
