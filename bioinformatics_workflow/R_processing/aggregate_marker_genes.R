# Marker gene analysis - look mom, no clustering!

library(ggplot2)
library(cowplot)
library(reshape2)
library(igraph)
library(RColorBrewer)

zscore <- function(counts){
  z <- (counts - mean(counts)) / sd(counts)
  return(z)
}

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=18, face="bold")
  )

### Load data (start with only one to save RAM and comment the rest out)
# Normalized read tables
snorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
snorm <- snorm[, which(colnames(snorm) != "GEODES014.nonrRNA" & colnames(snorm) != "GEODES033.nonrRNA")]
# tnorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
# tnorm <- tnorm[, which(colnames(tnorm) != "GEODES065.nonrRNA")]
#mnorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)
#mnorm <- mnorm[, which(colnames(mnorm) != "GEODES158.nonrRNA")]

# Gene keys
#mendota_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
spark_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
# trout_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# Choose your search term:
#searchterm <- c("Rubisco|RuBisCO|rubisco|bisphosphate carboxylase")
#searchterm <- c("citrate lyase|Citrate lyase")
#searchterm <- c("nitrogenase|Nitrogenase|NifH|NifD|NifK")
#searchterm <- c("Chitobiase|chitobiase|chitinase|Chitinase")
#searchterm <- c("glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase")
#searchterm <- c("sulfite reductase|Sulfite reductase")
#searchterm <- c("glucosaminidase|hexosaminidase")
#searchterm <- c("putrescine|spermidine|polymamine|Putrescine|Spermidine|Polyamine")
#searchterm <- c("cobalamin|thiamin|pyridoxal")
#searchterm <- c("rhamnulose-1-phosphate|Rhamnulose-1-phosphate|fuculose-phosphate|Fuculose phosphate|ribulokinase|mannose-6-phosphate")
#searchterm <- c("aminopeptidase|Aminopeptidase")
#searchterm <- c("alkaline phosphatase|Alkaline phosphatase")
#searchterm <- c("sugar transport|ose transport")
#searchterm <- c("rhodopsin|Rhodopsin")
#searchterm <- c("phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene")
#searchterm <- c("photosynth|Photosynth")
searchterm <- c("cellulase|cellulose")
#searchterm <- c("sulfatase")
#searchterm <- c("sulfite reductase|Sulfite reductase")
#searchterm <- c("peroxidase|peroxide|catalase")
#searchterm <- c("ammonia monooxygenase|hydroxylamine oxidoreductase")
#searchterm <- c("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase")
#searchterm <- c("nitrate reductase|Nitrate reductase|nitrate oxidoreductase|Nitrate oxidoreductase")
#searchterm <- c("urease")
#searchterm <- c("protease")

marker_genes <- spark_key[grep(searchterm, spark_key$Product), ]
marker_genes$Taxonomy <- gsub("Bacteria;", "", marker_genes$Taxonomy)
marker_genes$Taxonomy <- gsub("Eukaryota;", "", marker_genes$Taxonomy)
marker_genes$Taxonomy <- gsub("Proteobacteria;", "", marker_genes$Taxonomy)
marker_genes$Phylum <- sapply(strsplit(as.character(marker_genes$Taxonomy),";"), `[`, 1)
marker_genes$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", marker_genes$Phylum)
marker_genes$Phylum[grep("NO CLASS|Blank|unclassified|None", marker_genes$Phylum)] <- "Unclassified"
marker_genes$Phylum[which(is.na(marker_genes$Phylum) == T)] <- "Unclassified"
marker_genes <- marker_genes[which(marker_genes$Phylum != "Unclassified"), ]

marker_genes_snorm <- snorm[match(marker_genes$Gene, rownames(snorm)), ]
marker_genes_snorm <- marker_genes_snorm[which(rowSums(marker_genes_snorm) > (dim(marker_genes_snorm)[2] * 1000)), ]

#Aggregate by timepoint
marker_genes_snorm$Genes <- rownames(marker_genes_snorm)
marker_genes_snorm <- melt(marker_genes_snorm)
marker_genes_snorm$variable <- gsub(".nonrRNA", "", marker_genes_snorm$variable)
marker_genes_snorm$Timepoint <- metadata$Timepoint[match(marker_genes_snorm$variable, metadata$Sample)]
agg_marker_genes_snorm <- aggregate(value ~ Genes + Timepoint, data = marker_genes_snorm, mean)
new_marker_genes_snorm <- reshape(agg_marker_genes_snorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(new_marker_genes_snorm) <- new_marker_genes_snorm[, 1]
new_marker_genes_snorm <- new_marker_genes_snorm[, 2:dim(new_marker_genes_snorm)[2]]

# # Add taxonomy information and sum by unique taxonomy
# new_marker_genes_snorm$Genes <- rownames(new_marker_genes_snorm)
# new_marker_genes_snorm$Taxonomy <- marker_genes$Taxonomy[match(new_marker_genes_snorm$Genes, marker_genes$Gene)]
# new_marker_genes_snorm <- melt(new_marker_genes_snorm)
# new_marker_genes_snorm$Taxonomy <- gsub("Bacteria;", "", new_marker_genes_snorm$Taxonomy)
# new_marker_genes_snorm$Taxonomy <- gsub("Eukaryota;", "", new_marker_genes_snorm$Taxonomy)
# new_marker_genes_snorm$Taxonomy <- gsub("Proteobacteria;", "", new_marker_genes_snorm$Taxonomy)
# new_marker_genes_snorm$Phylum <- sapply(strsplit(as.character(new_marker_genes_snorm$Taxonomy),";"), `[`, 1)
# new_marker_genes_snorm$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", new_marker_genes_snorm$Phylum)
# #summed_markers <- aggregate(value ~ variable + Taxonomy, new_marker_genes_snorm, sum)
# summed_markers <- aggregate(value ~ variable + Phylum, new_marker_genes_snorm, sum)
# summed_markers <- summed_markers[grep("unclassified|NO CLASS|None|Blank", summed_markers$Phylum, invert = T), ]
# 
# # Plot total sums by phyla
# summed_phyla <- aggregate(value ~ Phylum, new_marker_genes_snorm, sum)
# summed_phyla <- summed_phyla[grep("unclassified|NO CLASS|None|Blank", summed_phyla$Phylum, invert = T), ]
# summed_phyla <- summed_phyla[which(summed_phyla$value > 1000000),]
# p2 <- ggplot(data = summed_phyla, aes(x = Phylum, y = value)) + geom_bar(stat = "identity") + coord_flip() + labs(x = NULL, y = "Summed Read Counts", title = "Marker Gene Reads Assigned")
# 
# summed_markers <- summed_markers[which(summed_markers$Phylum %in% summed_phyla$Phylum),]
# 
# #Change to fold change
# taxa <- unique(summed_markers$Phylum)
# for(i in 1:length(taxa)){
#   subset <- summed_markers[which(summed_markers$Phylum == taxa[i]), ]
#   ave <- mean(subset$value)
#   subset$value <- subset$value - ave
#   maxrange <- max(abs(subset$value))
#   subset$value <- subset$value * (10/maxrange)
#   summed_markers[which(summed_markers$Phylum == taxa[i]), ] <- subset
# }
# 
# p1 <- ggplot(data = summed_markers, aes(x = variable, y = Phylum, fill = value)) + geom_tile()+ scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM")) + theme(legend.position = "none") + labs(x = NULL, y = NULL, title = paste(searchterm))
# 
# full_plot <- plot_grid(p1, p2)
# full_plot
# 
# # correlation matrix
# phyla <- c("Cyanobacteria", "Actinobacteria", "Bacteroidetes")
# cor.matrix <- matrix(nrow = length(phyla), ncol = length(phyla), 0)
# rownames(cor.matrix) <- colnames(cor.matrix) <- phyla
# for(i in 1:length(phyla)){
#   thing1 <- summed_markers[which(summed_markers$Phylum == phyla[i]), ]
#     for(j in 1:length(phyla)){
#       thing2 <- summed_markers[which(summed_markers$Phylum == phyla[j]), ]
#       cor.matrix[i, j] <- cor(thing1$value, thing2$value)
#     }
# }
# cor.matrix
# 
# 
# # Don't aggregate by taxonomy
# new_marker_genes_snorm$Genes <- rownames(new_marker_genes_snorm)
# melted_markers <- melt(new_marker_genes_snorm)
# genes <- marker_genes$Gene
# for(i in 1:length(genes)){
#   subset <- melted_markers[which(melted_markers$Genes == genes[i]), ]
#   ave <- mean(subset$value)
#   subset$value <- subset$value - ave
#   maxrange <- max(abs(subset$value))
#   subset$value <- subset$value * (10/maxrange)
#   melted_markers[which(melted_markers$Genes == genes[i]), ] <- subset
# }
# 
# 
# cor.matrix <- matrix(nrow = length(genes), ncol = length(genes), 0)
# rownames(cor.matrix) <- colnames(cor.matrix) <- genes
# for(i in 1:length(genes)){
#   thing1 <- melted_markers[which(melted_markers$Genes == genes[i]), ]
#   for(j in 1:length(genes)){
#     thing2 <- melted_markers[which(melted_markers$Genes == genes[j]), ]
#     cor.matrix[i, j] <- cor(thing1$value, thing2$value)
#   }
# }


# cor.matrix2 <- cor.matrix
# cor.matrix2[which(cor.matrix2 < 0.8 & cor.matrix2 > -0.8)] <- 0
# net <- graph_from_incidence_matrix(cor.matrix2)
# l <- layout_on_sphere(net)
# plot(net, vertex.label = NA, vertex.size = 5)
# # Identify hubs? Not bother with network? Abundance threshold?

#Keep in long format
cor.data <- expand.grid(rownames(new_marker_genes_snorm), rownames(new_marker_genes_snorm))
cor.data <- cor.data[which((cor.data$Var1 == cor.data$Var2) == F),]
test.dups <- duplicated(t(apply(cor.data, 1, sort)))
cor.data <- cor.data[which(test.dups == F), ]

correlations <- c()
for(i in 1:dim(cor.data)[1]){
  thing1 <- new_marker_genes_snorm[which(rownames(new_marker_genes_snorm) == cor.data$Var1[i]), ]
  thing2 <- new_marker_genes_snorm[which(rownames(new_marker_genes_snorm) == cor.data$Var2[i]), ]
  correlations[i] <- cor(as.numeric(thing1), as.numeric(thing2))
  if(abs(correlations[i]) < 0.6){
    correlations[i] <- 0
  }
}
cor.data$edges <- correlations
cor.data <- cor.data[which(cor.data$edges != 0), ]
kept_genes <- unique(c(as.character(cor.data$Var1), as.character(cor.data$Var2)))
vert.data <- marker_genes[match(kept_genes, marker_genes$Gene),]

vert.data$gene_expression <- rowSums(new_marker_genes_snorm)[match(vert.data$Gene, rownames(new_marker_genes_snorm))]


phylum_colors <- brewer.pal(length(unique(vert.data$Phylum)), "Paired")
net <- graph_from_data_frame(cor.data, directed = F, vertices = vert.data)
V(net)$color <- phylum_colors[factor(V(net)$Phylum, levels = unique(V(net)$Phylum))]
V(net)$size <- log(vert.data$gene_expression)/2
l <- layout_with_lgl(net)
plot(net, vertex.label = NA, vertex.size = V(net)$size, layout = l, vertex.color=V(net)$color)
legend(x=-1.5, y=0, unique(V(net)$Phylum), pch=21, col="#777777", pt.bg=phylum_colors, pt.cex=2, cex=.8, bty="n", ncol=1)

ceb <- cluster_edge_betweenness(net) 
# dendPlot(ceb, mode="hclust")
# plot(ceb, net, vertex.label = NA, vertex.size = V(net)$size, layout = l)
# legend(x=-1.5, y=0, unique(V(net)$Phylum), pch=21, col="#777777", pt.bg=phylum_colors, pt.cex=2, cex=.8, bty="n", ncol=1) 
cluster.nums <- membership(ceb)
vert.data$Cluster <- cluster.nums[match(vert.data$Gene, names(cluster.nums))]
vert.data <- vert.data[order(vert.data$Cluster), ]

cluster.reads <- aggregate(gene_expression ~ Cluster, vert.data, sum)
cluster.reads <- cluster.reads[order(cluster.reads$gene_expression), ]
cluster.reads
# tol14rainbow=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
# cluster.colors <- tol14rainbow[factor(vert.data$Cluster, levels = unique(vert.data$Cluster))]
# plot(net, vertex.label = NA, vertex.size = V(net)$size, layout = l, vertex.color=cluster.colors)
# legend(x=-1.5, y=0, unique(vert.data$Cluster), pch=21, col="#777777", pt.bg=tol14rainbow, pt.cex=2, cex=.8, bty="n", ncol=1)

# Plot pie chart of comp and trend over time in selected clusters

c <- 4

cluster.table <- table(vert.data$Phylum[which(vert.data$Cluster == c)])
cluster.table <- data.frame(names(cluster.table), as.numeric(cluster.table))
colnames(cluster.table) <- c("Phylum", "Genes_Contributed")
ggplot(data = cluster.table, aes(x = "", y = Genes_Contributed, fill = Phylum)) + geom_bar(stat = "identity", width = 1, color = "black") + coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values = phylum_colors, name = NULL)

cluster.abun <- colSums(new_marker_genes_snorm[match(vert.data$Gene[which(vert.data$Cluster == c)], rownames(new_marker_genes_snorm)), ])
cluster.abun <- data.frame(abundance = as.numeric(cluster.abun), time = names(cluster.abun))
cluster.abun$time <- factor(cluster.abun$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
cluster.abun$abundance <- zscore(cluster.abun$abundance )
ggplot(data = cluster.abun, aes(x = time, y = abundance)) + geom_bar(stat = "identity", fill = "#5289C7") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM")) + labs(x = NULL, y = NULL, title = NULL)

# Also do this for all data
cluster.table <- table(vert.data$Phylum)
cluster.table <- data.frame(names(cluster.table), as.numeric(cluster.table))
colnames(cluster.table) <- c("Phylum", "Genes_Contributed")
ggplot(data = cluster.table, aes(x = "", y = Genes_Contributed, fill = Phylum)) + geom_bar(stat = "identity", width = 1, color = "black") + coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values = phylum_colors, name = NULL)

#vert.data <- vert.data[which(vert.data$Phylum == "Betaproteobacteria"), ]
cluster.abun <- colSums(new_marker_genes_snorm[match(vert.data$Gene, rownames(new_marker_genes_snorm)), ])
cluster.abun <- data.frame(abundance = as.numeric(cluster.abun), time = names(cluster.abun))
cluster.abun$abundance <- zscore(cluster.abun$abundance )
cluster.abun$time <- factor(cluster.abun$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
ggplot(data = cluster.abun, aes(x = time, y = abundance)) + geom_bar(stat = "identity", fill = "grey", color = "black") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM")) + labs(x = NULL, y = NULL, title = NULL)

# For when clustering breaks:
phylum_colors <- brewer.pal(12, "Paired")
cluster.table <- table(marker_genes$Phylum)
cluster.table <- data.frame(names(cluster.table), as.numeric(cluster.table))
colnames(cluster.table) <- c("Phylum", "Genes_Contributed")
ggplot(data = cluster.table, aes(x = "", y = Genes_Contributed, fill = Phylum)) + geom_bar(stat = "identity", width = 1, color = "black") + coord_polar("y", start=0) + blank_theme + theme(axis.text.x=element_blank()) + scale_fill_manual(values = phylum_colors, name = NULL)

#vert.data <- vert.data[which(vert.data$Phylum == "Betaproteobacteria"), ]
cluster.abun <- colSums(new_marker_genes_snorm)
cluster.abun <- data.frame(abundance = as.numeric(cluster.abun), time = names(cluster.abun))
cluster.abun$abundance <- zscore(cluster.abun$abundance )
cluster.abun$time <- factor(cluster.abun$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
ggplot(data = cluster.abun, aes(x = time, y = abundance)) + geom_bar(stat = "identity", fill = "grey", color = "black") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM")) + labs(x = NULL, y = NULL, title = NULL)
