# Plot higher resolution taxonomic group aggregated trends over the time series
library(ggplot2)
library(cowplot)
library(reshape2)
library(igraph)
zscore <- function(counts){
  z <- (counts - mean(counts)) / sd(counts)
  return(z)
}

### Load data (start with only one to save RAM and comment the rest out)
# Normalized read tables
snorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
tnorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mnorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

# Remove abnormally large metatranscriptomes
snorm <- snorm[, which(colnames(snorm) != "GEODES014.nonrRNA" & colnames(snorm) != "GEODES033.nonrRNA")]
tnorm <- tnorm[, which(colnames(tnorm) != "GEODES065.nonrRNA")]
mnorm <- mnorm[, which(colnames(mnorm) != "GEODES158.nonrRNA")]

# Gene keys
mendota_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
spark_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
trout_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

mendota_key$Taxonomy <- gsub("Bacteria;", "", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("Eukaryota;", "", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("NO CLASSIFICATION MH", "Unclassified", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("NO CLASSIFICATION LP", "Unclassified", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("None", "Unclassified", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("unclassified unclassified", "Unclassified", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("unclassified Oligohymenophorea", "Ciliophora", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("unclassified Pelagophyceae", "Ochrophyta", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("unclassified", "Unclassified", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("Unclassified ", "Unclassified", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("UnclassifiedIsochrysidales", "Haptophyta", mendota_key$Taxonomy)
mendota_key$Taxonomy[grep("Blank", mendota_key$Taxonomy)] <- "Unclassified"
mendota_key <- mendota_key[which(mendota_key$Taxonomy != "Unclassified" & mendota_key$Taxonomy != "internal standard" & mendota_key$Taxonomy != ""),]

keep <- match(rownames(mnorm), mendota_key$Gene)
mnorm <- mnorm[which(is.na(keep) == F),]

mnorm$Genes <- rownames(mnorm)
mnorm <- melt(mnorm)
mnorm$variable <- gsub(".nonrRNA", "", mnorm$variable)
mnorm$Time <- metadata$Timepoint[match(mnorm$variable, metadata$Sample)]
mnorm$Taxonomy <- mendota_key$Taxonomy[match(mnorm$Genes, mendota_key$Gene)]

agg_mnorm <- aggregate(value ~ Time + Taxonomy, mnorm, mean)

agg_mnorm$Time <- factor(agg_mnorm$Time, levels = c("0", "4", "8", "12", "16", "20", "24", "28", "32", "36", "40", "44"))

# convert to wide format for clustering
wide_agg_mnorm <- reshape(agg_mnorm, direction = "wide", timevar = "Time", idvar = "Taxonomy")
rownames(wide_agg_mnorm) <- wide_agg_mnorm$Taxonomy
wide_agg_mnorm <- wide_agg_mnorm[,2:dim(wide_agg_mnorm)[2]]

# Networks
cor.data <- expand.grid(rownames(wide_agg_mnorm), rownames(wide_agg_mnorm))
cor.data <- cor.data[which((cor.data$Var1 == cor.data$Var2) == F),]
test.dups <- duplicated(t(apply(cor.data, 1, sort)))
cor.data <- cor.data[which(test.dups == F), ]

correlations <- c()
for(i in 1:dim(cor.data)[1]){
  thing1 <- wide_agg_mnorm[which(rownames(wide_agg_mnorm) == cor.data$Var1[i]), ]
  thing2 <- wide_agg_mnorm[which(rownames(wide_agg_mnorm) == cor.data$Var2[i]), ]
  correlations[i] <- cor(as.numeric(thing1), as.numeric(thing2))
  if(abs(correlations[i]) < 0.95){
    correlations[i] <- 0
  }
}
cor.data$edges <- correlations
cor.data <- cor.data[which(cor.data$edges != 0), ]
kept_taxa <- unique(c(as.character(cor.data$Var1), as.character(cor.data$Var2)))
vert.data <- data.frame(taxa = kept_taxa)

vert.data$gene_expression <- rowSums(wide_agg_mnorm)[match(vert.data$taxa, rownames(wide_agg_mnorm))]
vert.data$Phylum <- sapply(strsplit(as.character(vert.data$taxa),";"), `[`, 1)
vert.data$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", vert.data$Phylum)
vert.data$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", vert.data$Phylum)
vert.data$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", vert.data$Phylum)
vert.data$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", vert.data$Phylum)
vert.data$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", vert.data$Phylum)
vert.data$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", vert.data$Phylum)

phylum_colors <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "grey")
net <- graph_from_data_frame(cor.data, directed = F, vertices = vert.data)
V(net)$color <- phylum_colors[factor(V(net)$Phylum, levels = unique(V(net)$Phylum))]
V(net)$size <- log(vert.data$gene_expression)
l <- layout_with_lgl(net)
plot(net, vertex.color=V(net)$color, vertex.size = V(net)$size, , vertex.label = NA)
legend(x=-1.75, y=0.5, unique(V(net)$Phylum), pch=21, col="#777777", pt.bg=phylum_colors, pt.cex=2, cex=.8, bty="n", ncol=1)
#

ceb <- cluster_edge_betweenness(net) 
# plot(ceb, net, vertex.color="lightgrey", vertex.size = V(net)$size, vertex.labels = NA)
cluster.nums <- membership(ceb)
vert.data$Cluster <- cluster.nums[match(vert.data$taxa, names(cluster.nums))]
vert.data <- vert.data[order(vert.data$Cluster), ]
cluster.reads <- aggregate(gene_expression ~ Cluster, vert.data, sum)
cluster.reads <- cluster.reads[order(cluster.reads$gene_expression), ]

# Plot cluster over time
c <- 18
cluster.abun <- colSums(wide_agg_mnorm[match(vert.data$taxa[which(vert.data$Cluster == c)], rownames(wide_agg_mnorm)), ])
cluster.abun <- data.frame(abundance = as.numeric(cluster.abun), time = names(cluster.abun))
cluster.abun$time <- factor(cluster.abun$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
cluster.abun$abundance <- zscore(cluster.abun$abundance )
ggplot(data = cluster.abun, aes(x = time, y = abundance)) + geom_bar(stat = "identity", fill = "darkgrey", color = "black") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM")) + labs(x = NULL, y = NULL, title = NULL)

#Output files for grepping genes
genes2grep <- data.frame(mendota_key$Gene[grep("Planctomycetes;$", mendota_key$Taxonomy)])
write.table(genes2grep, "C:/Users/Goose and Gander/Desktop/genes2grep/Mendota_Planctomycetes.txt", row.names = F, quote = F)


# Repeat with Sparkling
spark_key$Taxonomy <- gsub("Bacteria;", "", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("Eukaryota;", "", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("NO CLASSIFICATION MH", "Unclassified", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("NO CLASSIFICATION LP", "Unclassified", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("None", "Unclassified", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("unclassified unclassified", "Unclassified", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("unclassified Oligohymenophorea", "Ciliophora", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("unclassified Pelagophyceae", "Ochrophyta", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("unclassified", "Unclassified", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("Unclassified ", "Unclassified", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("UnclassifiedIsochrysidales", "Haptophyta", spark_key$Taxonomy)
spark_key$Taxonomy[grep("Blank", spark_key$Taxonomy)] <- "Unclassified"
spark_key <- spark_key[which(spark_key$Taxonomy != "Unclassified" & spark_key$Taxonomy != "internal standard" & spark_key$Taxonomy != ""),]

keep <- match(rownames(snorm), spark_key$Gene)
snorm <- snorm[which(is.na(keep) == F),]

snorm$Genes <- rownames(snorm)
snorm <- melt(snorm)
snorm$variable <- gsub(".nonrRNA", "", snorm$variable)
snorm$Time <- metadata$Timepoint[match(snorm$variable, metadata$Sample)]
snorm$Taxonomy <- spark_key$Taxonomy[match(snorm$Genes, spark_key$Gene)]

agg_snorm <- aggregate(value ~ Time + Taxonomy, snorm, mean)

agg_snorm$Time <- factor(agg_snorm$Time, levels = c("0", "4", "8", "12", "16", "20", "24", "28", "32", "36", "40", "44"))

# convert to wide format for clustering
wide_agg_snorm <- reshape(agg_snorm, direction = "wide", timevar = "Time", idvar = "Taxonomy")
rownames(wide_agg_snorm) <- wide_agg_snorm$Taxonomy
wide_agg_snorm <- wide_agg_snorm[,2:dim(wide_agg_snorm)[2]]

# Networks
cor.data <- expand.grid(rownames(wide_agg_snorm), rownames(wide_agg_snorm))
cor.data <- cor.data[which((cor.data$Var1 == cor.data$Var2) == F),]
test.dups <- duplicated(t(apply(cor.data, 1, sort)))
cor.data <- cor.data[which(test.dups == F), ]

correlations <- c()
for(i in 1:dim(cor.data)[1]){
  thing1 <- wide_agg_snorm[which(rownames(wide_agg_snorm) == cor.data$Var1[i]), ]
  thing2 <- wide_agg_snorm[which(rownames(wide_agg_snorm) == cor.data$Var2[i]), ]
  correlations[i] <- cor(as.numeric(thing1), as.numeric(thing2))
  if(abs(correlations[i]) < 0.95){
    correlations[i] <- 0
  }
}
cor.data$edges <- correlations
cor.data <- cor.data[which(cor.data$edges != 0), ]
kept_taxa <- unique(c(as.character(cor.data$Var1), as.character(cor.data$Var2)))
vert.data <- data.frame(taxa = kept_taxa)

vert.data$gene_expression <- rowSums(wide_agg_snorm)[match(vert.data$taxa, rownames(wide_agg_snorm))]
vert.data$Phylum <- sapply(strsplit(as.character(vert.data$taxa),";"), `[`, 1)
vert.data$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", vert.data$Phylum)

phylum_colors <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "grey")
net <- graph_from_data_frame(cor.data, directed = F, vertices = vert.data)
V(net)$color <- phylum_colors[factor(V(net)$Phylum, levels = unique(V(net)$Phylum))]
V(net)$size <- log(vert.data$gene_expression)
l <- layout_with_lgl(net)
plot(net, vertex.color=V(net)$color, vertex.size = V(net)$size, , vertex.label = NA)
legend(x=-1.75, y=0.5, unique(V(net)$Phylum), pch=21, col="#777777", pt.bg=phylum_colors, pt.cex=2, cex=.8, bty="n", ncol=1)
#

ceb <- cluster_edge_betweenness(net) 
plot(ceb, net, vertex.color="lightgrey", vertex.size = V(net)$size, vertex.label = membership(ceb))
cluster.nums <- membership(ceb)
vert.data$Cluster <- cluster.nums[match(vert.data$taxa, names(cluster.nums))]
vert.data <- vert.data[order(vert.data$Cluster), ]
cluster.reads <- aggregate(gene_expression ~ Cluster, vert.data, sum)
cluster.reads <- cluster.reads[order(cluster.reads$gene_expression), ]

# Plot cluster over time
c <- 17
cluster.abun <- colSums(wide_agg_snorm[match(vert.data$taxa[which(vert.data$Cluster == c)], rownames(wide_agg_snorm)), ])
cluster.abun <- data.frame(abundance = as.numeric(cluster.abun), time = names(cluster.abun))
cluster.abun$time <- factor(cluster.abun$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
cluster.abun$abundance <- zscore(cluster.abun$abundance )
ggplot(data = cluster.abun, aes(x = time, y = abundance)) + geom_bar(stat = "identity", fill = "darkgrey", color = "black") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM")) + labs(x = NULL, y = NULL, title = NULL)

#Output files for grepping genes
genes2grep <- data.frame(spark_key$Gene[grep("Microbacteriaceae", spark_key$Taxonomy)])
write.table(genes2grep, "C:/Users/Goose and Gander/Desktop/genes2grep/Spark_Microbacteriaceae.txt", row.names = F, quote = F)


# Repeat with Trout
trout_key$Taxonomy <- gsub("Bacteria;", "", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("Eukaryota;", "", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("NO CLASSIFICATION MH", "Unclassified", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("NO CLASSIFICATION LP", "Unclassified", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("None", "Unclassified", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("unclassified unclassified", "Unclassified", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("unclassified Oligohymenophorea", "Ciliophora", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("unclassified Pelagophyceae", "Ochrophyta", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("unclassified", "Unclassified", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("Unclassified ", "Unclassified", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("UnclassifiedIsochrysidales", "Haptophyta", trout_key$Taxonomy)
trout_key$Taxonomy[grep("Blank", trout_key$Taxonomy)] <- "Unclassified"
trout_key <- trout_key[which(trout_key$Taxonomy != "Unclassified" & trout_key$Taxonomy != "internal standard" & trout_key$Taxonomy != ""),]

keep <- match(rownames(tnorm), trout_key$Gene)
tnorm <- tnorm[which(is.na(keep) == F),]

tnorm$Genes <- rownames(tnorm)
tnorm <- melt(tnorm)
tnorm$variable <- gsub(".nonrRNA", "", tnorm$variable)
tnorm$Time <- metadata$Timepoint[match(tnorm$variable, metadata$Sample)]
tnorm$Taxonomy <- trout_key$Taxonomy[match(tnorm$Genes, trout_key$Gene)]

agg_tnorm <- aggregate(value ~ Time + Taxonomy, tnorm, mean)

agg_tnorm$Time <- factor(agg_tnorm$Time, levels = c("0", "4", "8", "12", "16", "20", "24", "28", "32", "36", "40", "44"))

# convert to wide format for clustering
wide_agg_tnorm <- reshape(agg_tnorm, direction = "wide", timevar = "Time", idvar = "Taxonomy")
rownames(wide_agg_tnorm) <- wide_agg_tnorm$Taxonomy
wide_agg_tnorm <- wide_agg_tnorm[,2:dim(wide_agg_tnorm)[2]]

# Networks
cor.data <- expand.grid(rownames(wide_agg_tnorm), rownames(wide_agg_tnorm))
cor.data <- cor.data[which((cor.data$Var1 == cor.data$Var2) == F),]
test.dups <- duplicated(t(apply(cor.data, 1, sort)))
cor.data <- cor.data[which(test.dups == F), ]

correlations <- c()
for(i in 1:dim(cor.data)[1]){
  thing1 <- wide_agg_tnorm[which(rownames(wide_agg_tnorm) == cor.data$Var1[i]), ]
  thing2 <- wide_agg_tnorm[which(rownames(wide_agg_tnorm) == cor.data$Var2[i]), ]
  correlations[i] <- cor(as.numeric(thing1), as.numeric(thing2))
  if(abs(correlations[i]) < 0.95){
    correlations[i] <- 0
  }
}
cor.data$edges <- correlations
cor.data <- cor.data[which(cor.data$edges != 0), ]
kept_taxa <- unique(c(as.character(cor.data$Var1), as.character(cor.data$Var2)))
vert.data <- data.frame(taxa = kept_taxa)

vert.data$gene_expression <- rowSums(wide_agg_tnorm)[match(vert.data$taxa, rownames(wide_agg_tnorm))]
vert.data$Phylum <- sapply(strsplit(as.character(vert.data$taxa),";"), `[`, 1)
vert.data$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", vert.data$Phylum)
vert.data$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", vert.data$Phylum)

phylum_colors <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "grey")
net <- graph_from_data_frame(cor.data, directed = F, vertices = vert.data)
V(net)$color <- phylum_colors[factor(V(net)$Phylum, levels = unique(V(net)$Phylum))]
V(net)$size <- log(vert.data$gene_expression)
l <- layout_with_lgl(net)
plot(net, vertex.color=V(net)$color, vertex.size = V(net)$size, , vertex.label = membership(ceb))
legend(x=-1.75, y=0.5, unique(V(net)$Phylum), pch=21, col="#777777", pt.bg=phylum_colors, pt.cex=2, cex=.8, bty="n", ncol=1)
#

ceb <- cluster_edge_betweenness(net) 
plot(ceb, net, vertex.color="lightgrey", vertex.size = V(net)$size, vertex.label = membership(ceb))
cluster.nums <- membership(ceb)
vert.data$Cluster <- cluster.nums[match(vert.data$taxa, names(cluster.nums))]
vert.data <- vert.data[order(vert.data$Cluster), ]
cluster.reads <- aggregate(gene_expression ~ Cluster, vert.data, sum)
cluster.reads <- cluster.reads[order(cluster.reads$gene_expression), ]

# Plot cluster over time
c <- 14
cluster.abun <- colSums(wide_agg_tnorm[match(vert.data$taxa[which(vert.data$Cluster == c)], rownames(wide_agg_tnorm)), ])
cluster.abun <- data.frame(abundance = as.numeric(cluster.abun), time = names(cluster.abun))
cluster.abun$time <- factor(cluster.abun$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
cluster.abun$abundance <- zscore(cluster.abun$abundance )
ggplot(data = cluster.abun, aes(x = time, y = abundance)) + geom_bar(stat = "identity", fill = "darkgrey", color = "black") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM")) + labs(x = NULL, y = NULL, title = NULL)

#Output files for grepping genes
genes2grep <- data.frame(trout_key$Gene[grep("Flavo-A2", trout_key$Taxonomy)])
write.table(genes2grep, "C:/Users/Goose and Gander/Desktop/genes2grep/Trout_Chromobacteriaceae.txt", row.names = F, quote = F)
