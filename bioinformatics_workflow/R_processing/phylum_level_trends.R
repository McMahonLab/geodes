# Plot phyla aggregated trends over a day/night cycle
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

# Filter by phyla level data

mendota_key$Taxonomy <- gsub("Bacteria;", "", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("Eukaryota;", "", mendota_key$Taxonomy)
mendota_key$Phylum <- sapply(strsplit(as.character(mendota_key$Taxonomy),";"), `[`, 1)
mendota_key$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("None", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified unclassified", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Unclassified ", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", mendota_key$Phylum)
mendota_key$Phylum[grep("Blank", mendota_key$Phylum)] <- "Unclassified"
mendota_key <- mendota_key[which(mendota_key$Phylum != "Unclassified" & mendota_key$Phylum != "internal standard"),]

spark_key$Taxonomy <- gsub("Bacteria;", "", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("Eukaryota;", "", spark_key$Taxonomy)
spark_key$Phylum <- sapply(strsplit(as.character(spark_key$Taxonomy),";"), `[`, 1)
spark_key$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", spark_key$Phylum)
spark_key$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", spark_key$Phylum)
spark_key$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", spark_key$Phylum)
spark_key$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", spark_key$Phylum)
spark_key$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", spark_key$Phylum)
spark_key$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("None", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified unclassified", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("Unclassified ", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", spark_key$Phylum)
spark_key$Phylum[grep("Blank", spark_key$Phylum)] <- "Unclassified"
spark_key <- spark_key[which(spark_key$Phylum != "Unclassified" & spark_key$Phylum != "internal standard"),]

trout_key$Taxonomy <- gsub("Bacteria;", "", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("Eukaryota;", "", trout_key$Taxonomy)
trout_key$Phylum <- sapply(strsplit(as.character(trout_key$Taxonomy),";"), `[`, 1)
trout_key$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", trout_key$Phylum)
trout_key$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", trout_key$Phylum)
trout_key$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", trout_key$Phylum)
trout_key$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", trout_key$Phylum)
trout_key$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", trout_key$Phylum)
trout_key$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("None", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified unclassified", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("Unclassified ", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", trout_key$Phylum)
trout_key$Phylum[grep("Blank", trout_key$Phylum)] <- "Unclassified"
trout_key <- trout_key[which(trout_key$Phylum != "Unclassified" & trout_key$Phylum != "internal standard"),]

keep <- match(rownames(mnorm), mendota_key$Gene)
mnorm <- mnorm[which(is.na(keep) == F),]

keep <- match(rownames(snorm), spark_key$Gene)
snorm <- snorm[which(is.na(keep) == F),]

keep <- match(rownames(tnorm), trout_key$Gene)
tnorm <- tnorm[which(is.na(keep) == F),]

# Melt data tables and add designation for phylum and time of daty
mnorm$Genes <- rownames(mnorm)
mnorm <- melt(mnorm)
mnorm$variable <- gsub(".nonrRNA", "", mnorm$variable)
mnorm$Time <- metadata$Timepoint[match(mnorm$variable, metadata$Sample)]
mnorm$Phylum <- mendota_key$Phylum[match(mnorm$Genes, mendota_key$Gene)]

snorm$Genes <- rownames(snorm)
snorm <- melt(snorm)
snorm$variable <- gsub(".nonrRNA", "", snorm$variable)
snorm$Time <- metadata$Timepoint[match(snorm$variable, metadata$Sample)]
snorm$Phylum <- spark_key$Phylum[match(snorm$Genes, spark_key$Gene)]

tnorm$Genes <- rownames(tnorm)
tnorm <- melt(tnorm)
tnorm$variable <- gsub(".nonrRNA", "", tnorm$variable)
tnorm$Time <- metadata$Timepoint[match(tnorm$variable, metadata$Sample)]
tnorm$Phylum <- trout_key$Phylum[match(tnorm$Genes, trout_key$Gene)]


#Aggregate by time and phylum

agg_mnorm <- aggregate(value ~ Time + Phylum, mnorm, mean)
agg_snorm <- aggregate(value ~ Time + Phylum, snorm, mean)
agg_tnorm <- aggregate(value ~ Time + Phylum, tnorm, mean)

agg_mnorm$Time <- factor(agg_mnorm$Time, levels = c("0", "4", "8", "12", "16", "20", "24", "28", "32", "36", "40", "44"))
agg_snorm$Time <- factor(agg_snorm$Time, levels = c("0", "4", "8", "12", "16", "20", "24", "28", "32", "36", "40", "44"))
agg_tnorm$Time <- factor(agg_tnorm$Time, levels = c("0", "4", "8", "12", "16", "20", "24", "28", "32", "36", "40", "44"))

#Change to fold change
phyla <- unique(agg_mnorm$Phylum)
for(i in 1:length(phyla)){
  subset <- agg_mnorm[which(agg_mnorm$Phylum == phyla[i]), ]
  ave <- mean(subset$value)
  subset$value <- subset$value - ave
  maxrange <- max(abs(subset$value))
  subset$value <- subset$value * (10/maxrange)
  agg_mnorm[which(agg_mnorm$Phylum == phyla[i]), ] <- subset
}

phyla <- unique(agg_tnorm$Phylum)
for(i in 1:length(phyla)){
  subset <- agg_tnorm[which(agg_tnorm$Phylum == phyla[i]), ]
  ave <- mean(subset$value)
  subset$value <- subset$value - ave
  maxrange <- max(abs(subset$value))
  subset$value <- subset$value * (10/maxrange)
  agg_tnorm[which(agg_tnorm$Phylum == phyla[i]), ] <- subset
}

phyla <- unique(agg_snorm$Phylum)
for(i in 1:length(phyla)){
  subset <- agg_snorm[which(agg_snorm$Phylum == phyla[i]), ]
  ave <- mean(subset$value)
  subset$value <- subset$value - ave
  maxrange <- max(abs(subset$value))
  subset$value <- subset$value * (10/maxrange)
  agg_snorm[which(agg_snorm$Phylum == phyla[i]), ] <- subset
}

# Pick the most expressed phyla to plot
# me_expressed <- c("Chloroflexi", "Heterokonta", "Cyanobacteria", "Bacteroidetes", "Gemmatimonadetes", "Cryptophyta", "Crenarchaeaota", "Actinobacteria", "Viruses", "Proteobacteria")
# agg_mnorm <- agg_mnorm[which(agg_mnorm$Phylum %in% me_expressed),]
# agg_mnorm$Phylum <- factor(agg_mnorm$Phylum, levels = rev(me_expressed))
# 
# sp_expressed <- c("Heterokonta", "Cyanobacteria", "TM7", "Armatimonadetes", "Bacteroidetes", "Proteobacteria", "Cryptophyta",  "Haptophyta", "Actinobacteria", "Elusimicrobia")
# agg_snorm <- agg_snorm[which(agg_snorm$Phylum %in% sp_expressed),]
# agg_snorm$Phylum <- factor(agg_snorm$Phylum, levels = rev(sp_expressed))
# 
# tb_expressed <- c("Cyanobacteria", "Armatimonadetes", "Streptophyta", "Cryptophyta", "Verrucomicrobia", "Deinococcus-Thermus", "Bacteroidetes", "Actinobacteria", "Proteobacteria", "Chlorophyta")
# agg_tnorm <- agg_tnorm[which(agg_tnorm$Phylum %in% tb_expressed),]
# agg_tnorm$Phylum <- factor(agg_tnorm$Phylum, levels = rev(tb_expressed))

ggplot(agg_mnorm, aes(x = Time, y = Phylum, fill = value)) + geom_tile() + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow")

ggplot(agg_snorm, aes(x = Time, y = Phylum, fill = value)) + geom_tile() + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow")

ggplot(agg_tnorm, aes(x = Time, y = Phylum, fill = value)) + geom_tile() + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow")



# convert to wide format for clustering
wide_agg_mnorm <- reshape(agg_mnorm, direction = "wide", timevar = "Time", idvar = "Phylum")
rownames(wide_agg_mnorm) <- wide_agg_mnorm$Phylum
wide_agg_mnorm <- wide_agg_mnorm[,2:dim(wide_agg_mnorm)[2]]

wide_agg_tnorm <- reshape(agg_tnorm, direction = "wide", timevar = "Time", idvar = "Phylum")
rownames(wide_agg_tnorm) <- wide_agg_tnorm$Phylum
wide_agg_tnorm <- wide_agg_tnorm[,2:dim(wide_agg_tnorm)[2]]

wide_agg_snorm <- reshape(agg_snorm, direction = "wide", timevar = "Time", idvar = "Phylum")
rownames(wide_agg_snorm) <- wide_agg_snorm$Phylum
wide_agg_snorm <- wide_agg_snorm[,2:dim(wide_agg_snorm)[2]]

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
  if(abs(correlations[i]) < 0.8){
    correlations[i] <- 0
  }
}
cor.data$edges <- correlations
cor.data <- cor.data[which(cor.data$edges != 0), ]
kept_phyla <- unique(c(as.character(cor.data$Var1), as.character(cor.data$Var2)))
vert.data <- data.frame(phyla = kept_phyla)

vert.data$gene_expression <- rowSums(wide_agg_mnorm)[match(vert.data$phyla, rownames(wide_agg_mnorm))]


phylum_colors <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122")
net <- graph_from_data_frame(cor.data, directed = F, vertices = vert.data)
V(net)$color <- phylum_colors[factor(V(net)$phyla, levels = unique(V(net)$phyla))]
V(net)$size <- log(vert.data$gene_expression)
l <- layout_with_lgl(net)
plot(net, vertex.color="lightgrey", vertex.size = V(net)$size)
legend(x=-1.75, y=0.5, kept_phyla, pch=21, col="#777777", pt.bg=phylum_colors, pt.cex=2, cex=.8, bty="n", ncol=1)
#

ceb <- cluster_edge_betweenness(net) 
# plot(ceb, net, vertex.color="lightgrey", vertex.size = V(net)$size)
cluster.nums <- membership(ceb)
vert.data$Cluster <- cluster.nums[match(vert.data$phyla, names(cluster.nums))]
vert.data <- vert.data[order(vert.data$Cluster), ]
cluster.reads <- aggregate(gene_expression ~ Cluster, vert.data, sum)
cluster.reads <- cluster.reads[order(cluster.reads$gene_expression), ]

# Repeat with Trout
cor.data <- expand.grid(rownames(wide_agg_tnorm), rownames(wide_agg_tnorm))
cor.data <- cor.data[which((cor.data$Var1 == cor.data$Var2) == F),]
test.dups <- duplicated(t(apply(cor.data, 1, sort)))
cor.data <- cor.data[which(test.dups == F), ]

correlations <- c()
for(i in 1:dim(cor.data)[1]){
  thing1 <- wide_agg_tnorm[which(rownames(wide_agg_tnorm) == cor.data$Var1[i]), ]
  thing2 <- wide_agg_tnorm[which(rownames(wide_agg_tnorm) == cor.data$Var2[i]), ]
  correlations[i] <- cor(as.numeric(thing1), as.numeric(thing2))
  if(abs(correlations[i]) < 0.8){
    correlations[i] <- 0
  }
}
cor.data$edges <- correlations
cor.data <- cor.data[which(cor.data$edges != 0), ]
kept_phyla <- unique(c(as.character(cor.data$Var1), as.character(cor.data$Var2)))
vert.data <- data.frame(phyla = kept_phyla)

vert.data$gene_expression <- rowSums(wide_agg_tnorm)[match(vert.data$phyla, rownames(wide_agg_tnorm))]


phylum_colors <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122")
net <- graph_from_data_frame(cor.data, directed = F, vertices = vert.data)
V(net)$color <- phylum_colors[factor(V(net)$phyla, levels = unique(V(net)$phyla))]
V(net)$size <- log(vert.data$gene_expression)
l <- layout_with_lgl(net)
plot(net, vertex.color="lightgrey", vertex.size = V(net)$size)
#legend(x=-1.75, y=0.5, kept_phyla, pch=21, col="#777777", pt.bg=phylum_colors, pt.cex=2, cex=.8, bty="n", ncol=1)
#

ceb <- cluster_edge_betweenness(net) 
plot(ceb, net, vertex.color="lightgrey", vertex.size = V(net)$size)
cluster.nums <- membership(ceb)
vert.data$Cluster <- cluster.nums[match(vert.data$phyla, names(cluster.nums))]
vert.data <- vert.data[order(vert.data$Cluster), ]
cluster.reads <- aggregate(gene_expression ~ Cluster, vert.data, sum)
cluster.reads <- cluster.reads[order(cluster.reads$gene_expression), ]

# Repeat with Sparkling
cor.data <- expand.grid(rownames(wide_agg_snorm), rownames(wide_agg_snorm))
cor.data <- cor.data[which((cor.data$Var1 == cor.data$Var2) == F),]
test.dups <- duplicated(t(apply(cor.data, 1, sort)))
cor.data <- cor.data[which(test.dups == F), ]

correlations <- c()
for(i in 1:dim(cor.data)[1]){
  thing1 <- wide_agg_snorm[which(rownames(wide_agg_snorm) == cor.data$Var1[i]), ]
  thing2 <- wide_agg_snorm[which(rownames(wide_agg_snorm) == cor.data$Var2[i]), ]
  correlations[i] <- cor(as.numeric(thing1), as.numeric(thing2))
  if(abs(correlations[i]) < 0.8){
    correlations[i] <- 0
  }
}
cor.data$edges <- correlations
cor.data <- cor.data[which(cor.data$edges != 0), ]
kept_phyla <- unique(c(as.character(cor.data$Var1), as.character(cor.data$Var2)))
vert.data <- data.frame(phyla = kept_phyla)

vert.data$gene_expression <- rowSums(wide_agg_snorm)[match(vert.data$phyla, rownames(wide_agg_snorm))]


phylum_colors <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122")
net <- graph_from_data_frame(cor.data, directed = F, vertices = vert.data)
V(net)$color <- phylum_colors[factor(V(net)$phyla, levels = unique(V(net)$phyla))]
V(net)$size <- log(vert.data$gene_expression)
l <- layout_with_lgl(net)
plot(net, vertex.color="lightgrey", vertex.size = V(net)$size)
#legend(x=-1.75, y=0.5, kept_phyla, pch=21, col="#777777", pt.bg=phylum_colors, pt.cex=2, cex=.8, bty="n", ncol=1)
#

ceb <- cluster_edge_betweenness(net) 
plot(ceb, net, vertex.color="lightgrey", vertex.size = V(net)$size)
cluster.nums <- membership(ceb)
vert.data$Cluster <- cluster.nums[match(vert.data$phyla, names(cluster.nums))]
vert.data <- vert.data[order(vert.data$Cluster), ]
cluster.reads <- aggregate(gene_expression ~ Cluster, vert.data, sum)
cluster.reads <- cluster.reads[order(cluster.reads$gene_expression), ]


c = 4
cluster.abun <- colSums(wide_agg_snorm[match(vert.data$phyla[which(vert.data$Cluster == c)], rownames(wide_agg_snorm)), ])
cluster.abun <- data.frame(abundance = as.numeric(cluster.abun), time = names(cluster.abun))
cluster.abun$time <- factor(cluster.abun$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
cluster.abun$abundance <- zscore(cluster.abun$abundance )
ggplot(data = cluster.abun, aes(x = time, y = abundance)) + geom_bar(stat = "identity", fill = "#5289C7") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM")) + labs(x = NULL, y = NULL, title = NULL)


# What genes are correlated in these clusters?

grouping <- mnorm[which(mnorm$Phylum == "Chlorobi" | mnorm$Phylum == "Armatimonadetes"), ]
agg_mnorm <- aggregate(value ~ Time + Genes, grouping, mean)
agg_mnorm$Time <- factor(agg_mnorm$Time, levels = c("0", "4", "8", "12", "16", "20", "24", "28", "32", "36", "40", "44"))

agg_mnorm <- reshape(agg_mnorm, idvar = "Genes", timevar = "Time", direction = "wide")
rownames(agg_mnorm) <- agg_mnorm[, 1]
agg_mnorm <- agg_mnorm[, 2:dim(agg_mnorm)[2]]
agg_mnorm <- agg_mnorm[which(rowSums(agg_mnorm) > 5000), ]

# cor.data <- expand.grid(rownames(agg_mnorm), rownames(agg_mnorm))
# cor.data <- cor.data[which((cor.data$Var1 == cor.data$Var2) == F),]
# test.dups <- duplicated(t(apply(cor.data, 1, sort)))
# cor.data <- cor.data[which(test.dups == F), ]
# 
# correlations <- c()
# for(i in 1:dim(cor.data)[1]){
#   thing1 <- agg_mnorm[which(rownames(agg_mnorm) == cor.data$Var1[i]), ]
#   thing2 <- agg_mnorm[which(rownames(agg_mnorm) == cor.data$Var2[i]), ]
#   correlations[i] <- cor(as.numeric(thing1), as.numeric(thing2))
#   if(abs(correlations[i]) < 0.9){
#     correlations[i] <- 0
#   }
# }
# cor.data$edges <- correlations
# cor.data <- cor.data[which(cor.data$edges != 0), ]
# kept_genes <- unique(c(as.character(cor.data$Var1), as.character(cor.data$Var2)))
# vert.data <- mendota_key[match(kept_genes, mendota_key$Gene),]
vert.data <- mendota_key[match(rownames(agg_mnorm), mendota_key$Gene),]
vert.data$gene_expression <- rowSums(agg_mnorm)[match(vert.data$Gene, rownames(agg_mnorm))]
vert.data <- vert.data[order(vert.data$gene_expression, decreasing = T), ]


phylum_colors <- c("green", "blue")
net <- graph_from_data_frame(cor.data, directed = F, vertices = vert.data)
V(net)$color <- phylum_colors[factor(V(net)$Phylum, levels = unique(V(net)$Phylum))]
V(net)$size <- log(vert.data$gene_expression)/2
l <- layout_with_lgl(net)
plot(net, vertex.label = NA, vertex.size = V(net)$size, layout = l, vertex.color=V(net)$color)
legend(x=-1.5, y=0, unique(V(net)$Phylum), pch=21, col="#777777", pt.bg=phylum_colors, pt.cex=2, cex=.8, bty="n", ncol=1)

ceb <- cluster_edge_betweenness(net) 
cluster.nums <- membership(ceb)
vert.data$Cluster <- cluster.nums[match(vert.data$Gene, names(cluster.nums))]
vert.data <- vert.data[order(vert.data$Cluster), ]

cluster.reads <- aggregate(gene_expression ~ Cluster, vert.data, sum)
cluster.reads <- cluster.reads[order(cluster.reads$gene_expression), ]
cluster.reads

# output list of expressed genes in Chlorobi and Armata

armati_genes <- data.frame(vert.data$Gene[which(vert.data$Phylum == "Armatimonadetes")])
chlorobi_genes <- data.frame(vert.data$Gene[which(vert.data$Phylum == "Chlorobi")])
write.table(chlorobi_genes, "C:/Users/Goose and Gander/Desktop/chlorobi_genes.txt", row.names = F, quote = F)
