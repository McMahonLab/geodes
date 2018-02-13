# Fourier + WGCNA

### Load packages
library(ggplot2)
library(cowplot)
library(reshape2)
library(GeneCycle)
library(WGCNA)
allowWGCNAThreads()
enableWGCNAThreads()

### Load data (start with only one to save RAM and comment the rest out)
# Normalized read tables
#snorm <- read.csv("D:/geodes_data_tables/Sparkling_normalized.csv", header = T, row.names = 1)
#tnorm <- read.csv("D:/geodes_data_tables/Trout_normalized.csv", header = T, row.names = 1)
mnorm <- read.csv("D:/geodes_data_tables/Mendota_normalized.csv", header = T, row.names = 1)

# Gene keys
mendota_key <- read.csv("D:/geodes_data_tables/Mendota_ID90_genekey.csv", header = T)
#spark_key <- read.csv("D:/geodes_data_tables/Sparkling_ID90_genekey.csv", header = T)
#trout_key <- read.csv("D:/geodes_data_tables/Trout_ID90_genekey.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Alex/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

### Abundance filtering
# Before I begin a network analysis, I want to select only genes with cyclic trends
# First screen is to remove genes with too few reads to analyze a trend
# I'm arbitrarily defining this as the number of samples * 1000

abun_mnorm <- mnorm[which(rowSums(mnorm) > (dim(mnorm)[2] * 1000)), ]
#abun_snorm <- snorm[which(rowSums(snorm) > (dim(snorm)[2] * 1000)), ]
#abun_tnorm <- tnorm[which(rowSums(tnorm) > (dim(tnorm)[2] * 1000)), ]

### Fourier transformations
# I only want to include genes in my network with significant cyclic trends
# Keep in mind that there's too many genes for multiple testing correction, so I'm well-aware that there will be many false positives. Hence, I'm using this as a screen for network analysis, not as results.

# 1st, Average read counts by timepoint
# My approach is to convert to long format, use aggregate(), and convert back to wide format
abun_mnorm$Genes <- rownames(abun_mnorm)
abun_mnorm <- melt(abun_mnorm)
abun_mnorm$Timepoint <- metadata$Timepoint[match(abun_mnorm$variable, metadata$Sample)]
agg_abun_mnorm <- aggregate(value ~ Genes + Timepoint, data = abun_mnorm, mean)
new_abun_mnorm <- reshape(agg_abun_mnorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(new_abun_mnorm) <- new_abun_mnorm[, 1]
new_abun_mnorm <- new_abun_mnorm[, 2:dim(new_abun_mnorm)[2]]

# abun_snorm$Genes <- rownames(abun_snorm)
# abun_snorm <- melt(abun_snorm)
# abun_snorm$Timepoint <- metadata$Timepoint[match(abun_snorm$variable, metadata$Sample)]
# agg_abun_snorm <- aggregate(value ~ Genes + Timepoint, data = abun_snorm, mean)
# new_abun_snorm <- reshape(agg_abun_snorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
# rownames(new_abun_snorm) <- new_abun_snorm[, 1]
# new_abun_snorm <- new_abun_snorm[, 2:dim(new_abun_snorm)[2]]

# abun_tnorm$Genes <- rownames(abun_tnorm)
# abun_tnorm <- melt(abun_tnorm)
# abun_tnorm$Timepoint <- metadata$Timepoint[match(abun_tnorm$variable, metadata$Sample)]
# agg_abun_tnorm <- aggregate(value ~ Genes + Timepoint, data = abun_tnorm, mean)
# new_abun_tnorm <- reshape(agg_abun_tnorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
# rownames(new_abun_tnorm) <- new_abun_tnorm[, 1]
# new_abun_tnorm <- new_abun_tnorm[, 2:dim(new_abun_tnorm)[2]]
# new_abun_tnorm <- new_abun_tnorm[, 1:7]

# Then run the significance test for cyclic trends
fdr.mendota <- fdrtool(fisher.g.test(t(new_abun_mnorm)), statistic = "pvalue")
sig.mendota <- t(new_abun_mnorm[which(fdr.mendota$pval < 0.05),])

# fdr.spark <- fdrtool(fisher.g.test(t(new_abun_snorm)), statistic = "pvalue")
# sig.spark <- t(new_abun_snorm[which(fdr.spark$pval < 0.05),])
# 
# fdr.trout <- fdrtool(fisher.g.test(t(new_abun_tnorm)), statistic = "pvalue")
# sig.trout <- t(new_abun_tnorm[which(fdr.trout$pval < 0.05),])

### WGCNA
# 1st, check the soft threshold in case removing samples changed it.

# mpowers <- c()
# for(i in 1:100){
#   lottery <- sample(1:dim(sig.mendota)[2], 2000, replace = F)
#   sft = pickSoftThreshold(sig.mendota[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
#   mpowers[i] <- sft$powerEstimate
# }
# hist(mpowers)
# mean(mpowers)
# median(mpowers)
# # Mendota Soft threshold = 8 in this case


# spowers <- c()
# for(i in 1:100){
#   lottery <- sample(1:dim(sig.spark)[2], 2000, replace = F)
#   sft = pickSoftThreshold(sig.spark[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
#   spowers[i] <- sft$powerEstimate
# }
# hist(spowers)
# mean(spowers)
# median(spowers)
# 
# # Sparkling soft threshold = 7

# tpowers <- c()
# for(i in 1:100){
#   lottery <- sample(1:dim(sig.trout)[2], 2000, replace = F)
#   sft = pickSoftThreshold(sig.trout[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
#   tpowers[i] <- sft$powerEstimate
# }
# hist(tpowers)
# mean(tpowers)
# median(tpowers)

# Trout soft threshold = 18

# Make the modules
mendota_net <- blockwiseModules(sig.mendota, maxBlockSize = 1000, power = 8, loadTOM = F, saveTOMs = F, networkType = "signed", minModuleSize = 30, numericLabels = T, nThreads = 8, verbose = 3)

# spark_net <- blockwiseModules(sig.spark, maxBlockSize = 1000, power = 7, loadTOM = F, saveTOMs = F, networkType = "signed", minModuleSize = 30, numericLabels = T, nThreads = 8, verbose = 3)

# trout_net <- blockwiseModules(sig.trout, maxBlockSize = 1000, power = 18, loadTOM = F, saveTOMs = F, networkType = "signed", minModuleSize = 30, numericLabels = T, nThreads = 8, verbose = 3)


### Look at what genes ended up in what clusters
# I need to extract relevant info for my genes with significant trends from the genekey file

sig.mendota.key <- mendota_key[match(colnames(sig.mendota), mendota_key$Gene), ]
sig.mendota.key$Product <- as.character(sig.mendota.key$Product)
sig.mendota.key$Cluster <- mendota_net$colors
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

# sig.mendota.key$Phylum <- sapply(strsplit(as.character(sig.mendota.key$Taxonomy),";"), `[`, 1)
# sig.mendota.key$Class<- sapply(strsplit(as.character(sig.mendota.key$Taxonomy),";"), `[`, 2)
# sig.mendota.key$Order <- sapply(strsplit(as.character(sig.mendota.key$Taxonomy),";"), `[`, 3)
# sig.mendota.key$Lineage <- sapply(strsplit(as.character(sig.mendota.key$Taxonomy),";"), `[`, 4)
# sig.mendota.key$Clade <- sapply(strsplit(as.character(sig.mendota.key$Taxonomy),";"), `[`, 5)
# sig.mendota.key$Tribe <- sapply(strsplit(as.character(sig.mendota.key$Taxonomy),";"), `[`, 6)

total_reads <- rowSums(mnorm)
sig.mendota.key$Totals <- total_reads[match(sig.mendota.key$Gene, names(total_reads))]
write.csv(sig.mendota.key, "D:/geodes_data_tables/WGCNA_mendota_results.csv", row.names = F)
write.csv(mendota_net$MEs, "D:/geodes_data_tables/WGCNA_mendota_eigenvectors.csv", row.names = T)

# Calculate correlation matrix
eigenvectors <- mendota_net$MEs
adjacency <- matrix(0, nrow = dim(eigenvectors)[2], ncol = dim(eigenvectors)[2])
rownames(adjacency) <- colnames(adjacency) <- colnames(eigenvectors)
for(i in 1:dim(eigenvectors)[2]){
  for(j in 1:dim(eigenvectors)[2]){
    value <- cor(eigenvectors[1:11,i], eigenvectors[2:12,j])
    if(value > 0.75 | value < -0.75){
      adjacency[i,j] <- value
    }
  }
}



# Build eigennetwork
net <- as.network(x = adjacency, directed = TRUE, loops = FALSE, matrix.type = "adjacency")
plot.network(net, displaylabels = T, vertex.cex = 3, label.pos = 5)

# Make pretty plot of eigengenes
eigenvectors$Timepoint <- rownames(eigenvectors)
long_eigenvectors <- melt(eigenvectors)
plot.colors <- NA
plot.colors[which(long_eigenvectors$value > 0)] <- "green"
plot.colors[which(long_eigenvectors$value < 0)] <- "red"
long_eigenvectors$Sign <- plot.colors
long_eigenvectors$Timepoint <- factor(long_eigenvectors$Timepoint, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
cluster = "ME26"
ggplot(data = long_eigenvectors[which(long_eigenvectors$variable == cluster), ], aes(x = Timepoint, y = value, fill = Sign)) + geom_bar(stat = "identity") + labs(title = cluster) + scale_fill_manual(values = c("green", "red")) + theme(legend.position = "none")

# Get genes and taxonomy from each cluster
x <- sig.mendota.key[which(sig.mendota.key$Cluster == 26),]
x <- x[order(x$Totals),]
x[(dim(x)[1] - 50): dim(x)[1],]
phyla_breakdown <- aggregate(Totals ~ Phylum, x, sum)
ggplot(phyla_breakdown[grep("Unclassified", phyla_breakdown$Phylum, invert = T), ], aes(x = Phylum, y = Totals)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
x


# sig.spark.key <- spark_key[match(colnames(sig.spark), spark_key$Gene), ]
# sig.spark.key$Product <- as.character(sig.spark.key$Product)
# sig.spark.key$Cluster <- spark_net$colors
# sig.spark.key$Taxonomy <- gsub("Bacteria;", "", sig.spark.key$Taxonomy)
# sig.spark.key$Phylum <- sapply(strsplit(as.character(sig.spark.key$Taxonomy),";"), `[`, 1)
# sig.spark.key$Class<- sapply(strsplit(as.character(sig.spark.key$Taxonomy),";"), `[`, 2)
# sig.spark.key$Order <- sapply(strsplit(as.character(sig.spark.key$Taxonomy),";"), `[`, 3)
# sig.spark.key$Lineage <- sapply(strsplit(as.character(sig.spark.key$Taxonomy),";"), `[`, 4)
# sig.spark.key$Clade <- sapply(strsplit(as.character(sig.spark.key$Taxonomy),";"), `[`, 5)
# sig.spark.key$Tribe <- sapply(strsplit(as.character(sig.spark.key$Taxonomy),";"), `[`, 6)
# 
# total_reads <- rowSums(snorm)
# sig.spark.key$Totals <- total_reads[match(sig.spark.key$Gene, names(total_reads))]

sig.trout.key <- trout_key[match(colnames(sig.trout), trout_key$Gene), ]
sig.trout.key$Product <- as.character(sig.trout.key$Product)
sig.trout.key$Cluster <- trout_net$colors
sig.trout.key$Taxonomy <- gsub("Bacteria;", "", sig.trout.key$Taxonomy)
sig.trout.key$Phylum <- sapply(strsplit(as.character(sig.trout.key$Taxonomy),";"), `[`, 1)
sig.trout.key$Class<- sapply(strsplit(as.character(sig.trout.key$Taxonomy),";"), `[`, 2)
sig.trout.key$Order <- sapply(strsplit(as.character(sig.trout.key$Taxonomy),";"), `[`, 3)
sig.trout.key$Lineage <- sapply(strsplit(as.character(sig.trout.key$Taxonomy),";"), `[`, 4)
sig.trout.key$Clade <- sapply(strsplit(as.character(sig.trout.key$Taxonomy),";"), `[`, 5)
sig.trout.key$Tribe <- sapply(strsplit(as.character(sig.trout.key$Taxonomy),";"), `[`, 6)

total_reads <- rowSums(tnorm)
sig.trout.key$Totals <- total_reads[match(sig.trout.key$Gene, names(total_reads))]

# For the top 10 clusters, what taxonomic groups are there? What products? What overall trend?
# cluster = 2
# barplot(mendota_net$MEs$ME26)
# table(sig.mendota.key$Phylum[which(sig.mendota.key$Cluster == cluster)])
# x <- as.data.frame(table(sig.mendota.key$Product[which(sig.mendota.key$Cluster == cluster)]))
# x <- x[order(x$Freq),]
# x[(dim(x)[1] - 11): dim(x)[1], ]

# How many reads are assigned to each cluster?
cluster_sums <- aggregate(Totals ~ Cluster, sig.trout.key, sum)
cluster_sums[order(cluster_sums$Totals),]

# Group by trends!
cluster = 2
barplot(mendota_net$MEs$ME26)
x <- sig.mendota.key[which(sig.mendota.key$Cluster == 26),]
x <- x[order(x$Totals),]
x[(dim(x)[1] - 50): dim(x)[1],]
phyla_breakdown <- aggregate(Totals ~ Phylum, x, sum)
ggplot(phyla_breakdown[grep("NO CLASSIFICATION", phyla_breakdown$Phylum, invert = T), ], aes(x = Phylum, y = Totals)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))

#That worked great! Repeat in Sparkling

# cluster = 2
# barplot(spark_net$MEs$ME40)
# x <- sig.spark.key[which(sig.spark.key$Cluster == 40),]
# x <- x[order(x$Totals),]
# x[(dim(x)[1] - 50): dim(x)[1],]
# phyla_breakdown <- aggregate(Totals ~ Phylum, x, sum)
# ggplot(phyla_breakdown[grep("NO CLASSIFICATION", phyla_breakdown$Phylum, invert = T), ], aes(x = Phylum, y = Totals)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))


# Trout Bog
barplot(trout_net$MEs$ME0)
x <- sig.trout.key[which(sig.trout.key$Cluster == 0),]
x <- x[order(x$Totals),]
x[(dim(x)[1] - 50): dim(x)[1],]
phyla_breakdown <- aggregate(Totals ~ Phylum, x, sum)
ggplot(phyla_breakdown[grep("NO CLASSIFICATION", phyla_breakdown$Phylum, invert = T), ], aes(x = Phylum, y = Totals)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust =1))

