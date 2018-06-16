# GEODES thesis chapter code

library(ggplot2)
library(cowplot)
library(reshape2)
library(GeneCycle)

zscore <- function(counts){
  z <- (counts - mean(counts)) / sd(counts)
  return(z)
}

#snorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
#tnorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mnorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

#snorm <- snorm[, which(colnames(snorm) != "GEODES014.nonrRNA" & colnames(snorm) != "GEODES033.nonrRNA")]
#tnorm <- tnorm[, which(colnames(tnorm) != "GEODES065.nonrRNA")]
mnorm <- mnorm[, which(colnames(mnorm) != "GEODES158.nonrRNA")]

mendota_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#spark_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#trout_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

#### Prove most things aren't diel

### Abundance filtering
# Before I begin a network analysis, I want to select only genes with cyclic trends
# First screen is to remove genes with too few reads to analyze a trend
# I'm arbitrarily defining this as the number of samples * 1000
#Note: change this to the minimum needed to differentiate between replicates

abun_mnorm <- mnorm[which(rowSums(mnorm) > (dim(mnorm)[2] * 50000)), ]
#abun_snorm <- snorm[which(rowSums(snorm) > (dim(snorm)[2] * 30000)), ]
#abun_tnorm <- tnorm[which(rowSums(tnorm) > (dim(tnorm)[2] * 30000)), ]

# rm(snorm)
# rm(tnorm)

# How much does abundance threshold of 10000 remove?
# Sparkling: 2,887,392 to 1,601,700 genes

# 20000?
#S: 3 mil to 900000
#T: 1 mil to 640,000
#M: 3 mil to 1,100,000

# Fourier transformations

abun_mnorm$Genes <- rownames(abun_mnorm)
abun_mnorm <- melt(abun_mnorm)
abun_mnorm$variable <- gsub(".nonrRNA", "", abun_mnorm$variable)
abun_mnorm$Timepoint <- metadata$Timepoint[match(abun_mnorm$variable, metadata$Sample)]
agg_abun_mnorm <- aggregate(value ~ Genes + Timepoint, data = abun_mnorm, mean)
new_abun_mnorm <- reshape(agg_abun_mnorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(new_abun_mnorm) <- new_abun_mnorm[, 1]
new_abun_mnorm <- new_abun_mnorm[, 2:dim(new_abun_mnorm)[2]]


# abun_snorm$Genes <- rownames(abun_snorm)
# abun_snorm <- melt(abun_snorm)
# abun_snorm$variable <- gsub(".nonrRNA", "", abun_snorm$variable)
# abun_snorm$Timepoint <- metadata$Timepoint[match(abun_snorm$variable, metadata$Sample)]
# agg_abun_snorm <- aggregate(value ~ Genes + Timepoint, data = abun_snorm, mean)
# new_abun_snorm <- reshape(agg_abun_snorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
# rownames(new_abun_snorm) <- new_abun_snorm[, 1]
# new_abun_snorm <- new_abun_snorm[, 2:dim(new_abun_snorm)[2]]


# abun_tnorm$Genes <- rownames(abun_tnorm)
# abun_tnorm <- melt(abun_tnorm)
# abun_tnorm$variable <- gsub(".nonrRNA", "", abun_tnorm$variable)
# abun_tnorm$Timepoint <- metadata$Timepoint[match(abun_tnorm$variable, metadata$Sample)]
# agg_abun_tnorm <- aggregate(value ~ Genes + Timepoint, data = abun_tnorm, mean)
# new_abun_tnorm <- reshape(agg_abun_tnorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
# rownames(new_abun_tnorm) <- new_abun_tnorm[, 1]
# new_abun_tnorm <- new_abun_tnorm[, 2:dim(new_abun_tnorm)[2]]
# new_abun_tnorm <- new_abun_tnorm[, 1:7]


# Then run the significance test for cyclic trends
# fdr.spark <- fdrtool(fisher.g.test(t(new_abun_snorm)), statistic = "pvalue")
# length(which(fdr.spark$pval < 0.05))/length(fdr.spark$pval)
# quantile(fdr.spark$qval)

# fdr.trout <- fdrtool(fisher.g.test(t(new_abun_tnorm)), statistic = "pvalue")
# length(which(fdr.trout$pval < 0.05))/length(fdr.trout$pval)
# quantile(fdr.trout$qval)

fdr.mendota <- fdrtool(fisher.g.test(t(new_abun_mnorm)), statistic = "pvalue")
length(which(fdr.mendota$pval < 0.05))/length(fdr.mendota$pval)
quantile(fdr.mendota$qval)

sig.mendota <- t(new_abun_mnorm[which(fdr.mendota$qval < 0.05),])
sig.mendota.key <- mendota_key[match(colnames(sig.mendota), mendota_key$Gene), ]
table(as.character(sig.mendota.key[,4]))


#### Ok, that's well and good. so is anything diel?
# Re-run code below for each gene and lake
#searchterm <- c("rhodopsin|Rhodopsin")
#searchterm <- c("phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene")
#searchterm <- c("photosynth|Photosynth")
#searchterm <- c("Rubisco|RuBisCO|rubisco|bisphosphate carboxylase")
#searchterm <- c("sugar transport|ose transport")
searchterm <- c("peroxidase|peroxide|catalase")


marker_genes <- mendota_key[grep(searchterm, mendota_key$Product), ]

marker_genes_mnorm <- mnorm[match(marker_genes$Gene, rownames(mnorm)), ]
marker_genes_mnorm <- marker_genes_mnorm[which(rowSums(marker_genes_mnorm) > (dim(marker_genes_mnorm)[2] * 1000)), ]

#Aggregate by timepoint
marker_genes_mnorm$Genes <- rownames(marker_genes_mnorm)
marker_genes_mnorm <- melt(marker_genes_mnorm)
marker_genes_mnorm$variable <- gsub(".nonrRNA", "", marker_genes_mnorm$variable)
marker_genes_mnorm$Timepoint <- metadata$Timepoint[match(marker_genes_mnorm$variable, metadata$Sample)]
agg_marker_genes_mnorm <- aggregate(value ~ Genes + Timepoint, data = marker_genes_mnorm, mean)
new_marker_genes_mnorm <- reshape(agg_marker_genes_mnorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(new_marker_genes_mnorm) <- new_marker_genes_mnorm[, 1]
new_marker_genes_mnorm <- new_marker_genes_mnorm[, 2:dim(new_marker_genes_mnorm)[2]]

MErhodopsin <- zscore(colSums(new_marker_genes_mnorm)) #NO
MErhodo_biosynth <- zscore(colSums(new_marker_genes_mnorm)) #YES
MEphoto <- zscore(colSums(new_marker_genes_mnorm)) #SORT OF
MErubisco <- zscore(colSums(new_marker_genes_mnorm)) #NO
MEsugar <- zscore(colSums(new_marker_genes_mnorm)) #NO
MEperoxidase <- zscore(colSums(new_marker_genes_mnorm)) #MOSTLY

#Plot as heatmaps
diel_mendota_genes <- rep(c("rhodopsin_biosynthesis", "photosynthesis", "reactive_oxygen_defense"), each = 12)
diel_mendota_counts <- c(MErhodo_biosynth, MEphoto, MEperoxidase)
diel_mendota_time <- c(names(MErhodo_biosynth), names(MEphoto), names(MEperoxidase))
diel_mendota <- data.frame(Genes = diel_mendota_genes, Counts = diel_mendota_counts, Time = diel_mendota_time)
ggplot(diel_mendota, aes(x = Time, y = Genes, fill = Counts)) + geom_tile() + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow")
