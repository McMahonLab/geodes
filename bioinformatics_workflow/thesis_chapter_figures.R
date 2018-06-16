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
searchterm <- c("phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene")
#searchterm <- c("photosynth|Photosynth")
#searchterm <- c("Rubisco|RuBisCO|rubisco|bisphosphate carboxylase")
#searchterm <- c("sugar transport|ose transport")
#searchterm <- c("peroxidase|peroxide|catalase")


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

###
# barplot(zscore(colSums(new_marker_genes_mnorm)))
# MErhodopsin <- zscore(colSums(new_marker_genes_mnorm)) #NO
# MErhodo_biosynth <- zscore(colSums(new_marker_genes_mnorm)) #YES
# MEphoto <- zscore(colSums(new_marker_genes_mnorm)) #SORT OF
# MErubisco <- zscore(colSums(new_marker_genes_mnorm)) #NO
# MEsugar <- zscore(colSums(new_marker_genes_mnorm)) #NO
# MEperoxidase <- zscore(colSums(new_marker_genes_mnorm)) #MOSTLY
# 
df <- data.frame(time = names(MErhodo_biosynth), counts = MErhodo_biosynth)
df$time <- factor(df$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
ME1 <- ggplot(df, aes(x = time, y = counts)) + geom_bar(stat = "identity", fill = "#5ab4ac", color = "black") + labs(x = NULL, y = "Z-score Normalized Transcripts/L") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM"))
# 
# df <- data.frame(time = names(MEphoto), counts = MEphoto)
# df$time <- factor(df$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
# ME2 <- ggplot(df, aes(x = time, y = counts)) + geom_bar(stat = "identity", fill = "#5ab4ac", color = "black") + labs(x = NULL, y = "Z-score Normalized Transcripts/L") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM"))
# 
# df <- data.frame(time = names(MEperoxidase), counts = MEperoxidase)
# df$time <- factor(df$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
# ME3 <- ggplot(df, aes(x = time, y = counts)) + geom_bar(stat = "identity", fill = "#5ab4ac", color = "black") + labs(x = NULL, y = "Z-score Normalized Transcripts/L") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM"))

###
# barplot(zscore(colSums(new_marker_genes_snorm)))
# SPrhodopsin <- zscore(colSums(new_marker_genes_snorm)) #NO
# SPrhodo_biosynth <- zscore(colSums(new_marker_genes_snorm)) #NO
# SPphoto <- zscore(colSums(new_marker_genes_snorm)) #YES BUT FLIPPED
# SPrubisco <- zscore(colSums(new_marker_genes_snorm)) #NO
# SPsugar <- zscore(colSums(new_marker_genes_snorm)) #NO
# SPperoxidase <- zscore(colSums(new_marker_genes_snorm)) #NO
# 
# df <- data.frame(time = names(SPrhodo_biosynth), counts = SPrhodo_biosynth)
# df$time <- factor(df$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
# SP1 <- ggplot(df, aes(x = time, y = counts)) + geom_bar(stat = "identity", fill = "#67a9cf", color = "black") + labs(x = NULL, y = "Z-score Normalized Transcripts/L") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM"))
# 
# df <- data.frame(time = names(SPphoto), counts = SPphoto)
# df$time <- factor(df$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
# SP2 <- ggplot(df, aes(x = time, y = counts)) + geom_bar(stat = "identity", fill = "#67a9cf", color = "black") + labs(x = NULL, y = "Z-score Normalized Transcripts/L") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM"))
# 
# df <- data.frame(time = names(SPperoxidase), counts = SPperoxidase)
# df$time <- factor(df$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
# SP3 <- ggplot(df, aes(x = time, y = counts)) + geom_bar(stat = "identity", fill = "#67a9cf", color = "black") + labs(x = NULL, y = "Z-score Normalized Transcripts/L") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM"))

###
barplot(zscore(colSums(new_marker_genes_tnorm)))
TBrhodopsin <- zscore(colSums(new_marker_genes_tnorm)) #NO
TBrhodo_biosynth <- zscore(colSums(new_marker_genes_tnorm)) #NO
TBphoto <- zscore(colSums(new_marker_genes_tnorm)) #YES BUT FLIPPED
TBrubisco <- zscore(colSums(new_marker_genes_tnorm)) #NO
TBsugar <- zscore(colSums(new_marker_genes_tnorm)) #NO
TBperoxidase <- zscore(colSums(new_marker_genes_tnorm)) #NO

df <- data.frame(time = names(TBrhodo_biosynth), counts = TBrhodo_biosynth)
df$time <- factor(df$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
TB1 <- ggplot(df[which(df$time != "value.44"),], aes(x = time, y = counts)) + geom_bar(stat = "identity", fill = "#8c510a", color = "black") + labs(x = NULL, y = "Z-score Normalized Transcripts/L") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM"))

df <- data.frame(time = names(TBphoto), counts = TBphoto)
df$time <- factor(df$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
TB2 <- ggplot(df[which(df$time != "value.44"),], aes(x = time, y = counts)) + geom_bar(stat = "identity", fill = "#8c510a", color = "black") + labs(x = NULL, y = "Z-score Normalized Transcripts/L") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM"))

df <- data.frame(time = names(TBperoxidase), counts = TBperoxidase)
df$time <- factor(df$time, levels = c("value.0", "value.4", "value.8", "value.12", "value.16", "value.20", "value.24", "value.28", "value.32", "value.36", "value.40", "value.44"))
TB3 <- ggplot(df[which(df$time != "value.44"),], aes(x = time, y = counts)) + geom_bar(stat = "identity", fill = "#8c510a", color = "black") + labs(x = NULL, y = "Z-score Normalized Transcripts/L") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM"))

# Put all the plots together

diel_trends <- plot_grid(ME1, SP1, TB1, ME2, SP2, TB2, ME3, SP3, TB3, ncol = 3)


#######
# What are the most expressed genes and taxa in each lake?

totals <- rowSums(snorm)
totals <- totals[order(totals, decreasing = T)]
topgenes <- names(totals[1:1000])
geneinfo <- spark_key[match(topgenes, spark_key$Gene),]
geneinfo <- geneinfo[grep("photo|ribulose-bisphosphate|chloro|hypothetical", geneinfo$Product, invert = T),]
geneinfo <- geneinfo[grep("Cyanobacteria", geneinfo$Taxonomy, invert = T),]
