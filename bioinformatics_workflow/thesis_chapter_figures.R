# GEODES thesis chapter code

library(ggplot2)
library(cowplot)
library(reshape2)
library(GeneCycle)
library(DESeq2)

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


# Now taxa. This one's a lot harder.

theme_set(theme_cowplot(font_size=18))
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# Metagenomic data
metaG_reads <- read.table("E:/geodes_data_tables/GEODES_metaG_ID90_2018-03-10.readcounts.txt", row.names = 1, sep = "\t")
colnames(metaG_reads) <- c("GEODES005", "GEODES006", "GEODES057", "GEODES058", "GEODES117", "GEODES118", "GEODES165", "GEODES166", "GEODES167", "GEODES168")
metaG_key <- read.table("E:/geodes_data_tables/GEODES_metaG_genekey_2018-03-12.txt", sep = "\t", quote = "")
colnames(metaG_key) <- c("Gene", "Genome", "Taxonomy", "Product")
lakekey <- c("Sparkling", "Sparkling", "Trout", "Trout", "Mendota", "Mendota", "Sparkling2009", "Sparkling2009", "Sparkling2009", "Sparkling2009")
metaG_reads <- sweep(metaG_reads, 2, colSums(metaG_reads), "/")

metaG_key$Taxonomy <- gsub("Bacteria;", "", metaG_key$Taxonomy)
metaG_key$Taxonomy <- gsub("Eukaryota;", "", metaG_key$Taxonomy)
metaG_key$Phylum <- sapply(strsplit(as.character(metaG_key$Taxonomy),";"), `[`, 1)

metaG_reads$Genes <- rownames(metaG_reads)
spark_metaG <- metaG_reads[,c(1,2, 11)]
trout_metaG <- metaG_reads[,c(3,4, 11)]
mendota_metaG <- metaG_reads[,c(5,6, 11)]


spark_metaG <- melt(spark_metaG)
trout_metaG <- melt(trout_metaG)
mendota_metaG <- melt(mendota_metaG)

mendota_metaG$Phylum <- metaG_key$Phylum[match(mendota_metaG$Genes, metaG_key$Gene)]
mendota_metaG$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified unclassified unclassified", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("None", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified unclassified", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Unclassified ", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("TM7", "Candidatus Saccharibacteria", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Ignavibacteriae", "Ignavibacteria", mendota_metaG$Phylum)
mendota_metaG$Phylum[which(is.na(mendota_metaG$Phylum) == T)] <- "Unclassified"

mendota_metaG_phyla <- aggregate(value ~ Phylum, data = mendota_metaG, mean)

spark_metaG$Phylum <- metaG_key$Phylum[match(spark_metaG$Genes, metaG_key$Gene)]
spark_metaG$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified unclassified unclassified", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("None", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified unclassified", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Unclassified ", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("TM7", "Candidatus Saccharibacteria", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Ignavibacteriae", "Ignavibacteria", spark_metaG$Phylum)
spark_metaG$Phylum[which(is.na(spark_metaG$Phylum))] <- "Unclassified"
spark_metaG_phyla <- aggregate(value ~ Phylum, data = spark_metaG, mean)

trout_metaG$Phylum <- metaG_key$Phylum[match(trout_metaG$Genes, metaG_key$Gene)]
trout_metaG$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified unclassified unclassified", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("None", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified unclassified", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Unclassified ", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("TM7", "Candidatus Saccharibacteria", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Ignavibacteriae", "Ignavibacteria", trout_metaG$Phylum)
trout_metaG$Phylum[which(is.na(trout_metaG$Phylum))] <- "Unclassified"
trout_metaG_phyla <- aggregate(value ~ Phylum, data = trout_metaG, mean)


# Remove unused metaG datasets to sve RAM
rm(metaG_reads)
rm(mendota_metaG)
rm(trout_metaG)
rm(spark_metaG)
rm(metaG_key)

# Read in mendota data
mnorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mnorm <- mnorm[, which(colnames(mnorm) != "GEODES158.nonrRNA")]

mendota_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)

mendota_key$Taxonomy <- gsub("Bacteria;", "", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("Eukaryota;", "", mendota_key$Taxonomy)
mendota_key$Phylum <- sapply(strsplit(as.character(mendota_key$Taxonomy),";"), `[`, 1)

mnorm$Genes <- rownames(mnorm)
mnorm <- melt(mnorm)
mnorm$variable <- gsub(".nonrRNA", "", mnorm$variable)
mnorm$Timepoint <- metadata$Timepoint[match(mnorm$variable, metadata$Sample)]
mnorm$Taxonomy <- mendota_key$Phylum[match(mnorm$Genes, mendota_key$Gene)]
mnorm$Taxonomy <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("NO CLASSIFICATION MH", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("NO CLASSIFICATION LP", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("None", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified unclassified", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified Oligohymenophorea", "Ciliophora", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified Pelagophyceae", "Ochrophyta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Unclassified ", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("UnclassifiedIsochrysidales", "Haptophyta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", mnorm$Taxonomy)
#Remove unclassifieds to save space
mnorm <- mnorm[which(mnorm$Taxonomy != "Unclassified" & mnorm$Taxonomy != "internal standard"),]
averaged_tax <- aggregate(value ~ Taxonomy + Timepoint, data = mnorm, mean)
wide_mnorm <- reshape(averaged_tax, idvar = "Taxonomy", timevar = "Timepoint", direction = "wide")
rownames(wide_mnorm) <- wide_mnorm$Taxonomy
wide_mnorm <- wide_mnorm[, 2:dim(wide_mnorm)[2]]
wide_mnorm <- wide_mnorm[which(rowSums(wide_mnorm) > 3000),]

mendota_phyla <- data.frame(Taxonomy = rownames(wide_mnorm), Sums = rowSums(wide_mnorm))

rm(mendota_key)
rm(mnorm)

# Choose the intersection of the top 10 expressed and abundant
keep <- c("Chloroflexi", "Actinobacteria", "Ignavibacteria", "Planctomycetes", "Crenarchaeaota", "Cyanobacteria", "Gemmatimonadetes", "Bacteroidetes", "Proteobacteria", "Fibrobacteres", "Cryptophyta", "Heterokonta", "Viruses")
mendota_intersect <- mendota_phyla[which(mendota_phyla$Taxonomy %in% keep), ]
mendota_intersect$metaG <- mendota_metaG_phyla$value[match(mendota_intersect$Taxonomy, mendota_metaG_phyla$Phylum)]
mendota_intersect$Taxonomy <- factor(mendota_intersect$Taxonomy, levels = mendota_intersect$Taxonomy[order(mendota_intersect$metaG, decreasing = F)])
mendota_intersect <- mendota_intersect[order(mendota_intersect$metaG, decreasing = F),]

mendota_intersect$Type <- c("Virus", "Algae", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Archaea", "Bacteria", "Bacteria", "Bacteria", "Bacteria")

p <- ggplot(mendota_intersect[which(mendota_intersect$Taxonomy != "Chloroflexi"), ], aes(x = metaG, y = Sums, color = Type)) + geom_point(size = 2.5) + geom_label_repel(aes(label = Taxonomy, color = Type), force = 5, size = 7.5) + scale_color_manual(values = c("limegreen", "lavenderblush4", "royalblue", "goldenrod")) + labs(x = "Proportion of metagenomic reads assigned", y = "Transcripts/L assigned", title = "Lake Mendota") + theme(legend.position = "none")

save_plot("C:/Users/Goose and Gander/Desktop/geodes/Plots/mendota_dna_vs_rna_no_chloroflexi.pdf", p, base_height = 6, base_aspect_ratio = 8/6)

p <- ggplot(mendota_intersect, aes(x = metaG, y = Sums, color = Type)) + geom_point(size = 2.5) + geom_label_repel(aes(label = Taxonomy, color = Type), force = 5, size = 7.5) + scale_color_manual(values = c("limegreen", "lavenderblush4", "royalblue", "goldenrod")) + labs(x = "Proportion of metagenomic reads assigned", y = "Transcripts/L assigned", title = "Lake Mendota") + theme(legend.position = "none")

save_plot("C:/Users/Goose and Gander/Desktop/geodes/Plots/mendota_dna_vs_rna.pdf", p, base_height = 6, base_aspect_ratio = 8/6)

# Sparkling

snorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
snorm <- snorm[, which(colnames(snorm) != "GEODES014.nonrRNA" & colnames(snorm) != "GEODES033.nonrRNA")]

spark_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)

spark_key$Taxonomy <- gsub("Bacteria;", "", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("Eukaryota;", "", spark_key$Taxonomy)
spark_key$Phylum <- sapply(strsplit(as.character(spark_key$Taxonomy),";"), `[`, 1)

snorm$Genes <- rownames(snorm)
snorm <- melt(snorm)
snorm$variable <- gsub(".nonrRNA", "", snorm$variable)
snorm$Timepoint <- metadata$Timepoint[match(snorm$variable, metadata$Sample)]
snorm$Taxonomy <- spark_key$Phylum[match(snorm$Genes, spark_key$Gene)]
snorm$Taxonomy <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("NO CLASSIFICATION MH", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("NO CLASSIFICATION LP", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("None", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified unclassified", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified Oligohymenophorea", "Ciliophora", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified Pelagophyceae", "Ochrophyta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Unclassified ", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("UnclassifiedIsochrysidales", "Haptophyta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("UnclassifiedUnclassified", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Candidatus Saccharibacteria", "TM7", snorm$Taxonomy)
snorm$Taxonomy <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", snorm$Taxonomy)
#Remove unclassifieds to save space
snorm <- snorm[which(snorm$Taxonomy != "Unclassified" & snorm$Taxonomy != "internal standard"),]
averaged_tax <- aggregate(value ~ Taxonomy + Timepoint, data = snorm, mean)
wide_snorm <- reshape(averaged_tax, idvar = "Taxonomy", timevar = "Timepoint", direction = "wide")
rownames(wide_snorm) <- wide_snorm$Taxonomy
wide_snorm <- wide_snorm[, 2:dim(wide_snorm)[2]]
wide_snorm <- wide_snorm[which(rowSums(wide_snorm) > 3000),]

spark_phyla <- data.frame(Taxonomy = rownames(wide_snorm), Sums = rowSums(wide_snorm))

rm(spark_key)
rm(snorm)

# Choose the intersection of the top 10 expressed and abundant
keep <- c("Firmicutes", "Heterokonta", "Armatimonadetes", "Haptophtya", "Chlamydiae", "Verrucomicrobia", "Actinobacteria", "Viruses", "Planctomycetes", "Phaeophyceae", "Cyanobacteria", "TM7", "Bacteroidetes", "Proteobacteria", "Cryptophyta")
sparkling_intersect <- spark_phyla[which(spark_phyla$Taxonomy %in% keep), ]
sparkling_intersect$metaG <- spark_metaG_phyla$value[match(sparkling_intersect$Taxonomy, spark_metaG_phyla$Phylum)]
sparkling_intersect$Taxonomy <- factor(sparkling_intersect$Taxonomy, levels = sparkling_intersect$Taxonomy[order(sparkling_intersect$metaG, decreasing = F)])
sparkling_intersect <- sparkling_intersect[order(sparkling_intersect$metaG, decreasing = F),]

sparkling_intersect$Type <- c("Bacteria", "Bacteria", "Bacteria", "Algae", "Algae", "Bacteria", "Viruses", "Bacteria", "Bacteria", "Bacteria", "Algae", "Bacteria", "Bacteria")
sparkling_intersect$metaG[which(is.na(sparkling_intersect$metaG) == T)] <- 0

p <- ggplot(sparkling_intersect, aes(x = metaG, y = Sums, color = Type)) + geom_point(size = 2.5) + geom_label_repel(aes(label = Taxonomy, color = Type), force = 5, size = 7.5) + scale_color_manual(values = c("limegreen", "royalblue", "goldenrod")) + labs(x = "Proportion of metagenomic reads assigned", y = "Transcripts/L assigned", title = "Sparkling Lake") + theme(legend.position = "none")
#export as 6x8in pdf

save_plot("C:/Users/Goose and Gander/Desktop/geodes/Plots/spark_dna_vs_rna.pdf", p, base_height = 6, base_aspect_ratio = 8/6)

rm(averaged_tax)
rm(sparkling_intersect)
rm(wide_snorm)
rm(spark_metaG_phyla)
rm(spark_phyla)

# Trout Bog
tnorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
tnorm <- tnorm[, which(colnames(tnorm) != "GEODES065.nonrRNA")]
trout_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

trout_key$Taxonomy <- gsub("Bacteria;", "", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("Eukaryota;", "", trout_key$Taxonomy)
trout_key$Phylum <- sapply(strsplit(as.character(trout_key$Taxonomy),";"), `[`, 1)

tnorm$Genes <- rownames(tnorm)
tnorm <- melt(tnorm)
tnorm$variable <- gsub(".nonrRNA", "", tnorm$variable)
tnorm$Timepoint <- metadata$Timepoint[match(tnorm$variable, metadata$Sample)]
tnorm$Taxonomy <- trout_key$Phylum[match(tnorm$Genes, trout_key$Gene)]
tnorm$Taxonomy <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("NO CLASSIFICATION MH", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("NO CLASSIFICATION LP", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("None", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified unclassified", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified Oligohymenophorea", "Ciliophora", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified Pelagophyceae", "Ochrophyta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Unclassified ", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("UnclassifiedIsochrysidales", "Haptophyta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("UnclassifiedUnclassified", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Candidatus Saccharibacteria", "TM7", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", tnorm$Taxonomy)
#Remove unclassifieds to save space
tnorm <- tnorm[which(tnorm$Taxonomy != "Unclassified" & tnorm$Taxonomy != "internal standard"),]
averaged_tax <- aggregate(value ~ Taxonomy + Timepoint, data = tnorm, mean)
wide_tnorm <- reshape(averaged_tax, idvar = "Taxonomy", timevar = "Timepoint", direction = "wide")
rownames(wide_tnorm) <- wide_tnorm$Taxonomy
wide_tnorm <- wide_tnorm[, 2:dim(wide_tnorm)[2]]
wide_tnorm <- wide_tnorm[which(rowSums(wide_tnorm) > 3000),]

trout_phyla <- data.frame(Taxonomy = rownames(wide_tnorm), Sums = rowSums(wide_tnorm))

rm(tnorm)
rm(trout_key)

# Choose the intersection of the top 10 expressed and abundant
keep <- c("Verrucomicrobia", "Armatimonadetes", "Viruses", "Chlorophyta", "Acidobacteria", "Heterokonta", "Streptophyta", "Actinobacteria", "Bacteroidetes", "Cryptophyta", "Cyanobacteria", "Deinococcus-Thermus")
trout_intersect <- trout_phyla[which(trout_phyla$Taxonomy %in% keep), ]
trout_intersect$metaG <- trout_metaG_phyla$value[match(trout_intersect$Taxonomy, trout_metaG_phyla$Phylum)]
trout_intersect$Taxonomy <- factor(trout_intersect$Taxonomy, levels = trout_intersect$Taxonomy[order(trout_intersect$metaG, decreasing = F)])
trout_intersect <- trout_intersect[order(trout_intersect$metaG, decreasing = F),]

trout_intersect$Type <- c("Bacteria", "Bacteria", "Algae", "Bacteria", "Bacteria", "Algae", "Algae", "Bacteria", "Algae", "Viruses", "Bacteria", "Bacteria")

p <- ggplot(trout_intersect, aes(x = metaG, y = Sums, color = Type)) + geom_point(size = 2.5) + geom_label_repel(aes(label = Taxonomy, color = Type), force = 15, size = 7.5) + scale_color_manual(values = c("limegreen", "royalblue", "goldenrod")) + labs(x = "Proportion of metagenomic reads assigned", y = "Transcripts/L assigned", title = "Trout Bog") + theme(legend.position = "none")
#export as 6x8in pdf

save_plot("C:/Users/Goose and Gander/Desktop/geodes/Plots/trout_dna_vs_rna.pdf", p, base_height = 6, base_aspect_ratio = 8/6)

##### DESeq differential expression testing

abun_mnorm <- mnorm[order(rowSums(mnorm), decreasing = T), ]
abun_mnorm <- abun_mnorm[1:20000,]

colnames(abun_mnorm) <- gsub(".nonrRNA", "", colnames(abun_mnorm))

input <- as.matrix(abun_mnorm)
input <- input/100
input <- round(input, digits = 0)


conditions <- metadata$Time[match(colnames(abun_mnorm), metadata$Sample)]
conditions[which(conditions == 9 | conditions == 13 | conditions == 17)] <- "day"
conditions[which(conditions == 5 | conditions == 21 | conditions == 1)] <- "night"
coldata <- data.frame(samples = colnames(abun_mnorm), conditions)

cds <- DESeqDataSetFromMatrix(countData = input,
                              colData = coldata,
                              design = ~ conditions)

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds) 


dds <- DESeq(cds)
res <- results(dds)
sig.res <- res[which(res$padj < 0.05), ]
sig.res.day <- sig.res[which(sig.res$log2FoldChange < 0), ]
sig.res.night <- sig.res[which(sig.res$log2FoldChange > 0), ]
sig.res.day.key <- mendota_key[match(rownames(sig.res.day), mendota_key$Gene),]
sig.res.night.key <- mendota_key[match(rownames(sig.res.night), mendota_key$Gene),]

# Make plots of this aggreagted by short taxa and groupd products

sig.res.day.key$Taxonomy <- gsub(";;;", "", sig.res.day.key$Taxonomy)
sig.res.day.key$Taxonomy <- gsub(";;", "", sig.res.day.key$Taxonomy)
sig.res.day.key$Taxonomy <- gsub(";$", "", sig.res.day.key$Taxonomy)
spl <- strsplit(as.character(sig.res.day.key$Taxonomy), ";")
sig.res.day.key$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")

sig.res.night.key$Taxonomy <- gsub(";;;", "", sig.res.night.key$Taxonomy)
sig.res.night.key$Taxonomy <- gsub(";;", "", sig.res.night.key$Taxonomy)
sig.res.night.key$Taxonomy <- gsub(";$", "", sig.res.night.key$Taxonomy)
spl <- strsplit(as.character(sig.res.night.key$Taxonomy), ";")
sig.res.night.key$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")

# Make column by certain key words 
sig.res.day.key$Category <- "None"
sig.res.day.key$Category[grep("photo|Photo", sig.res.day.key$Product)] <- "Photosynthesis"
sig.res.day.key$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", sig.res.day.key$Product)] <- "Rhodopsin"
sig.res.day.key$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate|ose transport", sig.res.day.key$Product)] <- "Sugar degradation"
sig.res.day.key$Category[grep("rbcL|ribulose-bisphosphate carboxylase", sig.res.day.key$Product)] <- "RuBisCO"
sig.res.day.key$Category[grep("putrescine|Putrescine|spermidine|Spermidine", sig.res.day.key$Product)] <- "Polyamines"
sig.res.day.key$Category[grep("citrate lyase|Citrate lyase", sig.res.day.key$Product)] <- "rTCA"
sig.res.day.key$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", sig.res.day.key$Product)] <- "Nitrogenase"
sig.res.day.key$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", sig.res.day.key$Product)] <- "Chitinase"
sig.res.day.key$Category[grep("glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase", sig.res.day.key$Product)] <- "Glycoside_Hydrolase"
sig.res.day.key$Category[grep("alkaline phosphatase|Alkaline phosphatase", sig.res.day.key$Product)] <- "Alkaline_phosphatase"
sig.res.day.key$Category[grep("cellulase|cellulose", sig.res.day.key$Product)] <- "Cellulase"
sig.res.day.key$Category[grep("peroxidase|peroxide|catalase", sig.res.day.key$Product)] <- "ROS"
sig.res.day.key$Category[grep("ammonia monooxygenase|methane monoxygenase", sig.res.day.key$Product)] <- "Methane/Ammonia"
sig.res.day.key$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", sig.res.day.key$Product)] <- "Nitrite_reduction"
sig.res.day.key$Category[grep("urease", sig.res.day.key$Product)] <- "Urease"
sig.res.day.key$Category[grep("protease", sig.res.day.key$Product)] <- "Protease"

sig.res.night.key$Category <- "None"
sig.res.night.key$Category[grep("photo|Photo", sig.res.night.key$Product)] <- "Photosynthesis"
sig.res.night.key$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", sig.res.night.key$Product)] <- "Rhodopsin"
sig.res.night.key$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate|ose transport", sig.res.night.key$Product)] <- "Sugar degradation"
sig.res.night.key$Category[grep("rbcL|ribulose-bisphosphate carboxylase", sig.res.night.key$Product)] <- "RuBisCO"
sig.res.night.key$Category[grep("putrescine|Putrescine|spermidine|Spermidine", sig.res.night.key$Product)] <- "Polyamines"
sig.res.night.key$Category[grep("citrate lyase|Citrate lyase", sig.res.night.key$Product)] <- "rTCA"
sig.res.night.key$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", sig.res.night.key$Product)] <- "Nitrogenase"
sig.res.night.key$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", sig.res.night.key$Product)] <- "Chitinase"
sig.res.night.key$Category[grep("glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase", sig.res.night.key$Product)] <- "Glycoside_Hydrolase"
sig.res.night.key$Category[grep("alkaline phosphatase|Alkaline phosphatase", sig.res.night.key$Product)] <- "Alkaline_phosphatase"
sig.res.night.key$Category[grep("cellulase|cellulose", sig.res.night.key$Product)] <- "Cellulase"
sig.res.night.key$Category[grep("peroxidase|peroxide|catalase", sig.res.night.key$Product)] <- "ROS"
sig.res.night.key$Category[grep("ammonia monooxygenase|methane monoxygenase", sig.res.night.key$Product)] <- "Methane/Ammonia"
sig.res.night.key$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", sig.res.night.key$Product)] <- "Nitrite_reduction"
sig.res.night.key$Category[grep("urease", sig.res.night.key$Product)] <- "Urease"
sig.res.night.key$Category[grep("protease", sig.res.night.key$Product)] <- "Protease"

sig.res.day.key$Condition <- "day"
sig.res.night.key$Condition <- "night"

sig.abun_mnorm_day <- abun_mnorm[match(sig.res.day.key$Gene, rownames(abun_mnorm)), which(conditions == "day")]
sig.res.day.key$totals <- rowSums(sig.abun_mnorm_day)
sig.abun_mnorm_night <- abun_mnorm[match(sig.res.night.key$Gene, rownames(abun_mnorm)), which(conditions == "night")]
sig.res.night.key$totals <- rowSums(sig.abun_mnorm_night)

sig.res.key <- rbind(sig.res.day.key, sig.res.night.key)
# Add in zero counts so bars are the same widths

dummy <- data.frame("Gene" = rep("nogene", 4), "Genome" = rep("nogenome", 4), "Taxonomy" = rep("notax", 4), "Product" = rep("noproduct", 4), "ShortTax" = rep("noshorttax", 4), "Category" = c("Alkaline_phosphatase", "Chitinase", "Nitrite_reduction", "Polyamines"), "Condition" = c("night", "day", "night", "day"), "totals" = c(0, 0, 0, 0))

sig.res.key <- rbind(sig.res.key, dummy)

#ggplot(data = sig.res.key[which(sig.res.key$Category != "None"), ], aes(x = Category, y = totals, fill = Condition)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "Mendota")

ggplot(data = sig.res.key[which(sig.res.key$Category != "None" & sig.res.key$Category != "Photosynthesis"), ], aes(x = Category, y = totals, fill = Condition)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "Mendota", x = NULL, y = "Transcripts/L") + scale_fill_manual(values = c("darkgoldenrod2", "royalblue4")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0, 0))

### Repeat with sparkling

abun_snorm <- snorm[order(rowSums(snorm), decreasing = T), ]
abun_snorm <- abun_snorm[1:20000,]

colnames(abun_snorm) <- gsub(".nonrRNA", "", colnames(abun_snorm))

input <- as.matrix(abun_snorm)
input <- input/100
input <- round(input, digits = 0)


conditions <- metadata$Time[match(colnames(abun_snorm), metadata$Sample)]
conditions[which(conditions == 9 | conditions == 13 | conditions == 17)] <- "day"
conditions[which(conditions == 5 | conditions == 21 | conditions == 1)] <- "night"
coldata <- data.frame(samples = colnames(abun_snorm), conditions)

cds <- DESeqDataSetFromMatrix(countData = input,
                              colData = coldata,
                              design = ~ conditions)

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds) 


dds <- DESeq(cds)
res <- results(dds)
sig.res <- res[which(res$padj < 0.05), ]
sig.res.day <- sig.res[which(sig.res$log2FoldChange < 0), ]
sig.res.night <- sig.res[which(sig.res$log2FoldChange > 0), ]
sig.res.day.key <- spark_key[match(rownames(sig.res.day), spark_key$Gene),]
sig.res.night.key <- spark_key[match(rownames(sig.res.night), spark_key$Gene),]

# Make plots of this aggreagted by short taxa and groupd products

sig.res.day.key$Taxonomy <- gsub(";;;", "", sig.res.day.key$Taxonomy)
sig.res.day.key$Taxonomy <- gsub(";;", "", sig.res.day.key$Taxonomy)
sig.res.day.key$Taxonomy <- gsub(";$", "", sig.res.day.key$Taxonomy)
spl <- strsplit(as.character(sig.res.day.key$Taxonomy), ";")
sig.res.day.key$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")

sig.res.night.key$Taxonomy <- gsub(";;;", "", sig.res.night.key$Taxonomy)
sig.res.night.key$Taxonomy <- gsub(";;", "", sig.res.night.key$Taxonomy)
sig.res.night.key$Taxonomy <- gsub(";$", "", sig.res.night.key$Taxonomy)
spl <- strsplit(as.character(sig.res.night.key$Taxonomy), ";")
sig.res.night.key$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")

# Make column by certain key words 
sig.res.day.key$Category <- "None"
sig.res.day.key$Category[grep("photo|Photo", sig.res.day.key$Product)] <- "Photosynthesis"
sig.res.day.key$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", sig.res.day.key$Product)] <- "Rhodopsin"
sig.res.day.key$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate|ose transport", sig.res.day.key$Product)] <- "Sugar degradation"
sig.res.day.key$Category[grep("rbcL|ribulose-bisphosphate carboxylase", sig.res.day.key$Product)] <- "RuBisCO"
sig.res.day.key$Category[grep("putrescine|Putrescine|spermidine|Spermidine", sig.res.day.key$Product)] <- "Polyamines"
sig.res.day.key$Category[grep("citrate lyase|Citrate lyase", sig.res.day.key$Product)] <- "rTCA"
sig.res.day.key$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", sig.res.day.key$Product)] <- "Nitrogenase"
sig.res.day.key$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", sig.res.day.key$Product)] <- "Chitinase"
sig.res.day.key$Category[grep("glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase", sig.res.day.key$Product)] <- "Glycoside_Hydrolase"
sig.res.day.key$Category[grep("alkaline phosphatase|Alkaline phosphatase", sig.res.day.key$Product)] <- "Alkaline_phosphatase"
sig.res.day.key$Category[grep("cellulase|cellulose", sig.res.day.key$Product)] <- "Cellulase"
sig.res.day.key$Category[grep("peroxidase|peroxide|catalase", sig.res.day.key$Product)] <- "ROS"
sig.res.day.key$Category[grep("ammonia monooxygenase|methane monoxygenase", sig.res.day.key$Product)] <- "Methane/Ammonia"
sig.res.day.key$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", sig.res.day.key$Product)] <- "Nitrite_reduction"
sig.res.day.key$Category[grep("urease", sig.res.day.key$Product)] <- "Urease"
sig.res.day.key$Category[grep("protease", sig.res.day.key$Product)] <- "Protease"

sig.res.night.key$Category <- "None"
sig.res.night.key$Category[grep("photo|Photo", sig.res.night.key$Product)] <- "Photosynthesis"
sig.res.night.key$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", sig.res.night.key$Product)] <- "Rhodopsin"
sig.res.night.key$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate|ose transport", sig.res.night.key$Product)] <- "Sugar degradation"
sig.res.night.key$Category[grep("rbcL|ribulose-bisphosphate carboxylase", sig.res.night.key$Product)] <- "RuBisCO"
sig.res.night.key$Category[grep("putrescine|Putrescine|spermidine|Spermidine", sig.res.night.key$Product)] <- "Polyamines"
sig.res.night.key$Category[grep("citrate lyase|Citrate lyase", sig.res.night.key$Product)] <- "rTCA"
sig.res.night.key$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", sig.res.night.key$Product)] <- "Nitrogenase"
sig.res.night.key$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", sig.res.night.key$Product)] <- "Chitinase"
sig.res.night.key$Category[grep("glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase", sig.res.night.key$Product)] <- "Glycoside_Hydrolase"
sig.res.night.key$Category[grep("alkaline phosphatase|Alkaline phosphatase", sig.res.night.key$Product)] <- "Alkaline_phosphatase"
sig.res.night.key$Category[grep("cellulase|cellulose", sig.res.night.key$Product)] <- "Cellulase"
sig.res.night.key$Category[grep("peroxidase|peroxide|catalase", sig.res.night.key$Product)] <- "ROS"
sig.res.night.key$Category[grep("ammonia monooxygenase|methane monoxygenase", sig.res.night.key$Product)] <- "Methane/Ammonia"
sig.res.night.key$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", sig.res.night.key$Product)] <- "Nitrite_reduction"
sig.res.night.key$Category[grep("urease", sig.res.night.key$Product)] <- "Urease"
sig.res.night.key$Category[grep("protease", sig.res.night.key$Product)] <- "Protease"

sig.res.day.key$Condition <- "day"
sig.res.night.key$Condition <- "night"

sig.abun_snorm_day <- abun_snorm[match(sig.res.day.key$Gene, rownames(abun_snorm)), which(conditions == "day")]
sig.res.day.key$totals <- rowSums(sig.abun_snorm_day)
sig.abun_snorm_night <- abun_snorm[match(sig.res.night.key$Gene, rownames(abun_snorm)), which(conditions == "night")]
sig.res.night.key$totals <- rowSums(sig.abun_snorm_night)

sig.res.key <- rbind(sig.res.day.key, sig.res.night.key)

dummy <- data.frame("Gene" = rep("nogene", 5), "Genome" = rep("nogenome", 5), "Taxonomy" = rep("notax", 5), "Product" = rep("noproduct", 5), "ShortTax" = rep("noshorttax", 5), "Category" = c("Nitrite_reduction", "Nitrogenase", "Protease", "rTCA", "RuBisCO"), "Condition" = c("night", "night", "night", "night", "day"), "totals" = c(0, 0, 0, 0, 0))

sig.res.key <- rbind(sig.res.key, dummy)


#ggplot(data = sig.res.key[which(sig.res.key$Category != "None"), ], aes(x = Category, y = totals, fill = Condition)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "sparkling")

ggplot(data = sig.res.key[which(sig.res.key$Category != "None" & sig.res.key$Category != "Photosynthesis"), ], aes(x = Category, y = totals, fill = Condition)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "Sparkling Lake", x = NULL, y = "Transcripts/L") + scale_fill_manual(values = c("darkgoldenrod2", "royalblue4")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0, 0))


### Repeat with troutling

abun_tnorm <- tnorm[order(rowSums(tnorm), decreasing = T), ]
abun_tnorm <- abun_tnorm[1:20000,]

colnames(abun_tnorm) <- gsub(".nonrRNA", "", colnames(abun_tnorm))

input <- as.matrix(abun_tnorm)
input <- input/500
input <- round(input, digits = 0)


conditions <- metadata$Time[match(colnames(abun_tnorm), metadata$Sample)]
conditions[which(conditions == 9 | conditions == 13 | conditions == 17)] <- "day"
conditions[which(conditions == 5 | conditions == 21 | conditions == 1)] <- "night"
coldata <- data.frame(samples = colnames(abun_tnorm), conditions)

cds <- DESeqDataSetFromMatrix(countData = input,
                              colData = coldata,
                              design = ~ conditions)

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds) 


dds <- DESeq(cds)
res <- results(dds)
sig.res <- res[which(res$padj < 0.05), ]
sig.res.day <- sig.res[which(sig.res$log2FoldChange < 0), ]
sig.res.night <- sig.res[which(sig.res$log2FoldChange > 0), ]
sig.res.day.key <- trout_key[match(rownames(sig.res.day), trout_key$Gene),]
sig.res.night.key <- trout_key[match(rownames(sig.res.night), trout_key$Gene),]

# Make plots of this aggreagted by short taxa and groupd products

sig.res.day.key$Taxonomy <- gsub(";;;", "", sig.res.day.key$Taxonomy)
sig.res.day.key$Taxonomy <- gsub(";;", "", sig.res.day.key$Taxonomy)
sig.res.day.key$Taxonomy <- gsub(";$", "", sig.res.day.key$Taxonomy)
spl <- strsplit(as.character(sig.res.day.key$Taxonomy), ";")
sig.res.day.key$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")

sig.res.night.key$Taxonomy <- gsub(";;;", "", sig.res.night.key$Taxonomy)
sig.res.night.key$Taxonomy <- gsub(";;", "", sig.res.night.key$Taxonomy)
sig.res.night.key$Taxonomy <- gsub(";$", "", sig.res.night.key$Taxonomy)
spl <- strsplit(as.character(sig.res.night.key$Taxonomy), ";")
sig.res.night.key$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")

# Make column by certain key words 
sig.res.day.key$Category <- "None"
sig.res.day.key$Category[grep("photo|Photo", sig.res.day.key$Product)] <- "Photosynthesis"
sig.res.day.key$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", sig.res.day.key$Product)] <- "Rhodopsin"
sig.res.day.key$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate|ose transport", sig.res.day.key$Product)] <- "Sugar degradation"
sig.res.day.key$Category[grep("rbcL|ribulose-bisphosphate carboxylase", sig.res.day.key$Product)] <- "RuBisCO"
sig.res.day.key$Category[grep("putrescine|Putrescine|spermidine|Spermidine", sig.res.day.key$Product)] <- "Polyamines"
sig.res.day.key$Category[grep("citrate lyase|Citrate lyase", sig.res.day.key$Product)] <- "rTCA"
sig.res.day.key$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", sig.res.day.key$Product)] <- "Nitrogenase"
sig.res.day.key$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", sig.res.day.key$Product)] <- "Chitinase"
sig.res.day.key$Category[grep("glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase", sig.res.day.key$Product)] <- "Glycoside_Hydrolase"
sig.res.day.key$Category[grep("alkaline phosphatase|Alkaline phosphatase", sig.res.day.key$Product)] <- "Alkaline_phosphatase"
sig.res.day.key$Category[grep("cellulase|cellulose", sig.res.day.key$Product)] <- "Cellulase"
sig.res.day.key$Category[grep("peroxidase|peroxide|catalase", sig.res.day.key$Product)] <- "ROS"
sig.res.day.key$Category[grep("ammonia monooxygenase|methane monoxygenase", sig.res.day.key$Product)] <- "Methane/Ammonia"
sig.res.day.key$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", sig.res.day.key$Product)] <- "Nitrite_reduction"
sig.res.day.key$Category[grep("urease", sig.res.day.key$Product)] <- "Urease"
sig.res.day.key$Category[grep("protease", sig.res.day.key$Product)] <- "Protease"

sig.res.night.key$Category <- "None"
sig.res.night.key$Category[grep("photo|Photo", sig.res.night.key$Product)] <- "Photosynthesis"
sig.res.night.key$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", sig.res.night.key$Product)] <- "Rhodopsin"
sig.res.night.key$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate|ose transport", sig.res.night.key$Product)] <- "Sugar degradation"
sig.res.night.key$Category[grep("rbcL|ribulose-bisphosphate carboxylase", sig.res.night.key$Product)] <- "RuBisCO"
sig.res.night.key$Category[grep("putrescine|Putrescine|spermidine|Spermidine", sig.res.night.key$Product)] <- "Polyamines"
sig.res.night.key$Category[grep("citrate lyase|Citrate lyase", sig.res.night.key$Product)] <- "rTCA"
sig.res.night.key$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", sig.res.night.key$Product)] <- "Nitrogenase"
sig.res.night.key$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", sig.res.night.key$Product)] <- "Chitinase"
sig.res.night.key$Category[grep("glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase", sig.res.night.key$Product)] <- "Glycoside_Hydrolase"
sig.res.night.key$Category[grep("alkaline phosphatase|Alkaline phosphatase", sig.res.night.key$Product)] <- "Alkaline_phosphatase"
sig.res.night.key$Category[grep("cellulase|cellulose", sig.res.night.key$Product)] <- "Cellulase"
sig.res.night.key$Category[grep("peroxidase|peroxide|catalase", sig.res.night.key$Product)] <- "ROS"
sig.res.night.key$Category[grep("ammonia monooxygenase|methane monoxygenase", sig.res.night.key$Product)] <- "Methane/Ammonia"
sig.res.night.key$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", sig.res.night.key$Product)] <- "Nitrite_reduction"
sig.res.night.key$Category[grep("urease", sig.res.night.key$Product)] <- "Urease"
sig.res.night.key$Category[grep("protease", sig.res.night.key$Product)] <- "Protease"

sig.res.day.key$Condition <- "day"
sig.res.night.key$Condition <- "night"

sig.abun_tnorm_day <- abun_tnorm[match(sig.res.day.key$Gene, rownames(abun_tnorm)), which(conditions == "day")]
sig.res.day.key$totals <- rowSums(sig.abun_tnorm_day)
sig.abun_tnorm_night <- abun_tnorm[match(sig.res.night.key$Gene, rownames(abun_tnorm)), which(conditions == "night")]
sig.res.night.key$totals <- rowSums(sig.abun_tnorm_night)

sig.res.key <- rbind(sig.res.day.key, sig.res.night.key)

dummy <- data.frame("Gene" = rep("nogene", 3), "Genome" = rep("nogenome", 3), "Taxonomy" = rep("notax", 3), "Product" = rep("noproduct", 3), "ShortTax" = rep("noshorttax", 3), "Category" = c("Chitinase", "Rhodopsin", "Urease"), "Condition" = c("day", "night", "day"), "totals" = c(0, 0, 0))

sig.res.key <- rbind(sig.res.key, dummy)


#ggplot(data = sig.res.key[which(sig.res.key$Category != "None"), ], aes(x = Category, y = totals, fill = Condition)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "troutling")

ggplot(data = sig.res.key[which(sig.res.key$Category != "None" & sig.res.key$Category != "Photosynthesis" & sig.res.key$Category != "RuBisCO"), ], aes(x = Category, y = totals, fill = Condition)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "Trout Bog", x = NULL, y = "Transcripts/L") + scale_fill_manual(values = c("darkgoldenrod2", "royalblue4")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0, 0))
