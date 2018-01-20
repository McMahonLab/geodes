### Summary statistics on GEODES

### Load packages
library(ggplot2)
library(cowplot)
library(reshape2)
#library(GeneCycle)

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

# How expressed is each phylum?
mendota_key$Taxonomy <- gsub("Bacteria;", "", mendota_key$Taxonomy)
mendota_key$Phylum <- sapply(strsplit(as.character(mendota_key$Taxonomy),";"), `[`, 1)

mnorm$Genes <- rownames(mnorm)
mnorm <- melt(mnorm)
mnorm$Timepoint <- metadata$Timepoint[match(mnorm$variable, metadata$Sample)]
mnorm$Taxonomy <- metadata$Phylum[match(mendota_key$Genes, mendota_key$Gene)]
averaged_tax <- aggregate(value ~ Taxonomy + Timepoint, data = mnorm, mean)
long_mnorm <- reshape(averaged, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(long_mnorm) <- long_mnorm$Genes
long_mnorm <- long_mnorm[, 2:dim(long_mnorm)[2]]
long_mnorm <- long_mnorm[which(rowSums(long_mnorm) > 3000),]


# What are the top expressed genes?

