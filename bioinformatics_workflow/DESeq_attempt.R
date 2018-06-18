library(ggplot2)
library(cowplot)
library(reshape2)
library(DESeq)


snorm <- read.csv("D:/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
tnorm <- read.csv("D:/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mnorm <- read.csv("D:/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

snorm <- snorm[, which(colnames(snorm) != "GEODES014.nonrRNA" & colnames(snorm) != "GEODES033.nonrRNA")]
tnorm <- tnorm[, which(colnames(tnorm) != "GEODES065.nonrRNA")]
mnorm <- mnorm[, which(colnames(mnorm) != "GEODES158.nonrRNA")]

mendota_key <- read.csv("D:/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
spark_key <- read.csv("D:/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
trout_key <- read.csv("D:/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Alex/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# Take the top quantile of each data table
# Combine into a long table and aggregate product name by sample



