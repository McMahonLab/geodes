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
#Note: change this to the minimum needed to differentiate between replicates

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

####Calculate the mean of standard error for each timepoint for each gene

std_err <- function(x) (sd(x)/sqrt(length(x)))/mean(x) * 100

# err_snorm <- aggregate(value ~ Genes + Timepoint, data = abun_snorm, std_err)
# err_snorm <- err_snorm[which(is.na(err_snorm$value) == F),]
# err_snorm <- aggregate(value ~ Genes, data = err_snorm, mean)

# err_tnorm <- aggregate(value ~ Genes + Timepoint, data = abun_tnorm, std_err)
# err_tnorm <- err_tnorm[which(is.na(err_tnorm$value) == F),]
# err_tnorm <- aggregate(value ~ Genes, data = err_tnorm, mean)

err_mnorm <- aggregate(value ~ Genes + Timepoint, data = abun_mnorm, std_err)
err_mnorm <- err_mnorm[which(is.na(err_mnorm$value) == F),]
err_mnorm <- aggregate(value ~ Genes, data = err_mnorm, mean)

# Then run the significance test for cyclic trends
fdr.mendota <- fdrtool(fisher.g.test(t(new_abun_mnorm)), statistic = "pvalue")
sig.mendota <- t(new_abun_mnorm[which(fdr.mendota$pval < 0.05),])
sig.mendota <- sig.mendota[, which(err_mnorm$value[match(colnames(sig.mendota), err_mnorm$Genes)] < 50)]

# fdr.spark <- fdrtool(fisher.g.test(t(new_abun_snorm)), statistic = "pvalue")
# sig.spark <- t(new_abun_snorm[which(fdr.spark$pval < 0.05),])
# sig.spark <- sig.spark[, which(err_snorm$value[match(colnames(sig.spark), err_snorm$Genes)] < 50)]

# fdr.trout <- fdrtool(fisher.g.test(t(new_abun_tnorm)), statistic = "pvalue")
# sig.trout <- t(new_abun_tnorm[which(fdr.trout$pval < 0.05),])
# sig.trout <- sig.trout[, which(err_tnorm$value[match(colnames(sig.trout), err_tnorm$Genes)] < 50)]

### WGCNA
# 1st, check the soft threshold in case removing samples changed it.

mpowers <- c()
for(i in 1:100){
  lottery <- sample(1:dim(sig.mendota)[2], 2000, replace = F)
  sft = pickSoftThreshold(sig.mendota[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
  mpowers[i] <- sft$powerEstimate
}
hist(mpowers)
mean(mpowers)
median(mpowers)
# Mendota Soft threshold = 8 in this case


# spowers <- c()
# for(i in 1:100){
#   lottery <- sample(1:dim(sig.spark)[2], 2000, replace = F)
#   sft = pickSoftThreshold(sig.spark[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
#   spowers[i] <- sft$powerEstimate
# }
# hist(spowers)
# mean(spowers)
# median(spowers)

# Sparkling soft threshold = 7

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

sig.mendota.key <- mendota_key[match(colnames(sig.mendota), mendota_key$Gene), ]
sig.mendota.key$Cluster <- mendota_net$colors
sig.mendota.key$Totals <- rowSums(new_abun_mnorm)[match(sig.mendota.key$Gene, rownames(new_abun_mnorm))]

# sig.spark.key <- spark_key[match(colnames(sig.spark), spark_key$Gene), ]
# sig.spark.key$Cluster <- spark_net$colors
# sig.spark.key$Totals <- rowSums(new_abun_snorm)[match(sig.spark.key$Gene, rownames(new_abun_snorm))]

# sig.trout.key <- trout_key[match(colnames(sig.trout), trout_key$Gene), ]
# sig.trout.key$Cluster <- trout_net$colors
# sig.trout.key$Totals <- rowSums(new_abun_tnorm)[match(sig.trout.key$Gene, rownames(new_abun_tnorm))]


# Save output
write.csv(sig.mendota.key, "D:/geodes_data_tables/WGCNA_mendota_results_2018-03-09.csv", row.names = F)
write.csv(mendota_net$MEs, "D:/geodes_data_tables/WGCNA_mendota_eigenvectors_2018-03-09.csv", row.names = T)

# write.csv(sig.spark.key, "D:/geodes_data_tables/WGCNA_sparkling_results_2018-03-09.csv", row.names = F)
# write.csv(spark_net$MEs, "D:/geodes_data_tables/WGCNA_sparkling_eigenvectors_2018-03-09.csv", row.names = T)

# write.csv(sig.trout.key, "D:/geodes_data_tables/WGCNA_trout_results_2018-03-09.csv", row.names = F)
# write.csv(trout_net$MEs, "D:/geodes_data_tables/WGCNA_trout_eigenvectors_2018-03-09.csv", row.names = T)


####The opposite: which genes are constant?

spark_constant <- is.constant(t(new_abun_snorm))
spark_constant_table <- new_abun_snorm[which(spark_constant == T),]
#NONE

trout_constant <- is.constant(t(new_abun_tnorm))
trout_constant_table <- new_abun_tnorm[which(trout_constant == T),]
write.csv(trout_constant_table, "D:/geodes_data_tables/trout_constant_genes_2018-03-09.csv", row.names = T)
#3564

mendota_constant <- is.constant(t(new_abun_mnorm))
mendota_constant_table <- new_abun_mnorm[which(mendota_constant == T),]
#NONE