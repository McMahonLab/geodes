##### Find specific products in nutrient cycling

library(ggplot2)
library(cowplot)
library(reshape2)
library(GeneCycle)

#snorm <- read.csv("D:/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
#tnorm <- read.csv("D:/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mnorm <- read.csv("D:/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

# Gene keys
mendota_key <- read.csv("D:/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#spark_key <- read.csv("D:/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#trout_key <- read.csv("D:/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Alex/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# RuBisCo

rubisco_genes <- mendota_key$Gene[grep("RuBisCO|ribulose 1,5-bisphosphate carboxylase|ribulose-bisphosphate carboxylase", mendota_key$Product)]
rubisco_reads <- mnorm[match(rubisco_genes, rownames(mnorm)),]
rubisco_reads$Genes <- rownames(rubisco_reads)
rubisco_reads <- melt(rubisco_reads)
rubisco_reads$variable <- gsub(".nonrRNA", "", rubisco_reads$variable)
rubisco_reads$Timepoint <- metadata$Timepoint[match(rubisco_reads$variable, metadata$Sample)]
averaged_rubisco <- aggregate(value ~ Genes + Timepoint, data = rubisco_reads, mean)

wide_rubisco<- reshape(averaged_rubisco, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(wide_rubisco) <- wide_rubisco$Genes
wide_rubisco <- wide_rubisco[, 2:dim(wide_rubisco)[2]]
colnames(wide_rubisco) <- c("T0", "T4", "T8", "T12", "T16", "T20", "T24", "T28", "T32", "T36", "T40", "T44")

fdr.rubisco <- fdrtool(fisher.g.test(t(wide_rubisco)), statistic = "pvalue")
sig.rubisco <- wide_rubisco[which(fdr.rubisco$pval < 0.05),]

barplot(colSums(sig.rubisco))

# multiple sugar transport
sugar_genes <- mendota_key$Gene[grep("multiple sugar transport", mendota_key$Product)]
sugar_reads <- mnorm[match(sugar_genes, rownames(mnorm)),]
sugar_reads$Genes <- rownames(sugar_reads)
sugar_reads <- melt(sugar_reads)
sugar_reads$variable <- gsub(".nonrRNA", "", sugar_reads$variable)
sugar_reads$Timepoint <- metadata$Timepoint[match(sugar_reads$variable, metadata$Sample)]
averaged_sugar <- aggregate(value ~ Genes + Timepoint, data = sugar_reads, mean)

wide_sugar<- reshape(averaged_sugar, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(wide_sugar) <- wide_sugar$Genes
wide_sugar <- wide_sugar[, 2:dim(wide_sugar)[2]]
colnames(wide_sugar) <- c("T0", "T4", "T8", "T12", "T16", "T20", "T24", "T28", "T32", "T36", "T40", "T44")

fdr.sugar <- fdrtool(fisher.g.test(t(wide_sugar)), statistic = "pvalue")
sig.sugar <- wide_sugar[which(fdr.sugar$pval < 0.05),]

#Test things!
test_genes <- mendota_key$Gene[grep("galactose", mendota_key$Product)]
test_reads <- mnorm[match(test_genes, rownames(mnorm)),]
test_reads$Genes <- rownames(test_reads)
test_reads <- melt(test_reads)
test_reads$variable <- gsub(".nonrRNA", "", test_reads$variable)
test_reads$Timepoint <- metadata$Timepoint[match(test_reads$variable, metadata$Sample)]
averaged_test <- aggregate(value ~ Genes + Timepoint, data = test_reads, mean)

wide_test<- reshape(averaged_test, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(wide_test) <- wide_test$Genes
wide_test <- wide_test[, 2:dim(wide_test)[2]]
colnames(wide_test) <- c("T0", "T4", "T8", "T12", "T16", "T20", "T24", "T28", "T32", "T36", "T40", "T44")

fdr.test <- fdrtool(fisher.g.test(t(wide_test)), statistic = "pvalue")
sig.test <- wide_test[which(fdr.test$pval < 0.05),]
barplot(colSums(sig.test))
