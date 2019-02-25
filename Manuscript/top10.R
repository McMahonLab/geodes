# Set up environment
path <- "/Users/Alex/Desktop/"
path2 <- "/Users/Alex/"

library(ggplot2)
library(cowplot)
library(reshape2)
library(ggrepel)
library(raster)


# Sample data
metadata <- read.csv(file = paste(path2, "Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", sep = ""), header = T)

mnorm <- read.csv(paste(path, "geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)
mnorm <- mnorm[, which(colnames(mnorm) != "GEODES158.nonrRNA")]
colnames(mnorm) <- gsub(".nonrRNA", "", colnames(mnorm))
mendota_key <- read.csv(paste(path, "geodes_data_tables/Mendota_ID90_genekey_geneclassifications_2018-11-28.csv", sep = ""), header = T)

tnorm <- read.csv(paste(path, "geodes_data_tables/Trout_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)
tnorm <- tnorm[, which(colnames(tnorm) != "GEODES065.nonrRNA")]
colnames(tnorm) <- gsub(".nonrRNA", "", colnames(tnorm))
trout_key <- read.csv(paste(path, "geodes_data_tables/Trout_ID90_genekey_geneclassifications_2018-11-28.csv", sep = ""), header = T)

snorm <- read.csv(paste(path, "geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)
snorm <- snorm[, which(colnames(snorm) != "GEODES014.nonrRNA" & colnames(snorm) != "GEODES033.nonrRNA")]
colnames(snorm) <- gsub(".nonrRNA", "", colnames(snorm))
spark_key <- read.csv(paste(path, "geodes_data_tables/Sparkling_ID90_genekey_geneclassifications_2018-11-28.csv", sep = ""), header = T)

totals <- rowSums(snorm)
top10 <- names(totals)[order(totals, decreasing = T)]
top10 <- top10[1:11]
top10 <- spark_key[match(top10, spark_key$Gene),]
top10$totals <- totals[match(top10$Gene, names(totals))]
top10 <- top10[which(top10$Product != "internal standard"), ]
top10$Phylum[which(is.na(top10$Phylum) == T)] <- "Unclassified"
top10$Gene <- factor(top10$Gene, levels = top10$Gene[order(top10$totals, decreasing = T)])
# This one has a mislabeled contig class - the gene itself is Cyano, which makes more sense
top10$Phylum[which(top10$Gene == "Ga0164294_100668996_GEODES005")] <- "Cyanobacteria"
p1 <- ggplot(top10, aes(x = Gene, y = log10(totals))) + geom_point(size = 5, color = "grey") + theme(axis.text.x = element_blank()) + geom_text_repel(aes(label = seq(1:10)), force = 10, size = 4, color = "black") + labs(x = "Gene Rank Abundance", y = "Log Read Counts", title = "A. Sparkling Lake") + theme(legend.position = c(0.8, 0.9))

totals.2 <- rowSums(mnorm)
top10.2 <- rownames(mnorm)[order(totals.2, decreasing = T)]
top10.2 <- top10.2[1:11]
top10.2 <- mendota_key[match(top10.2, mendota_key$Gene),]
top10.2$totals.2 <- totals.2[match(top10.2$Gene, names(totals.2))]
top10.2 <- top10.2[which(top10.2$Product != "internal standard"), ]
top10.2$Phylum[which(is.na(top10.2$Phylum) == T)] <- "Unclassified"
top10.2$Gene <- factor(top10.2$Gene, levels = top10.2$Gene[order(top10.2$totals.2, decreasing = T)])
#
top10.2$Phylum[which(top10.2$Gene  == "Ga0164292_101211313_GEODES117")] <- "Unclassified"
p2 <- ggplot(top10.2, aes(x = Gene, y = log10(totals.2))) + geom_point(size = 5, color = "grey") + theme(axis.text.x = element_blank()) + geom_text_repel(aes(label = seq(1:10)), force = 10, size = 4, color = "black") + labs(x = "Gene Rank Abundance", y = "Log Read Counts", title = "B. Lake Mendota")  + theme(legend.position = c(0.8, 0.9))

totals.3 <- rowSums(tnorm)
top10.3 <- rownames(tnorm)[order(totals.3, decreasing = T)]
top10.3 <- top10.3[1:11]
top10.3 <- trout_key[match(top10.3, trout_key$Gene),]
top10.3$totals.3 <- totals.3[match(top10.3$Gene, names(totals.3))]
top10.3 <- top10.3[which(top10.3$Product != "internal standard"), ]
top10.3$Phylum[which(is.na(top10.3$Phylum) == T)] <- "Unclassified"
top10.3$Gene <- factor(top10.3$Gene, levels = top10.3$Gene[order(top10.3$totals.3, decreasing = T)])
p3 <- ggplot(top10.3, aes(x = Gene, y = log10(totals.3))) + geom_point(size = 5, color = "grey") + theme(axis.text.x = element_blank()) + geom_text_repel(aes(label = seq(1:10)), force = 10, size = 4, color = "black") + labs(x = "Gene Rank Abundance", y = "Log Read Counts", title = "C. Trout Bog") + theme(legend.position = c(0.8, 0.9))


totals <- rowSums(snorm)
top10 <- names(totals)[order(totals, decreasing = T)]
top10 <- top10[1:100]
top10 <- spark_key[match(top10, spark_key$Gene),]
top10$totals <- totals[match(top10$Gene, names(totals))]
top10 <- top10[which(top10$Product != "internal standard"), ]
top10 <- top10[grep("photo|Photo|hypothetical|ribulose-bisphosphate", top10$Product, invert = T), ]
top10$Phylum[which(is.na(top10$Phylum) == T)] <- "Unclassified"
top10 <- top10[1:10,]
top10$Gene <- factor(top10$Gene, levels = top10$Gene[order(top10$totals, decreasing = T)])
# This one has a mislabeled contig class - the gene itself is Cyano, which makes more sense
p4 <- ggplot(top10, aes(x = Gene, y = log10(totals))) + geom_point(size = 5, color = "grey") + theme(axis.text.x = element_blank()) + geom_text_repel(aes(label = seq(1:10)), force = 10, size = 4, color = "black") + labs(x = "Gene Rank Abundance", y = "Log Read Counts", title = "D. Sparkling Lake")

totals.2 <- rowSums(mnorm)
top10.2 <- names(totals.2)[order(totals.2, decreasing = T)]
top10.2 <- top10.2[1:100]
top10.2 <- mendota_key[match(top10.2, mendota_key$Gene),]
top10.2$totals.2 <- totals.2[match(top10.2$Gene, names(totals.2))]
top10.2 <- top10.2[which(top10.2$Product != "internal standard"), ]
top10.2 <- top10.2[grep("photo|Photo|hypothetical|ribulose-bisphosphate", top10.2$Product, invert = T), ]
top10.2$Phylum[which(is.na(top10$Phylum) == T)] <- "Unclassified"
top10.2 <- top10.2[1:10,]
top10.2$Gene <- factor(top10.2$Gene, levels = top10.2$Gene[order(top10.2$totals.2, decreasing = T)])
#
p5 <- ggplot(top10.2, aes(x = Gene, y = log10(totals.2))) + geom_point(size = 5, color = "grey") + theme(axis.text.x = element_blank()) + geom_text_repel(aes(label = seq(1:10)), force = 10, size = 4, color = "black") + labs(x = "Gene Rank Abundance", y = "Log Read Counts", title = "E. Lake Mendota")


totals.3 <- rowSums(tnorm)
top10.3 <- names(totals.3)[order(totals.3, decreasing = T)]
top10.3 <- top10.3[1:100]
top10.3 <- trout_key[match(top10.3, trout_key$Gene),]
top10.3$totals.3 <- totals.3[match(top10.3$Gene, names(totals.3))]
top10.3 <- top10.3[which(top10.3$Product != "internal standard"), ]
top10.3 <- top10.3[grep("photo|Photo|hypothetical|ribulose-bisphosphate", top10.3$Product, invert = T), ]
top10.3$Phylum[which(is.na(top10$Phylum) == T)] <- "Unclassified"
top10.3 <- top10.3[1:10,]
top10.3$Gene <- factor(top10.3$Gene, levels = top10.3$Gene[order(top10.3$totals.3, decreasing = T)])
p6 <- ggplot(top10.3, aes(x = Gene, y = log10(totals.3))) + geom_point(size = 5, color = "grey") + theme(axis.text.x = element_blank()) + geom_text_repel(aes(label = seq(1:10)), force = 10, size = 4, color = "black") + labs(x = "Gene Rank Abundance", y = "Log Read Counts", title = "F. Trout Bog")

part1 <- plot_grid(p1, p2, p3, nrow = 3)
part2 <- plot_grid(p4, p5, p6, nrow = 3)

save_plot("/Users/Alex/Desktop/geodes/Manuscript/figures_and_tables/top10.pdf", part1, base_height = 8, base_aspect_ratio = 0.5)
save_plot("/Users/Alex/Desktop/geodes/Manuscript/figures_and_tables/top10_hetero.pdf", part2, base_height = 8, base_aspect_ratio = 0.5)
