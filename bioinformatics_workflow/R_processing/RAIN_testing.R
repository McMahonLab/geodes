# Test RAIN with and without internal std

# Set up environment
path <- "/Users/Alex/Desktop/"
path2 <- "/Users/Alex/"

library("rain")

# Metatranscriptome data
snorm <- read.csv(paste(path, "geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)
tnorm <- read.csv(paste(path, "geodes_data_tables/Trout_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)
mnorm <- read.csv(paste(path, "geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)

snorm <- snorm[, which(colnames(snorm) != "GEODES014.nonrRNA" & colnames(snorm) != "GEODES033.nonrRNA")]
tnorm <- tnorm[, which(colnames(tnorm) != "GEODES065.nonrRNA")]
mnorm <- mnorm[, which(colnames(mnorm) != "GEODES158.nonrRNA")]

colnames(mnorm) <- gsub(".nonrRNA", "", colnames(mnorm))

mendota_key <- read.csv(paste(path, "geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", sep = ""), header = T)
spark_key <- read.csv(paste(path, "geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", sep = ""), header = T)
trout_key <- read.csv(paste(path, "geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", sep = ""), header = T)

# Sample data
metadata <- read.csv(file = paste(path2, "Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", sep = ""), header = T)

# To run RAIN, I need to transpose my dataset so that timepoints are rows and "series" (genes) are columns.
# I'll also need to provide a vector saying how many reps are in each timepoint
results <- rain(t(mnorm), deltat = 4, period = 24, measure.sequence = table(metadata$Timepoint[match(colnames(mnorm), metadata$Sample)]), verbose = T, adjp.method = "Bonferroni")


# Starting with Mendota post-DESeq
photo.genes <- as.character(sig.res.key$Gene[which(sig.res.key$Category == "Photosynthesis")])
photo.table <- abun_mnorm[match(photo.genes, rownames(abun_mnorm)),]
photo.results <- rain(t(photo.table), deltat = 4, period = 24, measure.sequence = table(metadata$Timepoint[match(colnames(photo.table), metadata$Sample)]), verbose = T, adjp.method = "Bonferroni")
# Plot trends of all significant genes + a trend line
sig.photo.table <- photo.table[which(photo.results$pVal < 0.05 & photo.results$phase == 12), ]
sig.photo.table <- as.data.frame(t(apply(sig.photo.table, 1, zscore)))

sig.photo.table$Genes <- rownames(sig.photo.table)
melt.sig.photo.table <- melt(sig.photo.table)
melt.sig.photo.table$Timepoint <- metadata$Timepoint[match(melt.sig.photo.table$variable, metadata$Sample)]
ave.sig.photo.table <- aggregate(value ~ Genes + Timepoint, melt.sig.photo.table, mean)

day.night.boxes <- data.frame(x1 = c(0, 0.53, 15.35, 24.53, 39.35), x2 = c(0.53, 15.53, 24.53, 39.53, 44), y1 = rep(min(ave.sig.photo.table$value), 5), y2 = rep(max(ave.sig.photo.table$value), 5), labels = c("", "Sunrise", "Sunset", "Sunrise", "Sunset"))

ggplot()  + geom_rect(data = day.night.boxes, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = c("lightsteelblue", "lemonchiffon", "lightsteelblue", "lemonchiffon", "lightsteelblue"), alpha = 0.5) + geom_text(data = day.night.boxes, inherit.aes = FALSE, aes(x = x1 + 2, y = y2 - 0.25, label = labels)) + geom_line(data = ave.sig.photo.table, inherit.aes = F, aes(x = Timepoint, y = value, group = Genes), color = "grey") + stat_summary(data = ave.sig.photo.table, inherit.aes = F, aes(x = Timepoint, y = value), fun.y=mean, geom="line", colour="red") + scale_x_continuous(breaks = seq(0,44, by = 4)) + labs(x = "Hours into Time Series", y = "Z-score normalized reads", title = "Lake Mendota - Photosynthesis Expression") + background_grid(major = "xy", minor = "none")

cyclic.genes <- rownames(sig.photo.table)
cyclic.genes.key <- sig.res.key[match(cyclic.genes, sig.res.key$Gene), ]


ggplot(data = cyclic.genes.key, aes(x = Condition, y = totals, fill = Phylum)) + geom_bar(stat = "identity", position = "fill")  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0, 0)) + labs(x = NULL, y = "Proportion of reads assigned", title = "Lake Mendota") + facet_wrap(~ Category) + scale_fill_manual(values = as.character(colorkey$color[which((colorkey$phylum %in% unique(sig.res.key2$Phylum)) == T)]))
