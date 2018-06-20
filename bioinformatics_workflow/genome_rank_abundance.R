### Rank abundance of genomes by expression
### Expression should be relative to their size
### Bins count as genomes
library(ggplot2)
library(cowplot)

### Read in data
# metatranscriptome data:
snorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
tnorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mnorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

snorm <- snorm[, which(colnames(snorm) != "GEODES014.nonrRNA" & colnames(snorm) != "GEODES033.nonrRNA")]
tnorm <- tnorm[, which(colnames(tnorm) != "GEODES065.nonrRNA")]
mnorm <- mnorm[, which(colnames(mnorm) != "GEODES158.nonrRNA")]

mendota_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T, colClasses = c("character"))
spark_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T, colClasses = c("character"))
trout_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T, colClasses = c("character"))

# bin data:
bininfo <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/GEODES_bin_data.csv", header = T)
binned_contigs <- read.table("C:/Users/Goose and Gander/Documents/geodes_data_tables/GEODES_binned_contigs.txt", header = F, colClasses = c("character"))

# GEODES SAG list:
sags_mags <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/GEODES_MAGS_SAGS_info.csv", header = T, sep = '\t', colClasses = c("character", "character", "numeric", "character", "character", "character" , "character", "character", "character", "character"))
sags_mags$Taxonomy <- apply(sags_mags[,5:10], 1, paste, collapse=";")

# Algae info
algae <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/refseq_algae_lengths.csv", header = T, sep = ' ')

# Combine genome info
red.bininfo <- bininfo[,c(1,3)]
red.bininfo$Type <- "Bin"
red.sags_mags <- sags_mags[,c(1,3,4)]
colnames(red.bininfo) <- colnames(red.sags_mags) <- colnames(algae)

genome_info <- rbind(red.bininfo, red.sags_mags, algae)

# Step 1. Replace contig names with bin names where applicable. Remove unbinned contigs.

mendota_index <- which(mendota_key$Genome %in% binned_contigs$V1)
mendota_to_replace <- mendota_key$Genome[mendota_index]
mendota_bins <- binned_contigs$V2[match(mendota_to_replace, binned_contigs$V1)]
mendota_key$Genome[mendota_index] <- mendota_bins
remove <- substr(mendota_key$Genome, start = 1, stop = 2)
mendota_key <- mendota_key[which(remove != "Ga"), ]

trout_index <- which(trout_key$Genome %in% binned_contigs$V1)
trout_to_replace <- trout_key$Genome[trout_index]
trout_bins <- binned_contigs$V2[match(trout_to_replace, binned_contigs$V1)]
trout_key$Genome[trout_index] <- trout_bins
remove <- substr(trout_key$Genome, start = 1, stop = 2)
trout_key <- trout_key[which(remove != "Ga"), ]

spark_index <- which(spark_key$Genome %in% binned_contigs$V1)
spark_to_replace <- spark_key$Genome[spark_index]
spark_bins <- binned_contigs$V2[match(spark_to_replace, binned_contigs$V1)]
spark_key$Genome[spark_index] <- spark_bins
remove <- substr(spark_key$Genome, start = 1, stop = 2)
spark_key <- spark_key[which(remove != "Ga"), ]

# Step 2. Sum expression by genome
mnorm <- mnorm[match(mendota_key$Gene, rownames(mnorm)), ]
mendota_sums <- data.frame(Gene = mendota_key$Gene, Genome = mendota_key$Genome, Totals = rowSums(mnorm))
mendota_genome_sums <- aggregate(Totals ~ Genome, mendota_sums, sum)

tnorm <- tnorm[match(trout_key$Gene, rownames(tnorm)), ]
trout_sums <- data.frame(Gene = trout_key$Gene, Genome = trout_key$Genome, Totals = rowSums(tnorm))
trout_genome_sums <- aggregate(Totals ~ Genome, trout_sums, sum)

snorm <- snorm[match(spark_key$Gene, rownames(snorm)), ]
spark_sums <- data.frame(Gene = spark_key$Gene, Genome = spark_key$Genome, Totals = rowSums(snorm))
spark_genome_sums <- aggregate(Totals ~ Genome, spark_sums, sum)

# Clean up unnecessary files
rm(mnorm)
rm(snorm)
rm(tnorm)
rm(mendota_key)
rm(spark_key)
rm(trout_key)
rm(binned_contigs)
rm(red.bininfo)
rm(red.sags_mags)
rm(mendota_sums)
rm(spark_sums)
rm(trout_sums)
rm(mendota_bins)
rm(spark_bins)
rm(trout_bins)
rm(mendota_index)
rm(trout_index)
rm(spark_index)
rm(mendota_to_replace)
rm(spark_to_replace)
rm(trout_to_replace)
rm(remove)


# Step 3. Normalize expression by length
mendota_genome_sums$Length <- genome_info$Length[match(mendota_genome_sums$Genome, genome_info$Genome)]
spark_genome_sums$Length <- genome_info$Length[match(spark_genome_sums$Genome, genome_info$Genome)]
trout_genome_sums$Length <- genome_info$Length[match(trout_genome_sums$Genome, genome_info$Genome)]

mendota_genome_sums$NormTotal <- mendota_genome_sums$Totals/mendota_genome_sums$Length
spark_genome_sums$NormTotal <- spark_genome_sums$Totals/spark_genome_sums$Length
trout_genome_sums$NormTotal <- trout_genome_sums$Totals/trout_genome_sums$Length

# Step 4. Designate genome type
mendota_genome_sums$Type <- genome_info$Type[match(mendota_genome_sums$Genome, genome_info$Genome)]
spark_genome_sums$Type <- genome_info$Type[match(spark_genome_sums$Genome, genome_info$Genome)]
trout_genome_sums$Type <- genome_info$Type[match(trout_genome_sums$Genome, genome_info$Genome)]

# Step 5. Barplot of rank abundance with genome type colored

mendota_genome_sums$Genome <- factor(mendota_genome_sums$Genome, levels = mendota_genome_sums$Genome[order(mendota_genome_sums$NormTotal, decreasing = T)])
spark_genome_sums$Genome <- factor(spark_genome_sums$Genome, levels = spark_genome_sums$Genome[order(spark_genome_sums$NormTotal, decreasing = T)])
trout_genome_sums$Genome <- factor(trout_genome_sums$Genome, levels = trout_genome_sums$Genome[order(trout_genome_sums$NormTotal, decreasing = T)])

mendota_genome_sums <- mendota_genome_sums[order(mendota_genome_sums$NormTotal, decreasing = T), ]
spark_genome_sums <- spark_genome_sums[order(spark_genome_sums$NormTotal, decreasing = T), ]
trout_genome_sums <- trout_genome_sums[order(trout_genome_sums$NormTotal, decreasing = T), ]

a <- ggplot(mendota_genome_sums, aes(x = Genome, y = NormTotal, fill = Type)) + geom_bar(stat = "identity") + labs(y = "Expression/Number of Base Pairs", title = "Mendota", x = "Genome") + theme(axis.text.x = element_blank()) + scale_fill_brewer(palette = "Paired")

b <- ggplot(mendota_genome_sums[1:100,], aes(x = Genome, y = NormTotal, fill = Type)) + geom_bar(stat = "identity") + labs(y = "Expression/Number of Base Pairs", title = "Mendota Top 100", x = "Genome") + theme(axis.text.x = element_blank()) + scale_fill_brewer(palette = "Paired")

mendota <- plot_grid(a, b, nrow = 2)


a <- ggplot(trout_genome_sums, aes(x = Genome, y = NormTotal, fill = Type)) + geom_bar(stat = "identity") + labs(y = "Expression/Number of Base Pairs", title = "Trout", x = "Genome") + theme(axis.text.x = element_blank()) + scale_fill_brewer(palette = "Paired")

b <- ggplot(trout_genome_sums[1:100,], aes(x = Genome, y = NormTotal, fill = Type)) + geom_bar(stat = "identity") + labs(y = "Expression/Number of Base Pairs", title = "Trout Top 100", x = "Genome") + theme(axis.text.x = element_blank()) + scale_fill_brewer(palette = "Paired")

trout <- plot_grid(a, b, nrow = 2)


a <- ggplot(spark_genome_sums, aes(x = Genome, y = NormTotal, fill = Type)) + geom_bar(stat = "identity") + labs(y = "Expression/Number of Base Pairs", title = "Sparkling", x = "Genome") + theme(axis.text.x = element_blank()) + scale_fill_brewer(palette = "Paired")

b <- ggplot(spark_genome_sums[1:100,], aes(x = Genome, y = NormTotal, fill = Type)) + geom_bar(stat = "identity") + labs(y = "Expression/Number of Base Pairs", title = "Sparkling Top 100", x = "Genome") + theme(axis.text.x = element_blank()) + scale_fill_brewer(palette = "Paired")

spark <- plot_grid(a, b, nrow = 2)


# Trina also wanted this info as tables - do this and add taxonomy/completeness info where available

mendota_genome_sums$Taxonomy <- "Unclassified"
mendota_genome_sums$Taxonomy[which(mendota_genome_sums$Type == "Bin")] <- as.character(bininfo$phylodist_taxonomy)[match(mendota_genome_sums$Genome[which(mendota_genome_sums$Type == "Bin")], bininfo$bin)]
mendota_genome_sums$Taxonomy[which(mendota_genome_sums$Type == "Algae")] <- as.character(algae$Taxonomy)[match(mendota_genome_sums$Genome[which(mendota_genome_sums$Type == "Algae")], algae$Genome)]
mendota_genome_sums$Taxonomy[which(mendota_genome_sums$Type != "Algae" & mendota_genome_sums$Type != "Bin")] <- as.character(sags_mags$Taxonomy)[match(mendota_genome_sums$Genome[which(mendota_genome_sums$Type != "Algae" & mendota_genome_sums$Type != "Bin")], sags_mags$IMGOID)]
                                                                                                                                                       
spark_genome_sums$Taxonomy <- "Unclassified"
spark_genome_sums$Taxonomy[which(spark_genome_sums$Type == "Bin")] <- as.character(bininfo$phylodist_taxonomy)[match(spark_genome_sums$Genome[which(spark_genome_sums$Type == "Bin")], bininfo$bin)]
spark_genome_sums$Taxonomy[which(spark_genome_sums$Type == "Algae")] <- as.character(algae$Taxonomy)[match(spark_genome_sums$Genome[which(spark_genome_sums$Type == "Algae")], algae$Genome)]
spark_genome_sums$Taxonomy[which(spark_genome_sums$Type != "Algae" & spark_genome_sums$Type != "Bin")] <- as.character(sags_mags$Taxonomy)[match(spark_genome_sums$Genome[which(spark_genome_sums$Type != "Algae" & spark_genome_sums$Type != "Bin")], sags_mags$IMGOID)]

trout_genome_sums$Taxonomy <- "Unclassified"
trout_genome_sums$Taxonomy[which(trout_genome_sums$Type == "Bin")] <- as.character(bininfo$phylodist_taxonomy)[match(trout_genome_sums$Genome[which(trout_genome_sums$Type == "Bin")], bininfo$bin)]
trout_genome_sums$Taxonomy[which(trout_genome_sums$Type == "Algae")] <- as.character(algae$Taxonomy)[match(trout_genome_sums$Genome[which(trout_genome_sums$Type == "Algae")], algae$Genome)]
trout_genome_sums$Taxonomy[which(trout_genome_sums$Type != "Algae" & trout_genome_sums$Type != "Bin")] <- as.character(sags_mags$Taxonomy)[match(trout_genome_sums$Genome[which(trout_genome_sums$Type != "Algae" & trout_genome_sums$Type != "Bin")], sags_mags$IMGOID)]

# Output top 100 into tables

write.csv(mendota_genome_sums[1:100,], file = "C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_genomes.csv", row.names = F)
write.csv(spark_genome_sums[1:100,], file = "C:/Users/Goose and Gander/Documents/geodes_data_tables/Sparkling_genomes.csv", row.names = F)
write.csv(trout_genome_sums[1:100,], file = "C:/Users/Goose and Gander/Documents/geodes_data_tables/Trout_genomes.csv", row.names = F)
