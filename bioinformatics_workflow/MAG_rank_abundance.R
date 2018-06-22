# Check how ref MAG SAG mapping turned out

# Load data

genome_counts <- read.csv("C:/Users/Goose and Gander/Desktop/GEODES_refMAGsSAGs_mapping/GEODES2refMAGsSAGs_normalized_genomecounts.csv", header = T, row.names = 1)

genome_counts <- genome_counts[, which(colnames(genome_counts)  != "GEODES014.nonrRNA" & colnames(genome_counts)  != "GEODES033.nonrRNA" & colnames(genome_counts)  != "GEODES065.nonrRNA" & colnames(genome_counts)  != "GEODES158.nonrRNA")]

# mt_size <- read.table("C:/Users/Goose and Gander/Desktop/geodes/sample_data/MT_size.txt")
# 
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

genomeinfo <- read.csv(file = "C:/Users/Goose and Gander/Desktop/GEODES_refMAGsSAGs_mapping/genome_phylogeny.csv", header = T, row.names = 1, colClasses = c("character"))

sags_mags <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/GEODES_MAGS_SAGS_info.csv", header = T, sep = '\t', colClasses = c("character", "character", "numeric", "character", "character", "character" , "character", "character", "character", "character"))

genomeinfo$Phylum[which(genomeinfo$Phylum == "Proteobacteria")] <- genomeinfo$Class[which(genomeinfo$Phylum == "Proteobacteria")]

# read_counts <- read_counts[which(is.na(read_counts[,1]) == F), ]
# mapped_reads <- colSums(read_counts)
# 
# percent <- mapped_reads/mt_size$V2 * 100

# Rank abundance curves

# split by lake
colnames(genome_counts) <- gsub(".nonrRNA", "", colnames(genome_counts))
mendota_samples <- match(metadata$Sample[which(metadata$Lake == "Mendota")], colnames(genome_counts))
spark_samples <- match(metadata$Sample[which(metadata$Lake == "Sparkling")], colnames(genome_counts))
trout_samples <- match(metadata$Sample[which(metadata$Lake == "Trout")], colnames(genome_counts))

mendota_samples <- mendota_samples[which(is.na(mendota_samples) == F)]
trout_samples <- trout_samples[which(is.na(trout_samples) == F)]
spark_samples <- spark_samples[which(is.na(spark_samples) == F)]

mendota_counts <- genome_counts[, mendota_samples ]
spark_counts <- genome_counts[, spark_samples ]
trout_counts <- genome_counts[, trout_samples ]


mendota_plot <- data.frame("Genome" = rownames(mendota_counts), "Total" = rowSums(mendota_counts), "Phylum" = genomeinfo$Phylum[match(rownames(mendota_counts), rownames(genomeinfo))])
trout_plot <- data.frame("Genome" = rownames(trout_counts), "Total" = rowSums(trout_counts), "Phylum" = genomeinfo$Phylum[match(rownames(trout_counts), rownames(genomeinfo))])
spark_plot <- data.frame("Genome" = rownames(spark_counts), "Total" = rowSums(spark_counts), "Phylum" = genomeinfo$Phylum[match(rownames(spark_counts), rownames(genomeinfo))])

mendota_plot$Length <- sags_mags$Length[match(mendota_plot$Genome, sags_mags$IMGOID)]
trout_plot$Length <- sags_mags$Length[match(trout_plot$Genome, sags_mags$IMGOID)]
spark_plot$Length <- sags_mags$Length[match(spark_plot$Genome, sags_mags$IMGOID)]

mendota_plot$NormTotal <- mendota_plot$Total/mendota_plot$Length
trout_plot$NormTotal <- trout_plot$Total/trout_plot$Length
spark_plot$NormTotal <- spark_plot$Total/spark_plot$Length

mendota_plot <- mendota_plot[order(mendota_plot$NormTotal, decreasing = T), ]
trout_plot <- trout_plot[order(trout_plot$NormTotal, decreasing = T), ]
spark_plot <- spark_plot[order(spark_plot$NormTotal, decreasing = T), ]

mendota_plot$Genome <- factor(mendota_plot$Genome, levels = mendota_plot$Genome)
trout_plot$Genome <- factor(trout_plot$Genome, levels = trout_plot$Genome)
spark_plot$Genome <- factor(spark_plot$Genome, levels = spark_plot$Genome)

ggplot(data = mendota_plot[1:100,], aes(x = Genome, y = NormTotal, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = "Genome", y = "Coverage", title = "Mendota") + theme(axis.text.x = element_blank()) + scale_fill_brewer(palette = "Paired")
ggplot(data = trout_plot[1:100,], aes(x = Genome, y = NormTotal, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = "Genome", y = "Coverage", title = "Trout Bog") + theme(axis.text.x = element_blank()) + scale_fill_brewer(palette = "Paired")
ggplot(data = spark_plot[1:100,], aes(x = Genome, y = NormTotal, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = "Genome", y = "Coverage", title = "Sparkling Lake") + theme(axis.text.x = element_blank()) + scale_fill_brewer(palette = "Paired")

taxonomy <- apply(genomeinfo[,2:7], 1, paste, collapse=";")
mendota_plot$Taxonomy <- taxonomy[match(mendota_plot$Genome, rownames(genomeinfo))]
trout_plot$Taxonomy <- taxonomy[match(trout_plot$Genome, rownames(genomeinfo))]
spark_plot$Taxonomy <- taxonomy[match(spark_plot$Genome, rownames(genomeinfo))]

write.csv(mendota_plot[1:100,], file = "C:/Users/Goose and Gander/Documents/Mendota_rank_abundance_MAGs_SAGs.csv", row.names = F, quote = F)
write.csv(trout_plot[1:100,], file = "C:/Users/Goose and Gander/Documents/Trout_rank_abundance_MAGs_SAGs.csv", row.names = F, quote = F)
write.csv(spark_plot[1:100,], file = "C:/Users/Goose and Gander/Documents/Sparkling_rank_abundance_MAGs_SAGs.csv", row.names = F, quote = F)
