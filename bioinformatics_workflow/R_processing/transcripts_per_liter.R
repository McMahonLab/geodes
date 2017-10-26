# Normalize table to transcripts/L

readCounts <- read.csv("D:/GEODES2refMAGsSAGs_readcounts.txt", header = T, row.names = 1)
sample_data <- read.csv("C:/Users/Alex/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv")

# Remove samples with poor std
std <- readCounts[which(rownames(readCounts) == "pFN18A_DNA_transcript"),]
readCounts2 <- readCounts[1:248805, which(std > 50)]
# Calculate standard normalization factor

std2 <- readCounts2[which(rownames(readCounts2) == "pFN18A_DNA_transcript"),]
std_factor <- std2/614000000

# Calculate volume factor

new_sample_data <- sample_data[match(substr(start = 1, stop = 9, colnames(readCounts2)), sample_data$Sample[1:109]), ]
vol_factor <- new_sample_data$Vol_Filtered

norm_val <- std_factor * vol_factor
# Normalize to transcripts/L

for(i in 1:length(norm_val)){
  readCounts2[,i] <- readCounts2[,i]/as.numeric(norm_val[i])
}

# Save normalized table
write.csv(readCounts2, "D:/GEODES2refMAGsSAGs_normalized_readcounts.csv", quote = F)


# Make a table summed by genome
all_genomes <- substr(rownames(readCounts2), start = 1, stop = 10)
genomes <- unique(all_genomes)
genome_table <- readCounts2[1,]

for(i in 1:length(genomes)){
  genes <- readCounts2[which(all_genomes == genomes[i]), ]
  genome_row <- colSums(genes)
  genome_table <- rbind(genome_table, genome_row)
}

genome_table <- genome_table[2:267,]
rownames(genome_table) <- genomes

# save genome table

write.csv(genome_table, "D:/GEODES2refMAGsSAGs_normalized_genomecounts.csv", quote = F)

