# Overlay ANI data with expression data by lake
# This is from the "GEODES genomes" mapping, which is my semi-curated database of freshwater MAGs and SAGs

# Packages:
library(reshape2)
library(ggplot2)
library(cowplot)

# Read in data
readcounts <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/GEODES_genomes_normalized.csv", header = T, row.names = 1)
genekey <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/GEODES_genomes_genekey.csv", header = T, colClasses = c("character"))
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)
genomekey <- read.csv("C:/Users/Goose and Gander/Dropbox/genome_quality.csv", header = T, colClasses = c("character"))

# Some light data processing
# There's an extra column of outdated row numbers in genekey - remove
genekey$X <- NULL
# Swap contig names for bin names
swap <- grep("binned", genekey$Geneid)
split_names <- sapply(strsplit(genekey$Geneid,"_"), `[`, 3)
genekey$Genome[swap] <- split_names[swap]


# I know from my gene-centric analyses that there are several samples with abnormal results. Remove these:
colnames(readcounts) <- gsub(".nonrRNA", "", colnames(readcounts))
readcounts <- readcounts[, which(colnames(readcounts) != "GEODES014" & colnames(readcounts) != "GEODES033" & colnames(readcounts) != "GEODES065" & colnames(readcounts) != "GEODES158" & colnames(readcounts) != "GEODES053")]

# Write a function that takes a genome name, tells me its classification, and plots its summed expression in each lake. Also check number of genes to make sure it's not just one really high one

oracle <- function(genome){
  taxonomy <- genomekey$Classification[which(genomekey$Genome.Name == genome)]
  if(length(taxonomy) == 0){
    break
  }
  print(taxonomy)
  genes <- genekey$Geneid[which(genekey$Genome == genome)]
  if(length(genes) == 0){
    break
  }
  reads <- readcounts[match(genes, rownames(readcounts)), ]
  reads$Genes <- rownames(reads)
  reads <- melt(reads)
  reads$Lake <- metadata$Lake[match(reads$variable, metadata$Sample)]
  p <- ggplot(data = reads, aes(x = Lake, y = value)) + geom_bar(stat = "identity") + labs(title = genome)
  numgenes1 <- length(which(reads$value[which(reads$Lake == "Mendota")] > 0))
  numgenes2 <- length(which(reads$value[which(reads$Lake == "Sparkling")] > 0))
  numgenes3 <- length(which(reads$value[which(reads$Lake == "Trout")] > 0))
  print(paste("Mendota: ", numgenes1, ": Sparkling: ", numgenes2, ": Trout: ", numgenes3))
  print(p)
}

#####
# Start searching for transporters in the genomes I've identified
# My criteria: between 95-98% ANI, must be expressed in the same lake, must have at least ~500 genes expressed

# My pairs of genomes and their corresponding data:
# Betaproteobacteria: 2582580598, 2754412463 (MAG vs GEODES SAG - could be the same pop) ANI: 98%, 0.5 aligned, Methylophilales
# Correlation: 0.93

genome1 <- "GEODES057-binned.063"
#genome2 <- "2236876031"
lake <- "Trout"

# Calculate correlation over time 
lake_samples <- match(metadata$Sample[which(metadata$Lake == lake)], colnames(readcounts))
lake_samples <- lake_samples[which(is.na(lake_samples) == F)]
lake_counts <- readcounts[, lake_samples]

genome1_genekey <- genekey[which(genekey$Genome == genome1), ]
#genome2_genekey <- genekey[which(genekey$Genome == genome2), ]

genome1_counts <- lake_counts[match(genome1_genekey$Geneid, rownames(lake_counts)), ]
#genome2_counts <- lake_counts[match(genome2_genekey$Geneid, rownames(lake_counts)), ]

#cor(colSums(genome1_counts), colSums(genome2_counts))

# Trend over time
barplot(colSums(genome1_counts))


# Transporters ordered by expression

genome1_genekey$Sums <- rowSums(genome1_counts)
genome1_transporters <- genome1_genekey[grep("transport|porin|uptake|channel", genome1_genekey$Product), ]
genome1_transporters <- genome1_transporters[order(genome1_transporters$Sums, decreasing = T), ]


genome2_genekey$Sums <- rowSums(genome2_counts)
genome2_transporters <- genome2_genekey[grep("transport|porin|uptake|channel", genome2_genekey$Product), ]
genome2_transporters <- genome2_transporters[order(genome2_transporters$Sums, decreasing = T), ]

# Number of expressed transporters:
dim(genome1_transporters)[1]
dim(genome2_transporters)[1]


# Look up transporters in the sparse table to confirm that they have zero expression

readCounts <- read.table("C:/Users/Goose and Gander/Documents/geodes_data_tables/GEODES_genomes_2018-08-06.txt", header = T, row.names = 1, sep = "\t", fill = NA, quote = "")
sample_data <- read.csv("C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv")
genekey2 <- read.table("C:/Users/Goose and Gander/Documents/geodes_data_tables/genome_mapping.table.txt", sep = "\t", fill = NA, quote = "", header = T)
std <- readCounts[which(rownames(readCounts) == "pFN18A_DNA_transcript"),]
readCounts2 <- readCounts[, which(std > 50)]
std2 <- readCounts2[which(rownames(readCounts2) == "pFN18A_DNA_transcript"),]
std_factor <- std2/614000000
new_sample_data <- sample_data[match(substr(start = 1, stop = 9, colnames(readCounts2)), sample_data$Sample[1:109]), ]
vol_factor <- new_sample_data$Vol_Filtered
norm_val <- std_factor * vol_factor
for(i in 1:length(norm_val)){
  readCounts2[,i] <- readCounts2[,i]/as.numeric(norm_val[i])
}
readCounts2 <- readCounts2[, which(colnames(readCounts2) != "GEODES014" & colnames(readCounts2) != "GEODES033" & colnames(readCounts2) != "GEODES065" & colnames(readCounts2) != "GEODES158")]

genome <- "2264265094"
genomekey <- genekey2[which(genekey2$Genome == genome), ]
searchgenes <- genomekey$Geneid[grep("nitrate/sulfonate/bicarb", genomekey$Product)]
rowSums(readCounts2[match(searchgenes, rownames(readCounts2)), ])


# This isn't working particularly well. Make a rank abundance curve of all the genomes

readcounts$Genome <- sapply(strsplit(rownames(readcounts),"_"), `[`, 3)
readcounts2 <- melt(readcounts)
readcounts2$Lake <- metadata$Lake[match(readcounts2$variable, metadata$Sample)]

genomecounts <- aggregate(value ~ Genome + Lake, readcounts2, sum)
genomecounts$Phylum <- genomekey$Phylum[match(genomecounts$Genome, genomekey$Genome.Name)]

MEcounts <- genomecounts[which(genomecounts$Lake == "Mendota"), ]
SPcounts <- genomecounts[which(genomecounts$Lake == "Sparkling"), ]
TBcounts <- genomecounts[which(genomecounts$Lake == "Trout"), ]

MEcounts <- MEcounts[which(MEcounts$Phylum != "transcript"), ]
SPcounts <- SPcounts[which(SPcounts$Phylum != "transcript"), ]
TBcounts <- TBcounts[which(TBcounts$Phylum != "transcript"), ]

# Normalize reads by genome size
MEcounts$Lengths <- as.numeric(genomekey$Length[match(MEcounts$Genome, genomekey$Genome.Name)])
SPcounts$Lengths <- as.numeric(genomekey$Length[match(SPcounts$Genome, genomekey$Genome.Name)])
TBcounts$Lengths <- as.numeric(genomekey$Length[match(TBcounts$Genome, genomekey$Genome.Name)])

MEcounts$adj.value <- MEcounts$value/(MEcounts$Length/1000000)
SPcounts$adj.value <- SPcounts$value/(SPcounts$Length/1000000)
TBcounts$adj.value <- TBcounts$value/(TBcounts$Length/1000000)

# order by expression amount
MEcounts <- MEcounts[order(MEcounts$adj.value, decreasing = T), ]
SPcounts <- SPcounts[order(SPcounts$adj.value, decreasing = T), ]
TBcounts <- TBcounts[order(TBcounts$adj.value, decreasing = T), ]

MEcounts$Genome <- factor(MEcounts$Genome, levels = MEcounts$Genome)
SPcounts$Genome <- factor(SPcounts$Genome, levels = SPcounts$Genome)
TBcounts$Genome <- factor(TBcounts$Genome, levels = TBcounts$Genome)

ggplot(MEcounts[1:50,], aes(x = Genome, y = adj.value, fill = Phylum)) + geom_col() + labs(x = NULL, y = "Sum Transcripts/L per 1 million bp", title = "Mendota") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Paired")

ggplot(SPcounts[1:50,], aes(x = Genome, y = adj.value, fill = Phylum)) + geom_col() + labs(x = NULL, y = "Sum Transcripts/L per 1 million bp", title = "Sparkling") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Paired")

ggplot(TBcounts[1:50,], aes(x = Genome, y = adj.value, fill = Phylum)) + geom_col() + labs(x = NULL, y = "Sum Transcripts/L per 1 million bp", title = "Trout") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Paired")
