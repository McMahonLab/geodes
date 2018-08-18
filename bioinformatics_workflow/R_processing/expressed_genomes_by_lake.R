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
readcounts <- readcounts[, which(colnames(readcounts) != "GEODES014" & colnames(readcounts) != "GEODES033" & colnames(readcounts) != "GEODES065" & colnames(readcounts) != "GEODES158")]

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
