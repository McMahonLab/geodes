# Visualize the phylum level classifications of the contig classifications
library(ggplot2)
library(cowplot)

GEODES005 <- read.table("G:/geodes_MG/phylodist_classifications/GEODES005.contig.classification.perc70.minhit3.txt", header = F, fill = NA, colClasses = c("character"))
GEODES005 <- GEODES005[which(GEODES005$V1 != "CLASSIFICATION"),]
GEODES005$V2[which(GEODES005$V2 == "NO")] <- "NO;unclassified"
GEODES005_phyla <- sapply(strsplit(GEODES005$V2, ";"), '[', 2)
phyla005 <- table(GEODES005_phyla)
phyla005 <- phyla005/sum(phyla005)

GEODES006 <- read.table("G:/geodes_MG/phylodist_classifications/GEODES006.contig.classification.perc70.minhit3.txt", header = F, fill = NA, colClasses = c("character"))
GEODES006 <- GEODES006[which(GEODES006$V1 != "CLASSIFICATION"),]
GEODES006$V2[which(GEODES006$V2 == "NO")] <- "NO;unclassified"
GEODES006_phyla <- sapply(strsplit(GEODES006$V2, ";"), '[', 2)
phyla006 <- table(GEODES006_phyla)
phyla006 <- phyla006/sum(phyla006)

GEODES057 <- read.table("G:/geodes_MG/phylodist_classifications/GEODES057.contig.classification.perc70.minhit3.txt", header = F, fill = NA, colClasses = c("character"))
GEODES057 <- GEODES057[which(GEODES057$V1 != "CLASSIFICATION"),]
GEODES057$V2[which(GEODES057$V2 == "NO")] <- "NO;unclassified"
GEODES057_phyla <- sapply(strsplit(GEODES057$V2, ";"), '[', 2)
phyla057 <- table(GEODES057_phyla)
phyla057 <- phyla057/sum(phyla057)

GEODES058 <- read.table("G:/geodes_MG/phylodist_classifications/GEODES058.contig.classification.perc70.minhit3.txt", header = F, fill = NA, colClasses = c("character"))
GEODES058 <- GEODES058[which(GEODES058$V1 != "CLASSIFICATION"),]
GEODES058$V2[which(GEODES058$V2 == "NO")] <- "NO;unclassified"
GEODES058_phyla <- sapply(strsplit(GEODES058$V2, ";"), '[', 2)
phyla058 <- table(GEODES058_phyla)
phyla058 <- phyla058/sum(phyla058)

GEODES117 <- read.table("G:/geodes_MG/phylodist_classifications/GEODES117.contig.classification.perc70.minhit3.txt", header = F, fill = NA, colClasses = c("character"))
GEODES117 <- GEODES117[which(GEODES117$V1 != "CLASSIFICATION"),]
GEODES117$V2[which(GEODES117$V2 == "NO")] <- "NO;unclassified"
GEODES117_phyla <- sapply(strsplit(GEODES117$V2, ";"), '[', 2)
phyla117 <- table(GEODES117_phyla)
phyla117 <- phyla117/sum(phyla117)

GEODES118 <- read.table("G:/geodes_MG/phylodist_classifications/GEODES118.contig.classification.perc70.minhit3.txt", header = F, fill = NA, colClasses = c("character"))
GEODES118 <- GEODES118[which(GEODES118$V1 != "CLASSIFICATION"),]
GEODES118$V2[which(GEODES118$V2 == "NO")] <- "NO;unclassified"
GEODES118_phyla <- sapply(strsplit(GEODES118$V2, ";"), '[', 2)
phyla118 <- table(GEODES118_phyla)
phyla118 <- phyla118/sum(phyla118)

GEODES165 <- read.table("G:/geodes_MG/phylodist_classifications/GEODES165.contig.classification.perc70.minhit3.txt", header = F, fill = NA, colClasses = c("character"))
GEODES165 <- GEODES165[which(GEODES165$V1 != "CLASSIFICATION"),]
GEODES165$V2[which(GEODES165$V2 == "NO")] <- "NO;unclassified"
GEODES165_phyla <- sapply(strsplit(GEODES165$V2, ";"), '[', 2)
phyla165 <- table(GEODES165_phyla)
phyla165 <- phyla165/sum(phyla165)

GEODES166 <- read.table("G:/geodes_MG/phylodist_classifications/GEODES166.contig.classification.perc70.minhit3.txt", header = F, fill = NA, colClasses = c("character"))
GEODES166 <- GEODES166[which(GEODES166$V1 != "CLASSIFICATION"),]
GEODES166$V2[which(GEODES166$V2 == "NO")] <- "NO;unclassified"
GEODES166_phyla <- sapply(strsplit(GEODES166$V2, ";"), '[', 2)
phyla166 <- table(GEODES166_phyla)
phyla166 <- phyla166/sum(phyla166)

GEODES167 <- read.table("G:/geodes_MG/phylodist_classifications/GEODES167.contig.classification.perc70.minhit3.txt", header = F, fill = NA, colClasses = c("character"))
GEODES167 <- GEODES167[which(GEODES167$V1 != "CLASSIFICATION"),]
GEODES167$V2[which(GEODES167$V2 == "NO")] <- "NO;unclassified"
GEODES167_phyla <- sapply(strsplit(GEODES167$V2, ";"), '[', 2)
phyla167 <- table(GEODES167_phyla)
phyla167 <- phyla167/sum(phyla167)

GEODES168 <- read.table("G:/geodes_MG/phylodist_classifications/GEODES168.contig.classification.perc70.minhit3.txt", header = F, fill = NA, colClasses = c("character"))
GEODES168 <- GEODES168[which(GEODES168$V1 != "CLASSIFICATION"),]
GEODES168$V2[which(GEODES168$V2 == "NO")] <- "NO;unclassified"
GEODES168_phyla <- sapply(strsplit(GEODES168$V2, ";"), '[', 2)
phyla168 <- table(GEODES168_phyla)
phyla168 <- phyla168/sum(phyla168)

phyla_names <- c(names(phyla005), names(phyla006), names(phyla057), names(phyla058), names(phyla117), names(phyla118), names(phyla165), names(phyla166), names(phyla167), names(phyla168))
phyla_counts <- c(phyla005, phyla006, phyla057, phyla058, phyla117, phyla118, phyla165, phyla166, phyla167, phyla168)
sample_names <- c(rep("GEODES005", length(phyla005)), rep("GEODES006", length(phyla006)), rep("GEODES057", length(phyla057)), rep("GEODES058", length(phyla058)), rep("GEODES117", length(phyla117)), rep("GEODES118", length(phyla118)), rep("GEODES165", length(phyla165)), rep("GEODES166", length(phyla166)), rep("GEODES167", length(phyla167)), rep("GEODES168", length(phyla168)))
lake <- c()
lake[which(sample_names == "GEODES005" | sample_names == "GEODES006"| sample_names == "GEODES165"| sample_names == "GEODES166"| sample_names == "GEODES167"| sample_names == "GEODES168")] <- "Sparkling"
lake[which(sample_names == "GEODES057" | sample_names == "GEODES058")] <- "Trout"
lake[which(sample_names == "GEODES117" | sample_names == "GEODES118")] <- "Mendota"


all_data <- data.frame(phyla_names, phyla_counts, sample_names, lake)
colnames(all_data) <- c("Phylum", "Contigs", "Metagenome", "Lake")
all_data <- all_data[which(all_data$Phylum != "unclassified"),]
all_data$Contigs <- all_data$Contigs * 100

ggplot(all_data, aes(x = Phylum, y = Metagenome, color = Lake, size = Contigs)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  + background_grid(major = "x", minor = "none") + scale_color_manual(values=c("#35978f", "#74add1", "#bf812d")) + labs(title = "Contig Classifications")
