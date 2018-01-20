### Analyze expression by taxonomic group

### Load packages
library(ggplot2)
library(cowplot)
library(reshape2)
library(GeneCycle)

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

### Pull out genes that have a classification

genes2keep <- mendota_key[which(mendota_key$Taxonomy != "NO CLASSIFICATION MH" & mendota_key$Taxonomy != "NO CLASSIFICATION LP" & mendota_key$Taxonomy != "None" & mendota_key$Taxonomy != "Bacteria;"),]
#keep only genes classified to at least the lineage/genus level
lineage1 <- sapply(strsplit(as.character(genes2keep$Taxonomy),";"), `[`, 4)
lineage2 <- sapply(strsplit(as.character(genes2keep$Taxonomy),","), `[`, 4)

lineage1[which(is.na(lineage1) == T)] <- lineage2[which(is.na(lineage1) == T)]
lineage1[which(lineage1 == "")] <- lineage2[which(lineage1 == "")]
genes2keep$Lineage <- lineage1
genes2keep <- genes2keep[grep("unclassified", lineage1, invert = T), ]

classified_mnorm <- mnorm[match(genes2keep$Gene, rownames(mnorm)),]

#remove original datasets to save RAM
rm(mnorm)
rm(mendota_key)

### Aggregate read counts by sample and classification
#get the lowest classification for each gene
tribe <- sapply(strsplit(as.character(genes2keep$Taxonomy),";"), `[`, 6)
clade <- sapply(strsplit(as.character(genes2keep$Taxonomy),";"), `[`, 5)
lineage <- sapply(strsplit(as.character(genes2keep$Taxonomy),";"), `[`, 4)
algae <- sapply(strsplit(as.character(genes2keep$Taxonomy),","), `[`, 4)

lowest_tax <- tribe
lowest_tax[which(lowest_tax == "" | is.na(lowest_tax) == T)] <- clade[which(lowest_tax == "" | is.na(lowest_tax) == T)]
lowest_tax[which(lowest_tax == "" | is.na(lowest_tax) == T)] <- lineage[which(lowest_tax == "" | is.na(lowest_tax) == T)]
lowest_tax[which(lowest_tax == "" | is.na(lowest_tax) == T)] <- algae[which(lowest_tax == "" | is.na(lowest_tax) == T)]
genes2keep$Lowest_Tax <- lowest_tax

classified_mnorm$Genes <- rownames(classified_mnorm)
classified_mnorm <- melt(classified_mnorm)
classified_mnorm$Timepoint <- metadata$Timepoint[match(classified_mnorm$variable, metadata$Sample)]
classified_mnorm$Taxonomy <- genes2keep$Lowest_Tax[match(classified_mnorm$Genes, genes2keep$Gene)]
averaged <- aggregate(value ~ Genes + Timepoint, data = classified_mnorm, mean)
long_mnorm <- reshape(averaged, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(long_mnorm) <- long_mnorm$Genes
long_mnorm <- long_mnorm[, 2:dim(long_mnorm)[2]]
long_mnorm <- long_mnorm[which(rowSums(long_mnorm) > 3000),]

### For each gene, I want to know if it has a significant diel trend with a 3, 6, or 12 interval, and if so, when its max is
plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  #plot(plot.data, t="h", lwd=2, main="", 
       #xlab="Frequency (Hz)", ylab="Strength", 
       #xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
  
  return(plot.data)
}

fdr.mendota <- fdrtool(fisher.g.test(t(long_mnorm)), statistic = "pvalue")
sig.mnorm <- long_mnorm[which(fdr.mendota$pval < 0.05), ]

genes <- rownames(sig.mnorm)
interval <- c()
maxval <- c()
maxpoint <- c()

for(i in 1:dim(sig.mnorm)[1]){
  input_data <- as.numeric(sig.mnorm[i,])
  maxval[i] <- max(input_data)
  maxpoint[i] <- which(input_data == max(input_data))
  fourier <- fft(input_data)
  s <- plot.frequency.spectrum(fourier)
  int <- s[2:12,1][which(s[2:12,2] == max(s[2:12,2]))]*4
  interval[i] <- int[1]
}

taxa_plot <- data.frame(genes, interval, maxval, maxpoint)
taxa_plot$Taxonomy <- genes2keep$Lowest_Tax[match(taxa_plot$genes, genes2keep$Gene)]
taxa_plot$Product <- genes2keep$Product[match(taxa_plot$genes, genes2keep$Gene)]
taxa_plot <- taxa_plot[which(is.na(taxa_plot$Taxonomy) == F),]

ggplot(taxa_plot, aes(x = (maxpoint-1)*4, y = Taxonomy, size = maxval)) + geom_point()

hist(taxa_plot$maxpoint[which(taxa_plot$Taxonomy == "acI-B1")])
taxa_plot$Product[which(taxa_plot$Taxonomy == "acI-B1" & taxa_plot$maxpoint <= 3)]
taxa_plot$Product[which(taxa_plot$Taxonomy == "acI-B1" & taxa_plot$maxpoint > 3 & taxa_plot$maxpoint <= 6 )]
