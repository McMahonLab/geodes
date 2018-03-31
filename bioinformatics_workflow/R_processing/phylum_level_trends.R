# Plot phyla aggregated trends over a day/night cycle
library(ggplot2)
library(cowplot)
library(reshape2)

### Load data (start with only one to save RAM and comment the rest out)
# Normalized read tables
snorm <- read.csv("D:/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
tnorm <- read.csv("D:/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mnorm <- read.csv("D:/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

# Gene keys
mendota_key <- read.csv("D:/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
spark_key <- read.csv("D:/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
trout_key <- read.csv("D:/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Alex/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# Filter by phyla level data

mendota_key$Taxonomy <- gsub("Bacteria;", "", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("Eukaryota;", "", mendota_key$Taxonomy)
mendota_key$Phylum <- sapply(strsplit(as.character(mendota_key$Taxonomy),";"), `[`, 1)
mendota_key$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("None", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified unclassified", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", mendota_key$Phylum)
mendota_key$Phylum <- gsub("unclassified", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("Unclassified ", "Unclassified", mendota_key$Phylum)
mendota_key$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", mendota_key$Phylum)
mendota_key <- mendota_key[which(mendota_key$Phylum != "Unclassified"),]

spark_key$Taxonomy <- gsub("Bacteria;", "", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("Eukaryota;", "", spark_key$Taxonomy)
spark_key$Phylum <- sapply(strsplit(as.character(spark_key$Taxonomy),";"), `[`, 1)
spark_key$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", spark_key$Phylum)
spark_key$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", spark_key$Phylum)
spark_key$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", spark_key$Phylum)
spark_key$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", spark_key$Phylum)
spark_key$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", spark_key$Phylum)
spark_key$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("None", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified unclassified", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", spark_key$Phylum)
spark_key$Phylum <- gsub("unclassified", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("Unclassified ", "Unclassified", spark_key$Phylum)
spark_key$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", spark_key$Phylum)
spark_key <- spark_key[which(spark_key$Phylum != "Unclassified"),]

trout_key$Taxonomy <- gsub("Bacteria;", "", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("Eukaryota;", "", trout_key$Taxonomy)
trout_key$Phylum <- sapply(strsplit(as.character(trout_key$Taxonomy),";"), `[`, 1)
trout_key$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", trout_key$Phylum)
trout_key$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", trout_key$Phylum)
trout_key$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", trout_key$Phylum)
trout_key$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", trout_key$Phylum)
trout_key$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", trout_key$Phylum)
trout_key$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("None", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified unclassified", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", trout_key$Phylum)
trout_key$Phylum <- gsub("unclassified", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("Unclassified ", "Unclassified", trout_key$Phylum)
trout_key$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", trout_key$Phylum)
trout_key <- trout_key[which(trout_key$Phylum != "Unclassified"),]

keep <- match(rownames(mnorm), mendota_key$Gene)
mnorm <- mnorm[which(is.na(keep) == F),]

keep <- match(rownames(snorm), spark_key$Gene)
snorm <- snorm[which(is.na(keep) == F),]

keep <- match(rownames(tnorm), trout_key$Gene)
tnorm <- tnorm[which(is.na(keep) == F),]

# Melt data tables and add designation for phylum and time of daty
mnorm$Genes <- rownames(mnorm)
mnorm <- melt(mnorm)
mnorm$variable <- gsub(".nonrRNA", "", mnorm$variable)
mnorm$Time <- metadata$Time[match(mnorm$variable, metadata$Sample)]
mnorm$Phylum <- mendota_key$Phylum[match(mnorm$Genes, mendota_key$Gene)]

snorm$Genes <- rownames(snorm)
snorm <- melt(snorm)
snorm$variable <- gsub(".nonrRNA", "", snorm$variable)
snorm$Time <- metadata$Time[match(snorm$variable, metadata$Sample)]
snorm$Phylum <- spark_key$Phylum[match(snorm$Genes, spark_key$Gene)]

tnorm$Genes <- rownames(tnorm)
tnorm <- melt(tnorm)
tnorm$variable <- gsub(".nonrRNA", "", tnorm$variable)
tnorm$Time <- metadata$Time[match(tnorm$variable, metadata$Sample)]
tnorm$Phylum <- trout_key$Phylum[match(tnorm$Genes, trout_key$Gene)]


#Aggregate by time and phylum

agg_mnorm <- aggregate(value ~ Time + Phylum, mnorm, mean)
agg_snorm <- aggregate(value ~ Time + Phylum, snorm, mean)
agg_tnorm <- aggregate(value ~ Time + Phylum, tnorm, mean)

agg_mnorm$Time <- factor(agg_mnorm$Time, levels = c("5", "9", "13", "17", "21", "1"))
agg_snorm$Time <- factor(agg_snorm$Time, levels = c("5", "9", "13", "17", "21", "1"))
agg_tnorm$Time <- factor(agg_tnorm$Time, levels = c("5", "9", "13", "17", "21", "1"))

#Change to fold change
phyla <- unique(agg_mnorm$Phylum)
for(i in 1:length(phyla)){
  subset <- agg_mnorm[which(agg_mnorm$Phylum == phyla[i]), ]
  ave <- mean(subset$value)
  subset$value <- subset$value - ave
  maxrange <- max(abs(subset$value))
  subset$value <- subset$value * (10/maxrange)
  agg_mnorm[which(agg_mnorm$Phylum == phyla[i]), ] <- subset
}

phyla <- unique(agg_tnorm$Phylum)
for(i in 1:length(phyla)){
  subset <- agg_tnorm[which(agg_tnorm$Phylum == phyla[i]), ]
  ave <- mean(subset$value)
  subset$value <- subset$value - ave
  maxrange <- max(abs(subset$value))
  subset$value <- subset$value * (10/maxrange)
  agg_tnorm[which(agg_tnorm$Phylum == phyla[i]), ] <- subset
}

phyla <- unique(agg_snorm$Phylum)
for(i in 1:length(phyla)){
  subset <- agg_snorm[which(agg_snorm$Phylum == phyla[i]), ]
  ave <- mean(subset$value)
  subset$value <- subset$value - ave
  maxrange <- max(abs(subset$value))
  subset$value <- subset$value * (10/maxrange)
  agg_snorm[which(agg_snorm$Phylum == phyla[i]), ] <- subset
}

# Pick the most expressed phyla to plot
me_expressed <- c("Chloroflexi", "Heterokonta", "Cyanobacteria", "Bacteroidetes", "Gemmatimonadetes", "Cryptophyta", "Crenarchaeaota", "Actinobacteria", "Viruses", "Proteobacteria")
agg_mnorm <- agg_mnorm[which(agg_mnorm$Phylum %in% me_expressed),]
agg_mnorm$Phylum <- factor(agg_mnorm$Phylum, levels = rev(me_expressed))

sp_expressed <- c("Heterokonta", "Cyanobacteria", "TM7", "Armatimonadetes", "Bacteroidetes", "Proteobacteria", "Cryptophyta",  "Haptophyta", "Actinobacteria", "Elusimicrobia")
agg_snorm <- agg_snorm[which(agg_snorm$Phylum %in% sp_expressed),]
agg_snorm$Phylum <- factor(agg_snorm$Phylum, levels = rev(sp_expressed))

tb_expressed <- c("Cyanobacteria", "Armatimonadetes", "Streptophyta", "Cryptophyta", "Verrucomicrobia", "Deinococcus-Thermus", "Bacteroidetes", "Actinobacteria", "Proteobacteria", "Chlorophyta")
agg_tnorm <- agg_tnorm[which(agg_tnorm$Phylum %in% tb_expressed),]
agg_tnorm$Phylum <- factor(agg_tnorm$Phylum, levels = rev(tb_expressed))

ggplot(agg_mnorm, aes(x = Time, y = Phylum, fill = value)) + geom_tile() + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow")

ggplot(agg_snorm, aes(x = Time, y = Phylum, fill = value)) + geom_tile() + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow")

ggplot(agg_tnorm, aes(x = Time, y = Phylum, fill = value)) + geom_tile() + scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow")
