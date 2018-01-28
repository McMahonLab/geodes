### Summary statistics on GEODES

### Load packages
library(ggplot2)
library(cowplot)
library(reshape2)
#library(GeneCycle)

### Load data (start with only one to save RAM and comment the rest out)
# Normalized read tables
snorm <- read.csv("D:/geodes_data_tables/Sparkling_normalized.csv", header = T, row.names = 1)
#tnorm <- read.csv("D:/geodes_data_tables/Trout_normalized.csv", header = T, row.names = 1)
#mnorm <- read.csv("D:/geodes_data_tables/Mendota_normalized.csv", header = T, row.names = 1)

# Gene keys
#mendota_key <- read.csv("D:/geodes_data_tables/Mendota_ID90_genekey.csv", header = T)
spark_key <- read.csv("D:/geodes_data_tables/Sparkling_ID90_genekey.csv", header = T)
#trout_key <- read.csv("D:/geodes_data_tables/Trout_ID90_genekey.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Alex/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# How expressed is each phylum?
mendota_key$Taxonomy <- gsub("Bacteria;", "", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("Eukaryota;", "", mendota_key$Taxonomy)
mendota_key$Phylum <- sapply(strsplit(as.character(mendota_key$Taxonomy),";"), `[`, 1)

mnorm$Genes <- rownames(mnorm)
mnorm <- melt(mnorm)
mnorm$Timepoint <- metadata$Timepoint[match(mnorm$variable, metadata$Sample)]
mnorm$Taxonomy <- mendota_key$Phylum[match(mnorm$Genes, mendota_key$Gene)]
mnorm$Taxonomy <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("NO CLASSIFICATION MH", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("NO CLASSIFICATION LP", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("None", "Unclassified", mnorm$Taxonomy)
averaged_tax <- aggregate(value ~ Taxonomy + Timepoint, data = mnorm, mean)
wide_mnorm <- reshape(averaged_tax, idvar = "Taxonomy", timevar = "Timepoint", direction = "wide")
rownames(wide_mnorm) <- wide_mnorm$Taxonomy
wide_mnorm <- wide_mnorm[, 2:dim(wide_mnorm)[2]]
wide_mnorm <- wide_mnorm[which(rowSums(wide_mnorm) > 3000),]

mendota_phyla <- data.frame(Taxonomy = rownames(wide_mnorm), Sums = rowSums(wide_mnorm))
mendota_phyla$Taxonomy <- c("Unclassified", "Acidobacteria", "Actinobacteria", "Armatimonadetes", "Arthropoda", "Ascomycota", "Bacillariophyta", "Bacteroidetes", "Chlamydiae", "Chloroflexi", "Chlorophyta", "Chytridiomycota", "Cryptophyta", "Cyanobacteria", "Fibrobacteres", "Firmicutes", "Gemmatimonadetes", "Haptophyta", "Heterokonta", "Ignavibacteriae", "Planctomycetes", "Proteobacteria", "Spirochaetes", "Streptophyta", "Tenericutes", "Unclassified", "Unclassified", "Unclassified", "Oligohymenophorea", "Unclassified", "Isochyrysidales", "Perkinsida", "Verrucomicrobia", "Viruses")
unclassified <- mendota_phyla[which(mendota_phyla$Taxonomy == "Unclassified"), ]
mendota_phyla <- mendota_phyla[which(mendota_phyla$Taxonomy != "Unclassified"), ]
mendota_phyla <- rbind(mendota_phyla, c("Unclassified", sum(unclassified$Sums)))
mendota_phyla$Sums <- as.numeric(mendota_phyla$Sums)
mendota_phyla$Taxonomy <- factor(mendota_phyla$Taxonomy, levels = mendota_phyla$Taxonomy[order(mendota_phyla$Sums, decreasing = T)])
mendota_phyla$Type <- c("Bacteria", "Bacteria", "Bacteria", "Animals", "Fungi", "Algae", "Bacteria", "Bacteria", "Bacteria", "Algae", "Fungi", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Bacteria", "Protists", "Algae", "Protists", "Bacteria", "Viruses", "Unclassified")


p <- ggplot(mendota_phyla, aes(x = Taxonomy, y = Sums, fill = Type)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) + scale_fill_brewer(palette = "Set2") + labs(x = "", y = "Read Counts", title = "Lake Mendota Metatranscriptomes")

save_plot("C:/Users/Alex/Desktop/geodes/Plots/mendota_expression_by_phyla.pdf", p, base_height = 5, base_aspect_ratio = 1.5)

trout_key$Taxonomy <- gsub("Bacteria;", "", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("Eukaryota;", "", trout_key$Taxonomy)
trout_key$Phylum <- sapply(strsplit(as.character(trout_key$Taxonomy),";"), `[`, 1)

tnorm$Genes <- rownames(tnorm)
tnorm <- melt(tnorm)
tnorm$Timepoint <- metadata$Timepoint[match(tnorm$variable, metadata$Sample)]
tnorm$Taxonomy <- trout_key$Phylum[match(tnorm$Genes, trout_key$Gene)]
tnorm$Taxonomy <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("NO CLASSIFICATION MH", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("NO CLASSIFICATION LP", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("None", "Unclassified", tnorm$Taxonomy)
averaged_tax <- aggregate(value ~ Taxonomy + Timepoint, data = tnorm, mean)
wide_tnorm <- reshape(averaged_tax, idvar = "Taxonomy", timevar = "Timepoint", direction = "wide")
rownames(wide_tnorm) <- wide_tnorm$Taxonomy
wide_tnorm <- wide_tnorm[, 2:dim(wide_tnorm)[2]]
wide_tnorm <- wide_tnorm[which(rowSums(wide_tnorm) > 3000),]

trout_phyla <- data.frame(Taxonomy = rownames(wide_tnorm), Sums = rowSums(wide_tnorm))
trout_phyla$Taxonomy <- c("Acidobacteria", "Actinobacteria", "Armatimonadetes", "Arthropoda", "Bacillariophyta", "Bacteroidetes", "Chlorobi", "Chloroflexi", "Chlorophyta", "Cryptophyta", "Cyanobacteria", "Deinococcus-Thermus", "Firmicutes", "Haptophyta", "Heterokonta", "Ignavibacteriae", "Planctomycetes", "Porifera", "Proteobacteria", "Streptophyta", "Unclassified", "Unclassified", "Unclassified", "Isochyrysidales", "Verrucomicrobia", "Viruses")
unclassified <- trout_phyla[which(trout_phyla$Taxonomy == "Unclassified"), ]
trout_phyla <- trout_phyla[which(trout_phyla$Taxonomy != "Unclassified"), ]
trout_phyla <- rbind(trout_phyla, c("Unclassified", sum(unclassified$Sums)))
trout_phyla$Sums <- as.numeric(trout_phyla$Sums)
trout_phyla$Taxonomy <- factor(trout_phyla$Taxonomy, levels = trout_phyla$Taxonomy[order(trout_phyla$Sums, decreasing = T)])
trout_phyla$Type <- c("Bacteria", "Bacteria", "Bacteria", "Animals", "Fungi", "Bacteria", "Bacteria", "Bacteria", "Algae", "Fungi", "Bacteria", "Bacteria", "Bacteria", "Algae", "Algae", "Bacteria", "Bacteria", "Animals", "Bacteria", "Algae", "Algae", "Bacteria", "Viruses", "Unclassified")


p <- ggplot(trout_phyla, aes(x = Taxonomy, y = Sums, fill = Type)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) + scale_fill_brewer(palette = "Set2") + labs(x = "", y = "Read Counts", title = "Trout Bog Metatranscriptomes")

save_plot("C:/Users/Alex/Desktop/geodes/Plots/trout_expression_by_phyla.pdf", p, base_height = 5, base_aspect_ratio = 1.5)

spark_key$Taxonomy <- gsub("Bacteria;", "", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("Eukaryota;", "", spark_key$Taxonomy)
spark_key$Phylum <- sapply(strsplit(as.character(spark_key$Taxonomy),";"), `[`, 1)

snorm$Genes <- rownames(snorm)
snorm <- melt(snorm)
snorm$Timepoint <- metadata$Timepoint[match(snorm$variable, metadata$Sample)]
snorm$Taxonomy <- spark_key$Phylum[match(snorm$Genes, spark_key$Gene)]
snorm$Taxonomy <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("NO CLASSIFICATION MH", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("NO CLASSIFICATION LP", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("None", "Unclassified", snorm$Taxonomy)
averaged_tax <- aggregate(value ~ Taxonomy + Timepoint, data = snorm, mean)
wide_snorm <- reshape(averaged_tax, idvar = "Taxonomy", timevar = "Timepoint", direction = "wide")
rownames(wide_snorm) <- wide_snorm$Taxonomy
wide_snorm <- wide_snorm[, 2:dim(wide_snorm)[2]]
wide_snorm <- wide_snorm[which(rowSums(wide_snorm) > 3000),]

spark_phyla <- data.frame(Taxonomy = rownames(wide_snorm), Sums = rowSums(wide_snorm))
spark_phyla$Taxonomy <- c("Unclassified", "Actinobacteria", "Armatimonadetes", "Arthropoda", "Bacillariophyta", "Bacteroidetes", "Chlamydiae", "Chlorobi", "Chloroflexi", "Chlorophyta", "Chytridiomycota", "Cryptophyta", "Cyanobacteria", "Deinococcus-Thermus", "Firmicutes", "Gemmatimonadetes", "Haptophyta", "Heterokonta", "Nitrospirae", "Phaeophyceae", "Planctomycetes", "Proteobacteria", "Spirochaetes", "Streptophyta", "Unclassified", "Unclassified","Oligohymenophorea", "Unclassified", "Isochyrysidales","Unclassified", "Verrucomicrobia", "Viruses")
unclassified <- spark_phyla[which(spark_phyla$Taxonomy == "Unclassified"), ]
spark_phyla <- spark_phyla[which(spark_phyla$Taxonomy != "Unclassified"), ]
spark_phyla <- rbind(spark_phyla, c("Unclassified", sum(unclassified$Sums)))
spark_phyla$Sums <- as.numeric(spark_phyla$Sums)
spark_phyla$Taxonomy <- factor(spark_phyla$Taxonomy, levels = spark_phyla$Taxonomy[order(spark_phyla$Sums, decreasing = T)])
spark_phyla$Type <- c("Bacteria", "Bacteria", "Animals", "Fungi", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Fungi", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Algae", "Bacteria", "Algae", "Bacteria", "Bacteria", "Bacteria", "Algae", "Protists", "Algae", "Bacteria", "Viruses", "Unclassified")


p <- ggplot(spark_phyla, aes(x = Taxonomy, y = Sums, fill = Type)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) + scale_fill_brewer(palette = "Set2") + labs(x = "", y = "Read Counts", title = "Sparkling Lake Metatranscriptomes")

save_plot("C:/Users/Alex/Desktop/geodes/Plots/spark_expression_by_phyla.pdf", p, base_height = 5, base_aspect_ratio = 1.5)

# How does expression compare to abundance (ie the metagenomes?)
# No internal standard in the metagenomes, so normalized by library size (or report in % reads)
metaG_reads <- read.table("D:/geodes_data_tables/GEODES_metaG_2018-01-26.readcounts.txt", row.names = 1, sep = "\t")
colnames(metaG_reads) <- c("GEODES005", "GEODES006", "GEODES057", "GEODES058", "GEODES117", "GEODES118", "GEODES165", "GEODES166", "GEODES167", "GEODES168")
metaG_key <- read.table("D:/geodes_data_tables/GEODES_metaG_genekey.txt", sep = "\t", quote = "")
colnames(metaG_key) <- c("Gene", "Genome", "Taxonomy", "Product")
lakekey <- c("Sparkling", "Sparkling", "Trout", "Trout", "Mendota", "Mendota", "Sparkling2009", "Sparkling2009", "Sparkling2009", "Sparkling2009")
metaG_reads <- sweep(metaG_reads, 2, colSums(metaG_reads), "/")

metaG_key$Taxonomy <- gsub("Bacteria;", "", metaG_key$Taxonomy)
metaG_key$Taxonomy <- gsub("Eukaryota;", "", metaG_key$Taxonomy)
metaG_key$Phylum <- sapply(strsplit(as.character(metaG_key$Taxonomy),";"), `[`, 1)

metaG_reads$Genes <- rownames(metaG_reads)
spark_metaG <- metaG_reads[,c(1,2, 11)]
trout_metaG <- metaG_reads[,c(3,4, 11)]
mendota_metaG <- metaG_reads[,c(5,6, 11)]
spark2_metaG <- metaG_reads[,c(7,8,9,10, 11)]

spark_metaG <- melt(spark_metaG)
trout_metaG <- melt(trout_metaG)
mendota_metaG <- melt(mendota_metaG)
spark2_metaG <- melt(spark2_metaG)

spark_metaG$Phylum <- metaG_key$Phylum[match(spark_metaG$Genes, metaG_key$Gene)]
spark_metaG$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified unclassified unclassified", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("None", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified unclassified", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("unclassified", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Unclassified ", "Unclassified", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", spark_metaG$Phylum)
spark_metaG$Phylum[which(is.na(spark_metaG$Phylum))] <- "Unclassified"
spark_phyla <- aggregate(value ~ Phylum, data = spark_metaG, mean)

spark_phyla$Type <- c("Unclassified", "Bacteria", "Bacteria", "Bacteria", "Animals", "Fungi", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Fungi", "Protists", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Algae", "Bacteria", "Bacteria",  "Bacteria", "Bacteria", "Algae", "Protists", "Algae", "Bacteria", "Animals", "Bacteria", "Bacteria", "Algae", "Bacteria", "Unclassified", "Bacteria", "Viruses")

spark_phyla$Phylum <- factor(spark_phyla$Phylum, levels = spark_phyla$Phylum[order(spark_phyla$value, decreasing = T)])

p <- ggplot(spark_phyla, aes(x = Phylum, y = value, fill = Type)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) + scale_fill_brewer(palette = "Set2") + labs(x = "", y = "Read Counts", title = "Sparkling Lake Metagenomes")

save_plot("C:/Users/Alex/Desktop/geodes/Plots/spark_metagenome_by_phyla.pdf", p, base_height = 5, base_aspect_ratio = 1.6)

trout_metaG$Phylum <- metaG_key$Phylum[match(trout_metaG$Genes, metaG_key$Gene)]
trout_metaG$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified unclassified unclassified", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("None", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified unclassified", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("unclassified", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Unclassified ", "Unclassified", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", trout_metaG$Phylum)
trout_metaG$Phylum[which(is.na(trout_metaG$Phylum) == T)] <- "Unclassified"
trout_phyla <- aggregate(value ~ Phylum, data = trout_metaG, mean)

trout_phyla$Type <- c("Unclassified", "Bacteria", "Bacteria", "Bacteria", "Animals", "Fungi", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Fungi", "Protists", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Algae", "Bacteria", "Bacteria",  "Bacteria", "Bacteria", "Algae", "Protists", "Algae", "Bacteria", "Animals", "Bacteria", "Bacteria", "Algae", "Bacteria", "Unclassified", "Bacteria", "Viruses")

trout_phyla$Phylum <- factor(trout_phyla$Phylum, levels = trout_phyla$Phylum[order(trout_phyla$value, decreasing = T)])

p <- ggplot(trout_phyla, aes(x = Phylum, y = value, fill = Type)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) + scale_fill_brewer(palette = "Set2") + labs(x = "", y = "Read Counts", title = "Trout Bog Metagenomes")

save_plot("C:/Users/Alex/Desktop/geodes/Plots/trout_metagenome_by_phyla.pdf", p, base_height = 5, base_aspect_ratio = 1.6)

mendota_metaG$Phylum <- metaG_key$Phylum[match(mendota_metaG$Genes, metaG_key$Gene)]
mendota_metaG$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified unclassified unclassified", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("None", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified unclassified", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("unclassified", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Unclassified ", "Unclassified", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", mendota_metaG$Phylum)
mendota_metaG$Phylum[which(is.na(mendota_metaG$Phylum) == T)] <- "Unclassified"
mendota_phyla <- aggregate(value ~ Phylum, data = mendota_metaG, mean)

mendota_phyla$Type <- c("Unclassified", "Bacteria", "Bacteria", "Bacteria", "Animals", "Fungi", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Fungi", "Protists", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Algae", "Bacteria", "Bacteria",  "Bacteria", "Bacteria", "Algae", "Protists", "Algae", "Bacteria", "Animals", "Bacteria", "Bacteria", "Algae", "Bacteria", "Unclassified", "Bacteria", "Viruses")

mendota_phyla$Phylum <- factor(mendota_phyla$Phylum, levels = mendota_phyla$Phylum[order(mendota_phyla$value, decreasing = T)])

p <- ggplot(mendota_phyla, aes(x = Phylum, y = value, fill = Type)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) + scale_fill_brewer(palette = "Set2") + labs(x = "", y = "Read Counts", title = "Mendota Metagenomes")

save_plot("C:/Users/Alex/Desktop/geodes/Plots/mendota_metagenome_by_phyla.pdf", p, base_height = 5, base_aspect_ratio = 1.6)

spark2_metaG$Phylum <- metaG_key$Phylum[match(spark2_metaG$Genes, metaG_key$Gene)]
spark2_metaG$Phylum <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("unclassified unclassified unclassified unclassified", "Unclassified", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("unclassified unclassified unclassified", "Unclassified", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("NO CLASSIFICATION MH", "Unclassified", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("NO CLASSIFICATION LP", "Unclassified", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("None", "Unclassified", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("unclassified unclassified", "Unclassified", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("unclassified Oligohymenophorea", "Ciliophora", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("unclassified Pelagophyceae", "Ochrophyta", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("unclassified", "Unclassified", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("Unclassified ", "Unclassified", spark2_metaG$Phylum)
spark2_metaG$Phylum <- gsub("UnclassifiedIsochrysidales", "Haptophyta", spark2_metaG$Phylum)
spark2_metaG$Phylum[which(is.na(spark2_metaG$Phylum) == T)] <- "Unclassified"
spark2_phyla <- aggregate(value ~ Phylum, data = spark2_metaG, mean)

spark2_phyla$Type <- c("Unclassified", "Bacteria", "Bacteria", "Bacteria", "Animals", "Fungi", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Fungi", "Protists", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Algae", "Algae", "Bacteria", "Bacteria",  "Bacteria", "Bacteria", "Algae", "Protists", "Algae", "Bacteria", "Animals", "Bacteria", "Bacteria", "Algae", "Bacteria", "Unclassified", "Bacteria", "Viruses")

spark2_phyla$Phylum <- factor(spark2_phyla$Phylum, levels = spark2_phyla$Phylum[order(spark2_phyla$value, decreasing = T)])

p <- ggplot(spark2_phyla, aes(x = Phylum, y = value, fill = Type)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1)) + scale_fill_brewer(palette = "Set2") + labs(x = "", y = "Read Counts", title = "Sparkling Lake 2009 Metagenomes")

save_plot("C:/Users/Alex/Desktop/geodes/Plots/spark09_metagenome_by_phyla.pdf", p, base_height = 5, base_aspect_ratio = 1.6)

# What are the top expressed genes?

