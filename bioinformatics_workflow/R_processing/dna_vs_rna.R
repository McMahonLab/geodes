# Make plots of DNA vs RNA by lake by phylum
library(ggplot2)
library(cowplot)
library(reshape2)
library(ggrepel)
theme_set(theme_cowplot(font_size=18))
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# Metagenomic data
metaG_reads <- read.table("E:/geodes_data_tables/GEODES_metaG_ID90_2018-03-10.readcounts.txt", row.names = 1, sep = "\t")
colnames(metaG_reads) <- c("GEODES005", "GEODES006", "GEODES057", "GEODES058", "GEODES117", "GEODES118", "GEODES165", "GEODES166", "GEODES167", "GEODES168")
metaG_key <- read.table("E:/geodes_data_tables/GEODES_metaG_genekey_2018-03-12.txt", sep = "\t", quote = "")
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


spark_metaG <- melt(spark_metaG)
trout_metaG <- melt(trout_metaG)
mendota_metaG <- melt(mendota_metaG)

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
mendota_metaG$Phylum <- gsub("TM7", "Candidatus Saccharibacteria", mendota_metaG$Phylum)
mendota_metaG$Phylum <- gsub("Ignavibacteriae", "Ignavibacteria", mendota_metaG$Phylum)
mendota_metaG$Phylum[which(is.na(mendota_metaG$Phylum) == T)] <- "Unclassified"

mendota_metaG_phyla <- aggregate(value ~ Phylum, data = mendota_metaG, mean)

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
spark_metaG$Phylum <- gsub("TM7", "Candidatus Saccharibacteria", spark_metaG$Phylum)
spark_metaG$Phylum <- gsub("Ignavibacteriae", "Ignavibacteria", spark_metaG$Phylum)
spark_metaG$Phylum[which(is.na(spark_metaG$Phylum))] <- "Unclassified"
spark_metaG_phyla <- aggregate(value ~ Phylum, data = spark_metaG, mean)

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
trout_metaG$Phylum <- gsub("TM7", "Candidatus Saccharibacteria", trout_metaG$Phylum)
trout_metaG$Phylum <- gsub("Ignavibacteriae", "Ignavibacteria", trout_metaG$Phylum)
trout_metaG$Phylum[which(is.na(trout_metaG$Phylum))] <- "Unclassified"
trout_metaG_phyla <- aggregate(value ~ Phylum, data = trout_metaG, mean)


# Remove unused metaG datasets to sve RAM
rm(metaG_reads)
rm(mendota_metaG)
rm(trout_metaG)
rm(spark_metaG)
rm(metaG_key)
#Mendota
mnorm <- read.csv("E:/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mendota_key <- read.csv("E:/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)

# How expressed is each phylum?
mendota_key$Taxonomy <- gsub("Bacteria;", "", mendota_key$Taxonomy)
mendota_key$Taxonomy <- gsub("Eukaryota;", "", mendota_key$Taxonomy)
mendota_key$Phylum <- sapply(strsplit(as.character(mendota_key$Taxonomy),";"), `[`, 1)

mnorm$Genes <- rownames(mnorm)
mnorm <- melt(mnorm)
mnorm$variable <- gsub(".nonrRNA", "", mnorm$variable)
mnorm$Timepoint <- metadata$Timepoint[match(mnorm$variable, metadata$Sample)]
mnorm$Taxonomy <- mendota_key$Phylum[match(mnorm$Genes, mendota_key$Gene)]
mnorm$Taxonomy <- gsub("Cryptophyta,Cryptophyceae,Pyrenomonadales,Geminigeraceae,Guillardia,theta", "Cryptophyta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Haptophyta,Prymnesiophyceae,Isochrysidales,Noelaerhabdaceae,Emiliania,huxleyi", "Haptophyta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Heterokonta,Pelagophyceae,Pelagomonadales,Pelagomonadaceae,Aureococcus,anophagefferens", "Heterokonta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Heterokonta,Bacillariophyceae,Naviculales,Phaeodactylaceae,Phaeodactylum,tricornutum", "Heterokonta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Heterokonta,Ochrophyta,Eustigmataphyceae,Eustigmataceae,Nannochloropsis,gaditana", "Heterokonta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified unclassified unclassified unclassified unclassified", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("NO CLASSIFICATION MH", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("NO CLASSIFICATION LP", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("None", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified unclassified", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified Oligohymenophorea", "Ciliophora", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified Pelagophyceae", "Ochrophyta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("unclassified", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("Unclassified ", "Unclassified", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("UnclassifiedIsochrysidales", "Haptophyta", mnorm$Taxonomy)
mnorm$Taxonomy <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", mnorm$Taxonomy)
#Remove unclassifieds to save space
mnorm <- mnorm[which(mnorm$Taxonomy != "Unclassified" & mnorm$Taxonomy != "internal standard"),]
averaged_tax <- aggregate(value ~ Taxonomy + Timepoint, data = mnorm, mean)
wide_mnorm <- reshape(averaged_tax, idvar = "Taxonomy", timevar = "Timepoint", direction = "wide")
rownames(wide_mnorm) <- wide_mnorm$Taxonomy
wide_mnorm <- wide_mnorm[, 2:dim(wide_mnorm)[2]]
wide_mnorm <- wide_mnorm[which(rowSums(wide_mnorm) > 3000),]

mendota_phyla <- data.frame(Taxonomy = rownames(wide_mnorm), Sums = rowSums(wide_mnorm))

rm(mendota_key)
rm(mnorm)

# Choose the intersection of the top 10 expressed and abundant
keep <- c("Chloroflexi", "Actinobacteria", "Ignavibacteria", "Planctomycetes", "Crenarchaeaota", "Cyanobacteria", "Gemmatimonadetes", "Bacteroidetes", "Proteobacteria", "Fibrobacteres", "Cryptophyta", "Heterokonta", "Viruses")
mendota_intersect <- mendota_phyla[which(mendota_phyla$Taxonomy %in% keep), ]
mendota_intersect$metaG <- mendota_metaG_phyla$value[match(mendota_intersect$Taxonomy, mendota_metaG_phyla$Phylum)]
mendota_intersect$Taxonomy <- factor(mendota_intersect$Taxonomy, levels = mendota_intersect$Taxonomy[order(mendota_intersect$metaG, decreasing = F)])
mendota_intersect <- mendota_intersect[order(mendota_intersect$metaG, decreasing = F),]

mendota_intersect$Type <- c("Virus", "Algae", "Algae", "Bacteria", "Bacteria", "Bacteria", "Bacteria", "Archaea", "Bacteria", "Bacteria", "Bacteria", "Bacteria")

p <- ggplot(mendota_intersect[which(mendota_intersect$Taxonomy != "Chloroflexi"), ], aes(x = metaG, y = Sums, color = Type)) + geom_point(size = 2.5) + geom_label_repel(aes(label = Taxonomy, color = Type), force = 5, size = 7.5) + scale_color_manual(values = c("limegreen", "lavenderblush4", "royalblue", "goldenrod")) + labs(x = "Proportion of metagenomic reads assigned", y = "Transcripts/L assigned", title = "Lake Mendota") + theme(legend.position = "none")

save_plot("C:/Users/Goose and Gander/Desktop/geodes/Plots/mendota_dna_vs_rna_no_chloroflexi.pdf", p, base_height = 6, base_aspect_ratio = 8/6)

#Clear mendota data from workspace and start with sparkling

rm(averaged_tax)
rm(mendota_intersect)
rm(wide_mnorm)
rm(mendota_metaG_phyla)
rm(mendota_phyla)

snorm <- read.csv("E:/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
spark_key <- read.csv("E:/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)

spark_key$Taxonomy <- gsub("Bacteria;", "", spark_key$Taxonomy)
spark_key$Taxonomy <- gsub("Eukaryota;", "", spark_key$Taxonomy)
spark_key$Phylum <- sapply(strsplit(as.character(spark_key$Taxonomy),";"), `[`, 1)

snorm$Genes <- rownames(snorm)
snorm <- melt(snorm)
snorm$variable <- gsub(".nonrRNA", "", snorm$variable)
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
snorm$Taxonomy <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified unclassified", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified Oligohymenophorea", "Ciliophora", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified Pelagophyceae", "Ochrophyta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("unclassified", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Unclassified ", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("UnclassifiedIsochrysidales", "Haptophyta", snorm$Taxonomy)
snorm$Taxonomy <- gsub("UnclassifiedUnclassified", "Unclassified", snorm$Taxonomy)
snorm$Taxonomy <- gsub("Candidatus Saccharibacteria", "TM7", snorm$Taxonomy)
snorm$Taxonomy <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", snorm$Taxonomy)
#Remove unclassifieds to save space
snorm <- snorm[which(snorm$Taxonomy != "Unclassified" & snorm$Taxonomy != "internal standard"),]
averaged_tax <- aggregate(value ~ Taxonomy + Timepoint, data = snorm, mean)
wide_snorm <- reshape(averaged_tax, idvar = "Taxonomy", timevar = "Timepoint", direction = "wide")
rownames(wide_snorm) <- wide_snorm$Taxonomy
wide_snorm <- wide_snorm[, 2:dim(wide_snorm)[2]]
wide_snorm <- wide_snorm[which(rowSums(wide_snorm) > 3000),]

spark_phyla <- data.frame(Taxonomy = rownames(wide_snorm), Sums = rowSums(wide_snorm))

rm(spark_key)
rm(snorm)

# Choose the intersection of the top 10 expressed and abundant
keep <- c("Firmicutes", "Heterokonta", "Armatimonadetes", "Haptophtya", "Chlamydiae", "Verrucomicrobia", "Actinobacteria", "Viruses", "Planctomycetes", "Phaeophyceae", "Cyanobacteria", "TM7", "Bacteroidetes", "Proteobacteria", "Cryptophyta")
sparkling_intersect <- spark_phyla[which(spark_phyla$Taxonomy %in% keep), ]
sparkling_intersect$metaG <- spark_metaG_phyla$value[match(sparkling_intersect$Taxonomy, spark_metaG_phyla$Phylum)]
sparkling_intersect$Taxonomy <- factor(sparkling_intersect$Taxonomy, levels = sparkling_intersect$Taxonomy[order(sparkling_intersect$metaG, decreasing = F)])
sparkling_intersect <- sparkling_intersect[order(sparkling_intersect$metaG, decreasing = F),]

sparkling_intersect$Type <- c("Bacteria", "Bacteria", "Bacteria", "Algae", "Algae", "Bacteria", "Viruses", "Bacteria", "Bacteria", "Bacteria", "Algae", "Bacteria", "Bacteria")
sparkling_intersect$metaG[which(is.na(sparkling_intersect$metaG) == T)] <- 0

p <- ggplot(sparkling_intersect, aes(x = metaG, y = Sums, color = Type)) + geom_point(size = 2.5) + geom_label_repel(aes(label = Taxonomy, color = Type), force = 5, size = 7.5) + scale_color_manual(values = c("limegreen", "royalblue", "goldenrod")) + labs(x = "Proportion of metagenomic reads assigned", y = "Transcripts/L assigned", title = "Sparkling Lake") + theme(legend.position = "none")
#export as 6x8in pdf

save_plot("C:/Users/Goose and Gander/Desktop/geodes/Plots/spark_dna_vs_rna.pdf", p, base_height = 6, base_aspect_ratio = 8/6)

rm(averaged_tax)
rm(sparkling_intersect)
rm(wide_snorm)
rm(spark_metaG_phyla)
rm(spark_phyla)


#Finally Trout bog

tnorm <- read.csv("E:/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
trout_key <- read.csv("E:/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

trout_key$Taxonomy <- gsub("Bacteria;", "", trout_key$Taxonomy)
trout_key$Taxonomy <- gsub("Eukaryota;", "", trout_key$Taxonomy)
trout_key$Phylum <- sapply(strsplit(as.character(trout_key$Taxonomy),";"), `[`, 1)

tnorm$Genes <- rownames(tnorm)
tnorm <- melt(tnorm)
tnorm$variable <- gsub(".nonrRNA", "", tnorm$variable)
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
tnorm$Taxonomy <- gsub("NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified unclassified Perkinsida", "Perkinsozoa", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified unclassified", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified Oligohymenophorea", "Ciliophora", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified Pelagophyceae", "Ochrophyta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("unclassified", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Unclassified ", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("UnclassifiedIsochrysidales", "Haptophyta", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("UnclassifiedUnclassified", "Unclassified", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("Candidatus Saccharibacteria", "TM7", tnorm$Taxonomy)
tnorm$Taxonomy <- gsub("NO CLASSIFICATION BASED ON GIVEN PHYLODIST", "Unclassified", tnorm$Taxonomy)
#Remove unclassifieds to save space
tnorm <- tnorm[which(tnorm$Taxonomy != "Unclassified" & tnorm$Taxonomy != "internal standard"),]
averaged_tax <- aggregate(value ~ Taxonomy + Timepoint, data = tnorm, mean)
wide_tnorm <- reshape(averaged_tax, idvar = "Taxonomy", timevar = "Timepoint", direction = "wide")
rownames(wide_tnorm) <- wide_tnorm$Taxonomy
wide_tnorm <- wide_tnorm[, 2:dim(wide_tnorm)[2]]
wide_tnorm <- wide_tnorm[which(rowSums(wide_tnorm) > 3000),]

trout_phyla <- data.frame(Taxonomy = rownames(wide_tnorm), Sums = rowSums(wide_tnorm))

rm(tnorm)
rm(trout_key)

# Choose the intersection of the top 10 expressed and abundant
keep <- c("Verrucomicrobia", "Armatimonadetes", "Viruses", "Chlorophyta", "Acidobacteria", "Heterokonta", "Streptophyta", "Actinobacteria", "Bacteroidetes", "Cryptophyta", "Cyanobacteria", "Deinococcus-Thermus")
trout_intersect <- trout_phyla[which(trout_phyla$Taxonomy %in% keep), ]
trout_intersect$metaG <- trout_metaG_phyla$value[match(trout_intersect$Taxonomy, trout_metaG_phyla$Phylum)]
trout_intersect$Taxonomy <- factor(trout_intersect$Taxonomy, levels = trout_intersect$Taxonomy[order(trout_intersect$metaG, decreasing = F)])
trout_intersect <- trout_intersect[order(trout_intersect$metaG, decreasing = F),]

trout_intersect$Type <- c("Bacteria", "Bacteria", "Algae", "Bacteria", "Bacteria", "Algae", "Algae", "Bacteria", "Algae", "Viruses", "Bacteria", "Bacteria")

p <- ggplot(trout_intersect, aes(x = metaG, y = Sums, color = Type)) + geom_point(size = 2.5) + geom_label_repel(aes(label = Taxonomy, color = Type), force = 15, size = 7.5) + scale_color_manual(values = c("limegreen", "royalblue", "goldenrod")) + labs(x = "Proportion of metagenomic reads assigned", y = "Transcripts/L assigned", title = "Trout Bog") + theme(legend.position = "none")
#export as 6x8in pdf

save_plot("C:/Users/Goose and Gander/Desktop/geodes/Plots/trout_dna_vs_rna.pdf", p, base_height = 6, base_aspect_ratio = 8/6)

