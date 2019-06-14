# Set up environment
path <- "/Users/Alex/Desktop/"
path2 <- "/Users/Alex/"

library(ggplot2)
library(cowplot)
library(reshape2)
library(DESeq2)

zscore <- function(counts){
  z <- (counts - mean(counts)) / sd(counts)
  return(z)
}

# Sample data
metadata <- read.csv(file = paste(path2, "Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", sep = ""), header = T)

mnorm <- read.csv(paste(path, "geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)
mnorm <- mnorm[, which(colnames(mnorm) != "GEODES158.nonrRNA")]
colnames(mnorm) <- gsub(".nonrRNA", "", colnames(mnorm))
mendota_key <- read.csv(paste(path, "geodes_data_tables/Mendota_ID90_genekey_geneclassifications_2018-11-28.csv", sep = ""), header = T)

tnorm <- read.csv(paste(path, "geodes_data_tables/Trout_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)
tnorm <- tnorm[, which(colnames(tnorm) != "GEODES065.nonrRNA")]
colnames(tnorm) <- gsub(".nonrRNA", "", colnames(tnorm))
trout_key <- read.csv(paste(path, "geodes_data_tables/Trout_ID90_genekey_geneclassifications_2018-11-28.csv", sep = ""), header = T)

snorm <- read.csv(paste(path, "geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)
snorm <- snorm[, which(colnames(snorm) != "GEODES014.nonrRNA" & colnames(snorm) != "GEODES033.nonrRNA")]
colnames(snorm) <- gsub(".nonrRNA", "", colnames(snorm))
spark_key <- read.csv(paste(path, "geodes_data_tables/Sparkling_ID90_genekey_geneclassifications_2018-11-28.csv", sep = ""), header = T)

abun_mnorm <- mnorm[order(rowSums(mnorm), decreasing = T), ]
abun_mnorm <- abun_mnorm[1:20000,]

colnames(abun_mnorm) <- gsub(".nonrRNA", "", colnames(abun_mnorm))

input <- as.matrix(abun_mnorm)
input <- input/100
input <- round(input, digits = 0)

conditions <- metadata$Time[match(colnames(abun_mnorm), metadata$Sample)]
conditions[which(conditions == 9 | conditions == 13 | conditions == 17)] <- "day"
conditions[which(conditions == 5 | conditions == 21 | conditions == 1)] <- "night"
coldata <- data.frame(samples = colnames(abun_mnorm), conditions)

cds <- DESeqDataSetFromMatrix(countData = input,
                              colData = coldata,
                              design = ~ conditions)

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds) 


dds <- DESeq(cds)
res <- results(dds)

#res <- res[which(res$padj < 0.05), ]
reskey <- mendota_key[match(rownames(res), mendota_key$Gene), ]
reskey$log2FoldChange <- res$log2FoldChange[match(reskey$Gene, rownames(res))]

# Make column by certain key words 
reskey$Category <- "None"
reskey$Category[grep("photo|Photo", reskey$Product)] <- "Photosynthesis"
reskey$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", reskey$Product)] <- "Rhodopsin"
#reskey$Category[grep("rbcL|ribulose-bisphosphate carboxylase", reskey$Product)] <- "RuBisCO"
#reskey$Category[grep("citrate lyase|Citrate lyase", reskey$Product)] <- "rTCA"
#reskey$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", reskey$Product)] <- "Nitrogenase"
reskey$Category[grep("peroxidase|peroxide|catalase", reskey$Product)] <- "ROS"
#reskey$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", reskey$Product)] <- "Nitrite_reduction"
#reskey$Category[grep("urease", reskey$Product)] <- "Urease"
reskey$Category[grep("protease", reskey$Product)] <- "Protease"
reskey$Category[grep("sugar transport|carbohydrate ABC transport|Carbohydrate-selective porin", reskey$Product)] <- "Sugar transport"
reskey$Category[grep("raffinose/stachyose/melibiose transport", reskey$Product)] <- "R/S/M transport"
#reskey$Category[grep("chitobiose transport", reskey$Product)] <- "Chitobiose transport"
#reskey$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", reskey$Product)] <- "Chitinase"
reskey$Category[grep("ribose transport", reskey$Product)] <- "Ribose transport"

# Add day/night counts
gene_counts <- abun_mnorm[match(reskey$Gene, rownames(abun_mnorm)),]
reskey1 <- reskey
reskey2 <- reskey
reskey1$Counts <- rowSums(gene_counts[which(conditions == "day")])
reskey2$Counts <- rowSums(gene_counts[which(conditions == "night")])
reskey <- rbind(reskey1, reskey2)
reskey$Time <- c(rep("Day", dim(reskey1)[1]), rep("Night", dim(reskey2)[1]))
reskey$Counts <- as.numeric(reskey$Counts)

# Get totals:
formatC(sum(reskey$Counts[which(reskey$Category == "Sugar transport" & reskey$Time == "Day")]), format = "e", digits = 2)
formatC(sum(reskey$Counts[which(reskey$Category == "Sugar transport" & reskey$Time == "Night")]), format = "e", digits = 2)

# # Optional: use this code to get finer resolution taxonomy info
# library(dplyr)
# by_reskey <- reskey %>% group_by(Taxonomy, Category, Time)
# fine_taxa <- function(type, time){
#   x <- summarize(by_reskey, sum = sum(Counts))
#   y <- x[which(x$Category == type & x$Time == time), ]
#   return(y[order(y$sum, decreasing = T), ])
# }
# # example: fine_taxa("Sugar transport", "Day")

# Visually inspect phyla to keep
# ggplot(data = reskey, aes(x = Phylum, y = Counts, fill = Time)) + geom_bar(stat = "identity")  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0, 0)) + labs(x = NULL, y = "Proportion of reads", title = "Mendota") + facet_wrap(~ Category, ncol = 3, scales="free")

phyla2keep <- c("Actinobacteria", "Bacteroidetes", "Unclassified", "Cyanobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia", "Firmicutes", "Deltaproteobacteria", "Planctomycetes")
reskey <- reskey[which(reskey$Category != "None"), ] 
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
reskey$Phylum[which(is.na(reskey$Phylum) == T)] <- "Unclassified"
reskey$Phylum <- as.character(reskey$Phylum)
reskey$Phylum[which(reskey$Phylum %not in% phyla2keep)] <- "Other"

colorkey <- data.frame(phylum = c("Actinobacteria", "Alphaproteobacteria", "Armatimonadetes", "Bacteroidetes", "Betaproteobacteria", "Chlorophyta", "Cyanobacteria", "Deltaproteobacteria", "Eukaryota", "Firmicutes", "Gammaproteobacteria", "Other", "Planctomycetes", "Unclassified", "Verrucomicrobia"), color = c("chocolate1",  "hotpink2", "turquoise3", "goldenrod1", "dodgerblue", "palegreen4", "green3", "skyblue", "darkolivegreen3", "orchid3", "palegreen", "honeydew4", "lightcoral", "honeydew3", "lightslateblue"))

p1 <- ggplot(data = reskey, aes(x = Time, y = Counts, fill = Phylum)) + geom_bar(stat = "identity", position = "fill")  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0, 0)) + labs(x = NULL, y = "Proportion of reads", title = "A. Lake Mendota") + facet_wrap(~ Category, ncol = 3) + scale_fill_manual(values = as.character(colorkey$color[which(colorkey$phylum %in% unique(reskey$Phylum))])) + theme(legend.position = "none")
#save the mendota reskey to make the legend
MEreskey <- reskey

abun_tnorm <- tnorm[order(rowSums(tnorm), decreasing = T), ]
abun_tnorm <- abun_tnorm[1:20000,]

colnames(abun_tnorm) <- gsub(".nonrRNA", "", colnames(abun_tnorm))

input <- as.matrix(abun_tnorm)
input <- input/500
input <- round(input, digits = 0)

conditions <- metadata$Time[match(colnames(abun_tnorm), metadata$Sample)]
conditions[which(conditions == 9 | conditions == 13 | conditions == 17)] <- "day"
conditions[which(conditions == 5 | conditions == 21 | conditions == 1)] <- "night"
coldata <- data.frame(samples = colnames(abun_tnorm), conditions)

cds <- DESeqDataSetFromMatrix(countData = input,
                              colData = coldata,
                              design = ~ conditions)

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds) 


dds <- DESeq(cds)
res <- results(dds)
#res <- res[which(res$padj < 0.05), ]
reskey <- trout_key[match(rownames(res), trout_key$Gene), ]
reskey$log2FoldChange <- res$log2FoldChange[match(reskey$Gene, rownames(res))]

# Make column by certain key words 
reskey$Category <- "None"
reskey$Category[grep("photo|Photo", reskey$Product)] <- "Photosynthesis"
#reskey$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", reskey$Product)] <- "Rhodopsin"
reskey$Category[grep("rbcL|ribulose-bisphosphate carboxylase", reskey$Product)] <- "RuBisCO"
#reskey$Category[grep("citrate lyase|Citrate lyase", reskey$Product)] <- "rTCA"
#reskey$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", reskey$Product)] <- "Nitrogenase"
reskey$Category[grep("peroxidase|peroxide|catalase", reskey$Product)] <- "ROS"
#reskey$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", reskey$Product)] <- "Nitrite_reduction"
#reskey$Category[grep("urease", reskey$Product)] <- "Urease"
#reskey$Category[grep("protease", reskey$Product)] <- "Protease"
reskey$Category[grep("sugar transport|carbohydrate ABC transport|Carbohydrate-selective porin", reskey$Product)] <- "Sugar transport"
#reskey$Category[grep("R/S/M transport", reskey$Product)] <- "R/S/M transport"
#reskey$Category[grep("chitobiose transport", reskey$Product)] <- "Chitobiose transport"
#reskey$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", reskey$Product)] <- "Chitinase"
reskey$Category[grep("ribose transport", reskey$Product)] <- "Ribose transport"
#reskey$Category[grep("glucose/mannose transport", reskey$Product)] <- "Glucose/mannose transport"
#reskey$Category[grep("rhamnose transport", reskey$Product)] <- "Rhamnose transport"
reskey$Category[grep("xylose transport", reskey$Product)] <- "Xylose transport"
#reskey$Category[grep("fructose transport", reskey$Product)] <- "Fructose transport"
reskey$Category[grep("amino acid transport", reskey$Product)] <- "Amino acid transport"

# Add day/night counts
gene_counts <- abun_tnorm[match(reskey$Gene, rownames(abun_tnorm)),]
reskey1 <- reskey
reskey2 <- reskey
reskey1$Counts <- rowSums(gene_counts[which(conditions == "day")])
reskey2$Counts <- rowSums(gene_counts[which(conditions == "night")])
reskey <- rbind(reskey1, reskey2)
reskey$Time <- c(rep("Day", dim(reskey1)[1]), rep("Night", dim(reskey2)[1]))
reskey$Counts <- as.numeric(reskey$Counts)

formatC(sum(reskey$Counts[which(reskey$Category == "Xylose transport" & reskey$Time == "Day")]), format = "e", digits = 2)
formatC(sum(reskey$Counts[which(reskey$Category == "Xylose transport" & reskey$Time == "Night")]), format = "e", digits = 2)


# Optional: use this code to get finer resolution taxonomy info
# library(dplyr)
# by_reskey <- reskey %>% group_by(Taxonomy, Category, Time)
# fine_taxa <- function(type, time){
#   x <- summarize(by_reskey, sum = sum(Counts))
#   y <- x[which(x$Category == type & x$Time == time), ]
#   return(y[order(y$sum, decreasing = T), ])
# }
# example: fine_taxa("Sugar transport", "Day")

# Visually inspect phyla to keep
# ggplot(data = reskey, aes(x = Phylum, y = Counts, fill = Time)) + geom_bar(stat = "identity")  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0, 0)) + labs(x = NULL, y = "Proportion of reads", title = "Mendota") + facet_wrap(~ Category, ncol = 3, scales="free")

phyla2keep <- c("Actinobacteria", "Alphaproteobacteria", "Armatimonadetes", "Bacteroidetes", "Unclassified", "Cyanobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia", "Deltaproteobacteria", "Eukaryota", "Chlorophyta")
reskey <- reskey[which(reskey$Category != "None"), ] 
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
reskey$Phylum[which(is.na(reskey$Phylum) == T)] <- "Unclassified"
reskey$Phylum <- as.character(reskey$Phylum)
reskey$Phylum[which(reskey$Phylum %not in% phyla2keep)] <- "Other"

p2 <- ggplot(data = reskey, aes(x = Time, y = Counts, fill = Phylum)) + geom_bar(stat = "identity", position = "fill")  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0, 0)) + labs(x = NULL, y = "Proportion of reads", title = "B. Trout Bog") + facet_wrap(~ Category, ncol = 3) + scale_fill_manual(values = as.character(colorkey$color[which(colorkey$phylum %in% unique(reskey$Phylum))])) + theme(legend.position = "none")
TBreskey <- reskey

abun_snorm <- snorm[order(rowSums(snorm), decreasing = T), ]
abun_snorm <- abun_snorm[1:20000,]

colnames(abun_snorm) <- gsub(".nonrRNA", "", colnames(abun_snorm))

input <- as.matrix(abun_snorm)
input <- input/100
input <- round(input, digits = 0)

conditions <- metadata$Time[match(colnames(abun_snorm), metadata$Sample)]
conditions[which(conditions == 9 | conditions == 13 | conditions == 17)] <- "day"
conditions[which(conditions == 5 | conditions == 21 | conditions == 1)] <- "night"
coldata <- data.frame(samples = colnames(abun_snorm), conditions)

cds <- DESeqDataSetFromMatrix(countData = input,
                              colData = coldata,
                              design = ~ conditions)

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds) 


dds <- DESeq(cds)
res <- results(dds)

#res <- res[which(res$padj < 0.05), ]
reskey <- spark_key[match(rownames(res), spark_key$Gene), ]
reskey$log2FoldChange <- res$log2FoldChange[match(reskey$Gene, rownames(res))]

# Make column by certain key words 
reskey$Category <- "None"
reskey$Category[grep("photo|Photo", reskey$Product)] <- "Photosynthesis"
#reskey$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", reskey$Product)] <- "Rhodopsin"
#reskey$Category[grep("rbcL|ribulose-bisphosphate carboxylase", reskey$Product)] <- "RuBisCO"
#reskey$Category[grep("citrate lyase|Citrate lyase", reskey$Product)] <- "rTCA"
#reskey$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", reskey$Product)] <- "Nitrogenase"
#reskey$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", reskey$Product)] <- "Chitinase"
reskey$Category[grep("peroxidase|peroxide|catalase", reskey$Product)] <- "ROS"
#reskey$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", reskey$Product)] <- "Nitrite_reduction"
#reskey$Category[grep("urease", reskey$Product)] <- "Urease"
#reskey$Category[grep("protease", reskey$Product)] <- "Protease"
#reskey$Category[grep("sugar transport|carbohydrate ABC transport|Carbohydrate-selective porin", reskey$Product)] <- "Sugar transport"
reskey$Category[grep("raffinose/stachyose/melibiose transport", reskey$Product)] <- "R/S/M transport"
#reskey$Category[grep("chitobiose transport", reskey$Product)] <- "Chitobiose transport"
#reskey$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", reskey$Product)] <- "Chitinase"
#reskey$Category[grep("ribose transport", reskey$Product)] <- "Ribose transport"

# # Optional: use this code to get finer resolution taxonomy info
# library(dplyr)
# by_reskey <- reskey %>% group_by(Taxonomy, Category, Time)
# fine_taxa <- function(type, time){
#   x <- summarize(by_reskey, sum = sum(Counts))
#   y <- x[which(x$Category == type & x$Time == time), ]
#   return(y[order(y$sum, decreasing = T), ])
# }
# # example: fine_taxa("Sugar transport", "Day")

# Add day/night counts
gene_counts <- abun_snorm[match(reskey$Gene, rownames(abun_snorm)),]
reskey1 <- reskey
reskey2 <- reskey
reskey1$Counts <- rowSums(gene_counts[which(conditions == "day")])
reskey2$Counts <- rowSums(gene_counts[which(conditions == "night")])
reskey <- rbind(reskey1, reskey2)
reskey$Time <- c(rep("Day", dim(reskey1)[1]), rep("Night", dim(reskey2)[1]))
reskey$Counts <- as.numeric(reskey$Counts)

formatC(sum(reskey$Counts[which(reskey$Category == "ROS" & reskey$Time == "Day")]), format = "e", digits = 2)
formatC(sum(reskey$Counts[which(reskey$Category == "ROS" & reskey$Time == "Night")]), format = "e", digits = 2)

# Visually inspect phyla to keep
# ggplot(data = reskey, aes(x = Phylum, y = Counts, fill = Time)) + geom_bar(stat = "identity")  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0, 0)) + labs(x = NULL, y = "Proportion of reads", title = "Sparkling") + facet_wrap(~ Category, ncol = 3, scales="free")

# Optional: use this code to get finer resolution taxonomy info
# library(dplyr)
# by_reskey <- reskey %>% group_by(Taxonomy, Category, Time)
# fine_taxa <- function(type, time){
#   x <- summarize(by_reskey, sum = sum(Counts))
#   y <- x[which(x$Category == type & x$Time == time), ]
#   return(y[order(y$sum, decreasing = T), ])
# }
# example: fine_taxa("Sugar transport", "Day")


phyla2keep <- c("Actinobacteria", "Bacteroidetes", "Unclassified", "Cyanobacteria", "Betaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia", "Firmicutes", "Deltaproteobacteria", "Planctomycetes", "Eukaryota")
reskey <- reskey[which(reskey$Category != "None"), ] 
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))
reskey$Phylum[which(is.na(reskey$Phylum) == T)] <- "Unclassified"
reskey$Phylum <- as.character(reskey$Phylum)
reskey$Phylum[which(reskey$Phylum %not in% phyla2keep)] <- "Other"

p3 <- ggplot(data = reskey, aes(x = Time, y = Counts, fill = Phylum)) + geom_bar(stat = "identity", position = "fill")  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0, 0)) + labs(x = NULL, y = "Proportion of reads", title = "C. Sparkling Lake") + facet_wrap(~ Category, ncol = 3) + scale_fill_manual(values = as.character(colorkey$color[which(colorkey$phylum %in% unique(reskey$Phylum))])) + theme(legend.position = "none")
SPreskey <- reskey

legend_reskey <- rbind(MEreskey, TBreskey, SPreskey)
legend_plot <- ggplot(data = legend_reskey, aes(x = Time, y = Counts, fill = Phylum)) + geom_bar(stat = "identity", position = "fill")  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand = c(0, 0)) + labs(x = NULL, y = "Proportion of reads", title = "C. Sparkling") + facet_wrap(~ Category, ncol = 3) + scale_fill_manual(values = as.character(colorkey$color[which(colorkey$phylum %in% unique(legend_reskey$Phylum))]))
legend <- get_legend(legend_plot)

figures <- plot_grid(p1, p2, p3, nrow = 3, rel_heights = c(4, 4, 2.5))
to_plot <- plot_grid(figures, legend, nrow = 1, rel_widths = c(3, 1))
to_plot
save_plot("/Users/Alex/Desktop/geodes/Manuscript/figures_and_tables/taxonomy.pdf", to_plot, base_height = 8, base_aspect_ratio = 1)
