# Recreate time of maximal expression plot for Mendota

# Set up environment
path <- "/Users/Alex/Desktop/"
path2 <- "/Users/Alex/"

library(cowplot)
library(reshape2)
library(ggrepel)
library(raster)
library(DESeq2)
library(rain)
library(tidyverse)

zscore <- function(counts){
  z <- (counts - mean(counts)) / sd(counts)
  return(z)
}

taxa <- c("Actinobacteria", "Alphaproteobacteria", "Armatimonadetes", "Bacteroidetes", "Betaproteobacteria", "Chloroflexi", "Cyanobacteria", "Deltaproteobacteria", "Eukaryota", "Fusobacteria", "Gammaproteobacteria", "Planctomycetes", "Unclassified", "Verrucomicrobia", "Viruses", "Other")
colors <- c("firebrick2", "skyblue2", "lightcoral", "goldenrod", "dodgerblue", "darkgreen", "green", "steelblue", "forestgreen", "rosybrown", "royalblue4", "violet", "grey", "darkorchid", "yellow", "dimgrey")
# 
colorkey <- data.frame(taxa, colors)
# 
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

placeholder <- data.frame(genes = "x", time = c(5, 9, 13, 17, 21, 1, rep(5, 12)), max = 10, Tax = "Other", Product = "a thing", Category = c("Xylose transport", "RuBisCO", "ROS", "RNA polymerase", "Opsin", "Ribose transport", "Rhamnose transport", "Respiration", "R/S/M transport", "Protease", "Polyamines", "Photosynthesis", "Methane/Ammonia", "Lactose/arabinose transport", "Glucose/mannose transport", "General sugar transport", "Fructose transport", "Amino acid transport"), cyclic = "sort of")
# 
# Keep the top 20000 genes only
abun_mnorm <- mnorm[order(rowSums(mnorm), decreasing = T), ]
abun_mnorm <- abun_mnorm[1:20000,]
colnames(abun_mnorm) <- gsub(".nonrRNA", "", colnames(abun_mnorm))

abun_mnorm$Gene <- rownames(abun_mnorm)
abun_mnorm2 <- melt(abun_mnorm)
abun_mnorm2$Time <- metadata$Time[match(abun_mnorm2$variable, metadata$Sample)]

abun_tnorm <- tnorm[order(rowSums(tnorm), decreasing = T), ]
abun_tnorm <- abun_tnorm[1:20000,]
colnames(abun_tnorm) <- gsub(".nonrRNA", "", colnames(abun_tnorm))
abun_tnorm <- abun_tnorm[,which(colnames(abun_tnorm) %in% metadata$Sample[which(metadata$Lake == "Trout" & metadata$Timepoint <= 24)])]

abun_tnorm$Gene <- rownames(abun_tnorm)
abun_tnorm2 <- melt(abun_tnorm)
abun_tnorm2$Time <- metadata$Time[match(abun_tnorm2$variable, metadata$Sample)]


abun_snorm <- snorm[order(rowSums(snorm), decreasing = T), ]
abun_snorm <- abun_snorm[1:20000,]
colnames(abun_snorm) <- gsub(".nonrRNA", "", colnames(abun_snorm))

abun_snorm$Gene <- rownames(abun_snorm)
abun_snorm2 <- melt(abun_snorm)
abun_snorm2$Time <- metadata$Time[match(abun_snorm2$variable, metadata$Sample)]

# Aggregate by time point

aggdata <- aggregate(abun_mnorm2$value, by=list(abun_mnorm2$Gene, abun_mnorm2$Time), FUN=mean, na.rm=TRUE)

# For each gene, note its time of maximal expression
genelist <- unique(aggdata$Group.1)
cv <- c()
for(i in 1:length(genelist)){
  gene_subset <- aggdata[which(aggdata$Group.1 == genelist[i]), ]
  cv[i] <- sd(gene_subset$x)/mean(gene_subset$x)
}
genelist <- genelist[which(cv > 0.2)]

maxtime <- c()
max <- c()
for(i in 1:length(genelist)){
  gene_subset <- aggdata[which(aggdata$Group.1 == genelist[i]), ]
  maxtime[i] <- gene_subset$Group.2[which(gene_subset$x == max(gene_subset$x))]
  max[i] <- gene_subset$x[which(gene_subset$x == max(gene_subset$x))]
}

gene_data <- data.frame(genes = genelist, time = maxtime, max = max)
gene_data$Tax <- mendota_key$Phylum[match(gene_data$genes, mendota_key$Gene)]

#ggplot(gene_data, aes(x = time, y = Tax, size = max)) + geom_point()

gene_data$Product <- mendota_key$Product[match(gene_data$genes, mendota_key$Gene)]
gene_data$Category <- "None"
gene_data$Category[grep("photosystem|Photosystem|photosynth|Photosynth", gene_data$Product)] <- "Photosynthesis"
gene_data$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", gene_data$Product)] <- "Opsin"
# gene_data$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate|ose transport", gene_data$Product)] <- "Sugar degradation"
gene_data$Category[grep("rbcL|ribulose-bisphosphate carboxylase", gene_data$Product)] <- "RuBisCO"
gene_data$Category[grep("putrescine|Putrescine|spermidine|Spermidine", gene_data$Product)] <- "Polyamines"
#gene_data$Category[grep("citrate lyase|Citrate lyase", gene_data$Product)] <- "rTCA"
#gene_data$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", gene_data$Product)] <- "Nitrogenase"
#gene_data$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", gene_data$Product)] <- "Chitinase"
#gene_data$Category[grep("glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase", gene_data$Product)] <- "Glycoside_Hydrolase"
#gene_data$Category[grep("alkaline phosphatase|Alkaline phosphatase", gene_data$Product)] <- "Alkaline_phosphatase"
#gene_data$Category[grep("cellulase|cellulose", gene_data$Product)] <- "Cellulase"
gene_data$Category[grep("peroxidase|peroxide|catalase|photolyase|Photolyase", gene_data$Product)] <- "ROS"
gene_data$Category[grep("ammonia monooxygenase|methane monoxygenase", gene_data$Product)] <- "Methane/Ammonia"
#gene_data$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", gene_data$Product)] <- "Nitrite_reduction"
#gene_data$Category[grep("urease", gene_data$Product)] <- "Urease"
gene_data$Category[grep("protease", gene_data$Product)] <- "Protease"
#gene_data$Category[grep("carboxylate transport", gene_data$Product)] <- "Carboxylate transport"
#gene_data$Category[grep("ribose transport", gene_data$Product)] <- "Ribose transport"
gene_data$Category[grep("sugar transport|carbohydrate ABC transport|Carbohydrate-selective porin", gene_data$Product)] <- "General sugar transport"
gene_data$Category[grep("raffinose/stachyose/melibiose transport", gene_data$Product)] <- "R/S/M transport"
gene_data$Category[grep("glucose/mannose transport", gene_data$Product)] <- "Glucose/mannose transport"
gene_data$Category[grep("rhamnose transport", gene_data$Product)] <- "Rhamnose transport"
gene_data$Category[grep("xylose transport", gene_data$Product)] <- "Xylose transport"
gene_data$Category[grep("fructose transport", gene_data$Product)] <- "Fructose transport"
#gene_data$Category[grep("chitobiose transport", gene_data$Product)] <- "Chitobiose transport"
gene_data$Category[grep("lactose/L-arabinose transport", gene_data$Product)] <- "Lactose/arabinose transport"
gene_data$Category[grep("RNA polymerase", gene_data$Product)] <- "RNA polymerase"
#gene_data$Category[grep("DNA polymerase", gene_data$Product)] <- "DNA polymerase"
gene_data$Category[grep("cytochrome", gene_data$Product)] <- "Respiration"
#gene_data$Category[grep("sigma factor", gene_data$Product)] <- "Sigma factors"
gene_data$Category[grep("amino acid transport", gene_data$Product)] <- "Amino acid transport"
gene_data$Category[grep("eaction center M", gene_data$Product)] <- "AAP"
gene_data$Category[grep("carboxylate transport", gene_data$Product)] <- "Carboxylate transport"


#ggplot(gene_data[which(gene_data$Category != "None" & gene_data$Tax == "Cyanobacteria"), ], aes(x = time, y = Category)) + geom_jitter()

# Test cyclic trends without aggregation

abun_mnorm$Gene <- NULL
p.table <- abun_mnorm
p.results <- rain(t(p.table), deltat = 4, period = 24, measure.sequence = table(metadata$Timepoint[match(colnames(p.table), metadata$Sample)]), verbose = T, adjp.method = "Bonferroni")

sig.p.results <- p.results[which(p.results$pVal < 0.05), ]

# Add cyclic data to the max time plots
cyclic <- c()
for(i in 1:dim(gene_data)[1]){
  g <- gene_data$genes[i]
  cyclic[i] <- p.results$pVal[which(row.names(p.results) == g)] < 0.05
}
gene_data$cyclic <- cyclic

gene_data$time <- factor(gene_data$time, levels = c("5", "9", "13", "17", "21", "1"))
gene_data <- rbind(gene_data, placeholder)


p1 <- ggplot(gene_data[which(gene_data$Category != "None"), ], aes(x = time, y = Category, color = cyclic)) + geom_jitter(alpha = 0.3, size = 0.75) + scale_x_discrete(labels = c("5" = "5:00", "9" = "9:00", "13" = "13:00", "17" = "17:00", "21" = "21:00", "1" = "1:00")) + scale_color_manual(values = c("black", "white", "goldenrod")) + geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5), color = "grey") + geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), color = "grey") + labs(y = "", x = "", title = "Mendota") + theme(legend.position = "none")

# Make taxonomy plots of key regulators
aggdata$Phylum <- as.character(mendota_key$Phylum[match(aggdata$Group.1, mendota_key$Gene)])
aggdata$Phylum[which(is.na(aggdata$Phylum) == T)] <- "Unclassified"
aggdata$Category <- gene_data$Category[match(aggdata$Group.1, gene_data$genes)]
aggdata$Group.2 <- factor(aggdata$Group.2, levels = c("1", "5", "9", "13", "17", "21"))

rna_gene_data <- aggdata[which(aggdata$Category == "RNA polymerase"), ]
tax_counts <- table(rna_gene_data$Phylum)
grey_out <- names(tax_counts)[which(tax_counts < 15)]
toss <- match(rna_gene_data$Phylum, grey_out)
rna_gene_data$Phylum[which(is.na(toss) == F)] <- "Other"
rna_legend <- rna_gene_data
tax1.1 <- ggplot(rna_gene_data, aes(x = Group.2, y = x, fill = Phylum)) + geom_bar(stat = "identity")+ labs(x = NULL, y = NULL) + theme(legend.position = "none") + scale_fill_manual(values = as.character(colorkey$colors[which(colorkey$taxa %in% unique(rna_gene_data$Phylum))]))

rna_gene_data <- aggdata[which(aggdata$Category == "Photosynthesis"), ]
tax_counts <- table(rna_gene_data$Phylum)
grey_out <- names(tax_counts)[which(tax_counts < 15)]
toss <- match(rna_gene_data$Phylum, grey_out)
rna_gene_data$Phylum[which(is.na(toss) == F)] <- "Other"
rna_legend <- rbind(rna_legend, rna_gene_data)
tax1.2 <- ggplot(rna_gene_data, aes(x = Group.2, y = x, fill = Phylum)) + geom_bar(stat = "identity")+ labs(x = NULL, y = NULL) + theme(legend.position = "none") + scale_fill_manual(values = as.character(colorkey$colors[which(colorkey$taxa %in% unique(rna_gene_data$Phylum))]))

rna_gene_data <- aggdata[which(aggdata$Category == "General sugar transport"), ]
tax_counts <- table(rna_gene_data$Phylum)
grey_out <- names(tax_counts)[which(tax_counts < 15)]
toss <- match(rna_gene_data$Phylum, grey_out)
rna_gene_data$Phylum[which(is.na(toss) == F)] <- "Other"
rna_legend <- rbind(rna_legend, rna_gene_data)
tax1.3 <- ggplot(rna_gene_data, aes(x = Group.2, y = x, fill = Phylum)) + geom_bar(stat = "identity")+ labs(x = NULL, y = NULL) + theme(legend.position = "none") + scale_fill_manual(values = as.character(colorkey$colors[which(colorkey$taxa %in% unique(rna_gene_data$Phylum))]))

# Repeat for other lakes

aggdata <- aggregate(abun_snorm2$value, by=list(abun_snorm2$Gene, abun_snorm2$Time), FUN=mean, na.rm=TRUE)

# For each gene, note its time of maximal expression
genelist <- unique(aggdata$Group.1)
cv <- c()
for(i in 1:length(genelist)){
  gene_subset <- aggdata[which(aggdata$Group.1 == genelist[i]), ]
  cv[i] <- sd(gene_subset$x)/mean(gene_subset$x)
}
genelist <- genelist[which(cv > 0.2)]
maxtime <- c()
max <- c()
for(i in 1:length(genelist)){
  gene_subset <- aggdata[which(aggdata$Group.1 == genelist[i]), ]
  maxtime[i] <- gene_subset$Group.2[which(gene_subset$x == max(gene_subset$x))]
  max[i] <- gene_subset$x[which(gene_subset$x == max(gene_subset$x))]
}

gene_data <- data.frame(genes = genelist, time = maxtime, max = max)
gene_data$Tax <- spark_key$Phylum[match(gene_data$genes, spark_key$Gene)]

#ggplot(gene_data, aes(x = time, y = Tax, size = max)) + geom_point()

gene_data$Product <- spark_key$Product[match(gene_data$genes, spark_key$Gene)]
gene_data$Category <- "None"
gene_data$Category[grep("photosystem|Photosystem|photosynth|Photosynth", gene_data$Product)] <- "Photosynthesis"
gene_data$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", gene_data$Product)] <- "Opsin"
# gene_data$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate|ose transport", gene_data$Product)] <- "Sugar degradation"
gene_data$Category[grep("rbcL|ribulose-bisphosphate carboxylase", gene_data$Product)] <- "RuBisCO"
gene_data$Category[grep("putrescine|Putrescine|spermidine|Spermidine", gene_data$Product)] <- "Polyamines"
#gene_data$Category[grep("citrate lyase|Citrate lyase", gene_data$Product)] <- "rTCA"
#gene_data$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", gene_data$Product)] <- "Nitrogenase"
#gene_data$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", gene_data$Product)] <- "Chitinase"
#gene_data$Category[grep("glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase", gene_data$Product)] <- "Glycoside_Hydrolase"
#gene_data$Category[grep("alkaline phosphatase|Alkaline phosphatase", gene_data$Product)] <- "Alkaline_phosphatase"
#gene_data$Category[grep("cellulase|cellulose", gene_data$Product)] <- "Cellulase"
gene_data$Category[grep("peroxidase|peroxide|catalase|photolyase|Photolyase", gene_data$Product)] <- "ROS"
gene_data$Category[grep("ammonia monooxygenase|methane monoxygenase", gene_data$Product)] <- "Methane/Ammonia"
#gene_data$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", gene_data$Product)] <- "Nitrite_reduction"
#gene_data$Category[grep("urease", gene_data$Product)] <- "Urease"
gene_data$Category[grep("protease", gene_data$Product)] <- "Protease"
#gene_data$Category[grep("carboxylate transport", gene_data$Product)] <- "Carboxylate transport"
gene_data$Category[grep("ribose transport", gene_data$Product)] <- "Ribose transport"
gene_data$Category[grep("sugar transport|carbohydrate ABC transport|Carbohydrate-selective porin", gene_data$Product)] <- "General sugar transport"
gene_data$Category[grep("raffinose/stachyose/melibiose transport", gene_data$Product)] <- "R/S/M transport"
gene_data$Category[grep("glucose/mannose transport", gene_data$Product)] <- "Glucose/mannose transport"
gene_data$Category[grep("rhamnose transport", gene_data$Product)] <- "Rhamnose transport"
gene_data$Category[grep("xylose transport", gene_data$Product)] <- "Xylose transport"
gene_data$Category[grep("fructose transport", gene_data$Product)] <- "Fructose transport"
#gene_data$Category[grep("chitobiose transport", gene_data$Product)] <- "Chitobiose transport"
gene_data$Category[grep("lactose/L-arabinose transport", gene_data$Product)] <- "Lactose/arabinose transport"
gene_data$Category[grep("RNA polymerase", gene_data$Product)] <- "RNA polymerase"
#gene_data$Category[grep("DNA polymerase", gene_data$Product)] <- "DNA polymerase"
gene_data$Category[grep("cytochrome", gene_data$Product)] <- "Respiration"
#gene_data$Category[grep("sigma factor", gene_data$Product)] <- "Sigma factors"
gene_data$Category[grep("amino acid transport", gene_data$Product)] <- "Amino acid transport"
gene_data$Category[grep("eaction center M", gene_data$Product)] <- "AAP"
gene_data$Category[grep("carboxylate transport", gene_data$Product)] <- "Carboxylate transport"

# Test cyclic trends without aggregation

abun_snorm$Gene <- NULL
p.table <- abun_snorm
p.results <- rain(t(p.table), deltat = 4, period = 24, measure.sequence = table(metadata$Timepoint[match(colnames(p.table), metadata$Sample)]), verbose = T, adjp.method = "Bonferroni")

sig.p.results <- p.results[which(p.results$pVal < 0.05), ]

# Add cyclic data to the max time plots
cyclic <- c()
for(i in 1:dim(gene_data)[1]){
  g <- gene_data$genes[i]
  cyclic[i] <- p.results$pVal[which(row.names(p.results) == g)] < 0.05
}
gene_data$cyclic <- cyclic

gene_data$time <- factor(gene_data$time, levels = c("5", "9", "13", "17", "21", "1"))
gene_data <- rbind(gene_data, placeholder)

p2 <- ggplot(gene_data[which(gene_data$Category != "None"), ], aes(x = time, y = Category, color = cyclic)) + geom_jitter(alpha = 0.3, size = 0.75) + scale_x_discrete(labels = c("5" = "5:00", "9" = "9:00", "13" = "13:00", "17" = "17:00", "21" = "21:00", "1" = "1:00")) + scale_color_manual(values = c("black", "white", "goldenrod")) + geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5), color = "grey") + geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), color = "grey") + labs(title = "Sparkling", x = "", y = "") + theme(legend.position = "none")

# Make taxonomy plots of key regulators
aggdata$Phylum <- as.character(spark_key$Phylum[match(aggdata$Group.1, spark_key$Gene)])
aggdata$Phylum[which(is.na(aggdata$Phylum) == T)] <- "Unclassified"
aggdata$Category <- gene_data$Category[match(aggdata$Group.1, gene_data$genes)]
aggdata$Group.2 <- factor(aggdata$Group.2, levels = c("1", "5", "9", "13", "17", "21"))

rna_gene_data <- aggdata[which(aggdata$Category == "RNA polymerase"), ]
tax_counts <- table(rna_gene_data$Phylum)
grey_out <- names(tax_counts)[which(tax_counts < 15)]
toss <- match(rna_gene_data$Phylum, grey_out)
rna_gene_data$Phylum[which(is.na(toss) == F)] <- "Other"
rna_legend <- rbind(rna_legend, rna_gene_data)
tax2.1 <- ggplot(rna_gene_data, aes(x = Group.2, y = x, fill = Phylum)) + geom_bar(stat = "identity")+ labs(x = NULL, y = NULL) + theme(legend.position = "none") + scale_fill_manual(values = as.character(colorkey$colors[which(colorkey$taxa %in% unique(rna_gene_data$Phylum))]))

rna_gene_data <- aggdata[which(aggdata$Category == "Photosynthesis"), ]
tax_counts <- table(rna_gene_data$Phylum)
grey_out <- names(tax_counts)[which(tax_counts < 15)]
toss <- match(rna_gene_data$Phylum, grey_out)
rna_gene_data$Phylum[which(is.na(toss) == F)] <- "Other"
rna_legend <- rbind(rna_legend, rna_gene_data)
tax2.2 <- ggplot(rna_gene_data, aes(x = Group.2, y = x, fill = Phylum)) + geom_bar(stat = "identity")+ labs(x = NULL, y = NULL) + theme(legend.position = "none") + scale_fill_manual(values = as.character(colorkey$colors[which(colorkey$taxa %in% unique(rna_gene_data$Phylum))]))

rna_gene_data <- aggdata[which(aggdata$Category == "General sugar transport"), ]
tax_counts <- table(rna_gene_data$Phylum)
grey_out <- names(tax_counts)[which(tax_counts < 15)]
toss <- match(rna_gene_data$Phylum, grey_out)
rna_gene_data$Phylum[which(is.na(toss) == F)] <- "Other"
rna_legend <- rbind(rna_legend, rna_gene_data)
tax2.3 <- ggplot(rna_gene_data, aes(x = Group.2, y = x, fill = Phylum)) + geom_bar(stat = "identity")+ labs(x = NULL, y = NULL) + theme(legend.position = "none") + scale_fill_manual(values = as.character(colorkey$colors[which(colorkey$taxa %in% unique(rna_gene_data$Phylum))]))

# Repeat for other lakes

aggdata <- aggregate(abun_tnorm2$value, by=list(abun_tnorm2$Gene, abun_tnorm2$Time), FUN=mean, na.rm=TRUE)

cv <- c()
for(i in 1:length(genelist)){
  gene_subset <- aggdata[which(aggdata$Group.1 == genelist[i]), ]
  cv[i] <- sd(gene_subset$x)/mean(gene_subset$x)
}
genelist <- genelist[which(cv > 0.2)]

# For each gene, note its time of maximal expression
genelist <- unique(aggdata$Group.1)
maxtime <- c()
max <- c()
for(i in 1:length(genelist)){
  gene_subset <- aggdata[which(aggdata$Group.1 == genelist[i]), ]
  maxtime[i] <- gene_subset$Group.2[which(gene_subset$x == max(gene_subset$x))]
  max[i] <- gene_subset$x[which(gene_subset$x == max(gene_subset$x))]
}

gene_data <- data.frame(genes = genelist, time = maxtime, max = max)
gene_data$Tax <- mendota_key$Phylum[match(gene_data$genes, mendota_key$Gene)]

#ggplot(gene_data, aes(x = time, y = Tax, size = max)) + geom_point()

gene_data$Product <- mendota_key$Product[match(gene_data$genes, mendota_key$Gene)]
gene_data$Category <- "None"
gene_data$Category[grep("photosystem|Photosystem|photosynth|Photosynth", gene_data$Product)] <- "Photosynthesis"
gene_data$Category[grep("rhodopsin|Rhodopsin|phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene", gene_data$Product)] <- "Opsin"
# gene_data$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate|ose transport", gene_data$Product)] <- "Sugar degradation"
gene_data$Category[grep("rbcL|ribulose-bisphosphate carboxylase", gene_data$Product)] <- "RuBisCO"
gene_data$Category[grep("putrescine|Putrescine|spermidine|Spermidine", gene_data$Product)] <- "Polyamines"
#gene_data$Category[grep("citrate lyase|Citrate lyase", gene_data$Product)] <- "rTCA"
#gene_data$Category[grep("nitrogenase|Nitrogenase|NifH|NifD|NifK", gene_data$Product)] <- "Nitrogenase"
#gene_data$Category[grep("Chitobiase|chitobiase|chitinase|Chitinase", gene_data$Product)] <- "Chitinase"
#gene_data$Category[grep("glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase", gene_data$Product)] <- "Glycoside_Hydrolase"
#gene_data$Category[grep("alkaline phosphatase|Alkaline phosphatase", gene_data$Product)] <- "Alkaline_phosphatase"
#gene_data$Category[grep("cellulase|cellulose", gene_data$Product)] <- "Cellulase"
gene_data$Category[grep("peroxidase|peroxide|catalase|photolyase|Photolyase", gene_data$Product)] <- "ROS"
gene_data$Category[grep("ammonia monooxygenase|methane monoxygenase", gene_data$Product)] <- "Methane/Ammonia"
#gene_data$Category[grep("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase", gene_data$Product)] <- "Nitrite_reduction"
#gene_data$Category[grep("urease", gene_data$Product)] <- "Urease"
gene_data$Category[grep("protease", gene_data$Product)] <- "Protease"
#gene_data$Category[grep("carboxylate transport", gene_data$Product)] <- "Carboxylate transport"
gene_data$Category[grep("ribose transport", gene_data$Product)] <- "Ribose transport"
gene_data$Category[grep("sugar transport|carbohydrate ABC transport|Carbohydrate-selective porin", gene_data$Product)] <- "General sugar transport"
gene_data$Category[grep("raffinose/stachyose/melibiose transport", gene_data$Product)] <- "R/S/M transport"
gene_data$Category[grep("glucose/mannose transport", gene_data$Product)] <- "Glucose/mannose transport"
gene_data$Category[grep("rhamnose transport", gene_data$Product)] <- "Rhamnose transport"
gene_data$Category[grep("xylose transport", gene_data$Product)] <- "Xylose transport"
gene_data$Category[grep("fructose transport", gene_data$Product)] <- "Fructose transport"
#gene_data$Category[grep("chitobiose transport", gene_data$Product)] <- "Chitobiose transport"
gene_data$Category[grep("lactose/L-arabinose transport", gene_data$Product)] <- "Lactose/arabinose transport"
gene_data$Category[grep("RNA polymerase", gene_data$Product)] <- "RNA polymerase"
#gene_data$Category[grep("DNA polymerase", gene_data$Product)] <- "DNA polymerase"
gene_data$Category[grep("cytochrome", gene_data$Product)] <- "Respiration"
#gene_data$Category[grep("sigma factor", gene_data$Product)] <- "Sigma factors"
gene_data$Category[grep("amino acid transport", gene_data$Product)] <- "Amino acid transport"
gene_data$Category[grep("eaction center M", gene_data$Product)] <- "AAP"
gene_data$Category[grep("carboxylate transport", gene_data$Product)] <- "Carboxylate transport"

# Test cyclic trends without aggregation

abun_tnorm$Gene <- NULL
p.table <- abun_tnorm
p.results <- rain(t(p.table), deltat = 4, period = 24, measure.sequence = table(metadata$Timepoint[match(colnames(p.table), metadata$Sample)]), verbose = T, adjp.method = "Bonferroni")

sig.p.results <- p.results[which(p.results$pVal < 0.05), ]

# Add cyclic data to the max time plots
cyclic <- c()
for(i in 1:dim(gene_data)[1]){
  g <- gene_data$genes[i]
  cyclic[i] <- p.results$pVal[which(row.names(p.results) == g)] < 0.05
}
gene_data$cyclic <- cyclic

gene_data$time <- factor(gene_data$time, levels = c("5", "9", "13", "17", "21", "1"))

gene_data <- rbind(gene_data, placeholder)

p3 <- ggplot(gene_data[which(gene_data$Category != "None"), ], aes(x = time, y = Category, color = cyclic)) + geom_jitter(alpha = 0.3, size = 0.75) + scale_x_discrete(labels = c("5" = "5:00", "9" = "9:00", "13" = "13:00", "17" = "17:00", "21" = "21:00", "1" = "1:00")) + scale_color_manual(values = c("black", "white", "goldenrod")) + geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5), color = "grey") + geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), color = "grey") + labs(title = "Trout", x = "", y = "") + theme(legend.position = "none")

dummydata <- data.frame(genes = "x", time = c(17, 17, 21), max = 10, Tax = "Zed", Product = "a thing", Category = c("RNA polymerase", "General sugar degradation", "General sugar degradation"), cyclic = "sort of")
# Make taxonomy plots of key regulators
aggdata$Phylum <- as.character(trout_key$Phylum[match(aggdata$Group.1, trout_key$Gene)])
aggdata$Phylum[which(is.na(aggdata$Phylum) == T)] <- "Unclassified"
aggdata$Category <- gene_data$Category[match(aggdata$Group.1, gene_data$genes)]
aggdata$Group.2 <- factor(aggdata$Group.2, levels = c("1", "5", "9", "13", "17", "21"))

rna_gene_data <- aggdata[which(aggdata$Category == "RNA polymerase"), ]
tax_counts <- table(rna_gene_data$Phylum)
grey_out <- names(tax_counts)[which(tax_counts < 15)]
toss <- match(rna_gene_data$Phylum, grey_out)
rna_gene_data$Phylum[which(is.na(toss) == F)] <- "Other"
rna_legend <- rbind(rna_legend, rna_gene_data)
tax3.1 <- ggplot(rna_gene_data, aes(x = Group.2, y = x, fill = Phylum)) + geom_bar(stat = "identity")+ labs(x = NULL, y = NULL) + theme(legend.position = "none") + scale_fill_manual(values = as.character(colorkey$colors[which(colorkey$taxa %in% unique(rna_gene_data$Phylum))]))

rna_gene_data <- aggdata[which(aggdata$Category == "Photosynthesis"), ]
tax_counts <- table(rna_gene_data$Phylum)
grey_out <- names(tax_counts)[which(tax_counts < 15)]
toss <- match(rna_gene_data$Phylum, grey_out)
rna_gene_data$Phylum[which(is.na(toss) == F)] <- "Other"
rna_legend <- rbind(rna_legend, rna_gene_data)
tax3.2 <- ggplot(rna_gene_data, aes(x = Group.2, y = x, fill = Phylum)) + geom_bar(stat = "identity")+ labs(x = NULL, y = NULL) + theme(legend.position = "none") + scale_fill_manual(values = as.character(colorkey$colors[which(colorkey$taxa %in% unique(rna_gene_data$Phylum))]))

rna_gene_data <- aggdata[which(aggdata$Category == "General sugar transport"), ]
tax_counts <- table(rna_gene_data$Phylum)
grey_out <- names(tax_counts)[which(tax_counts < 15)]
toss <- match(rna_gene_data$Phylum, grey_out)
rna_gene_data$Phylum[which(is.na(toss) == F)] <- "Other"
rna_legend <- rbind(rna_legend, rna_gene_data)
tax3.3 <- ggplot(rna_gene_data, aes(x = Group.2, y = x, fill = Phylum)) + geom_bar(stat = "identity")+ labs(x = NULL, y = NULL) + theme(legend.position = "none") + scale_fill_manual(values = as.character(colorkey$colors[which(colorkey$taxa %in% unique(rna_gene_data$Phylum))]))



to_plot <- plot_grid(p1, p2, p3, nrow = 1)
#save_plot(paste(path,"geodes/Manuscript/figures_and_tables/maxtime.pdf", sep = ""),  to_plot, base_aspect_ratio = 4, base_height = 4)

to_plot2 <- plot_grid(tax1.1, tax1.2, tax1.3, tax2.1, tax2.2, tax2.3, tax3.1, tax3.2, tax3.3, nrow = 3)

legendplot <- ggplot(rna_legend, aes(x = Group.2, y = x, fill = Phylum)) + geom_bar(stat = "identity") + labs(x = NULL, y = NULL)  + scale_fill_manual(values = as.character(colorkey$colors[which(colorkey$taxa %in% unique(rna_legend$Phylum))]))
legend <- get_legend(legendplot)
to_plot2 <- plot_grid(to_plot2, legend, nrow = 1, rel_widths = c(5, 1))
to_plot2

#save_plot(paste(path,"geodes/Manuscript/figures_and_tables/taxatime.pdf", sep = ""),  to_plot2, base_aspect_ratio = 3, base_height = 4)
