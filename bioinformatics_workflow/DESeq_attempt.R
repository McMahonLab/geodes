library(ggplot2)
library(cowplot)
library(reshape2)
library(DESeq)
library(GeneCycle)

#path <- "D:/"
path <- "C:/Users/Goose and Gander/Documents/"

snorm <- read.csv(paste(path, "geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)
tnorm <- read.csv(paste(path, "geodes_data_tables/Trout_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)
mnorm <- read.csv(paste(path, "geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)

snorm <- snorm[, which(colnames(snorm) != "GEODES014.nonrRNA" & colnames(snorm) != "GEODES033.nonrRNA")]
tnorm <- tnorm[, which(colnames(tnorm) != "GEODES065.nonrRNA")]
mnorm <- mnorm[, which(colnames(mnorm) != "GEODES158.nonrRNA")]

mendota_key <- read.csv(paste(path, "geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", sep = ""), header = T)
spark_key <- read.csv(paste(path, "geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", sep = ""), header = T)
trout_key <- read.csv(paste(path, "geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", sep = ""), header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# Take the top 20000 of each data table

abun_mnorm <- mnorm[order(rowSums(mnorm), decreasing = T), ]
#Nophoto
mendota_key <- mendota_key[match(rownames(abun_mnorm), mendota_key$Gene), ]
abun_mnorm <- abun_mnorm[grep("Cyano", mendota_key$Taxonomy, invert = T), ]
mendota_key <- mendota_key[match(rownames(abun_mnorm), mendota_key$Gene), ]
abun_mnorm <- abun_mnorm[grep("photo|Photo|chlorophyll|rhodopsin", mendota_key$Product, invert = T), ]
abun_mnorm <- abun_mnorm[1:20000,]


abun_mnorm$Genes <- rownames(abun_mnorm)
abun_mnorm <- melt(abun_mnorm)
abun_mnorm$variable <- gsub(".nonrRNA", "", abun_mnorm$variable)
abun_mnorm$Timepoint <- metadata$Timepoint[match(abun_mnorm$variable, metadata$Sample)]
agg_abun_mnorm <- aggregate(value ~ Genes + Timepoint, data = abun_mnorm, mean)
new_abun_mnorm <- reshape(agg_abun_mnorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(new_abun_mnorm) <- new_abun_mnorm[, 1]
new_abun_mnorm <- new_abun_mnorm[, 2:dim(new_abun_mnorm)[2]]


# abun_snorm$Genes <- rownames(abun_snorm)
# abun_snorm <- melt(abun_snorm)
# abun_snorm$variable <- gsub(".nonrRNA", "", abun_snorm$variable)
# abun_snorm$Timepoint <- metadata$Timepoint[match(abun_snorm$variable, metadata$Sample)]
# agg_abun_snorm <- aggregate(value ~ Genes + Timepoint, data = abun_snorm, mean)
# new_abun_snorm <- reshape(agg_abun_snorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
# rownames(new_abun_snorm) <- new_abun_snorm[, 1]
# new_abun_snorm <- new_abun_snorm[, 2:dim(new_abun_snorm)[2]]


# abun_tnorm$Genes <- rownames(abun_tnorm)
# abun_tnorm <- melt(abun_tnorm)
# abun_tnorm$variable <- gsub(".nonrRNA", "", abun_tnorm$variable)
# abun_tnorm$Timepoint <- metadata$Timepoint[match(abun_tnorm$variable, metadata$Sample)]
# agg_abun_tnorm <- aggregate(value ~ Genes + Timepoint, data = abun_tnorm, mean)
# new_abun_tnorm <- reshape(agg_abun_tnorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
# rownames(new_abun_tnorm) <- new_abun_tnorm[, 1]
# new_abun_tnorm <- new_abun_tnorm[, 2:dim(new_abun_tnorm)[2]]
# new_abun_tnorm <- new_abun_tnorm[, 1:7]


# Then run the significance test for cyclic trends
# fdr.spark <- fdrtool(fisher.g.test(t(new_abun_snorm)), statistic = "pvalue")
# length(which(fdr.spark$pval < 0.05))/length(fdr.spark$pval)
# quantile(fdr.spark$qval)

# fdr.trout <- fdrtool(fisher.g.test(t(new_abun_tnorm)), statistic = "pvalue")
# length(which(fdr.trout$pval < 0.05))/length(fdr.trout$pval)
# quantile(fdr.trout$qval)

fdr.mendota <- fdrtool(fisher.g.test(t(new_abun_mnorm)), statistic = "pvalue")
length(which(fdr.mendota$pval < 0.05))/length(fdr.mendota$pval)
quantile(fdr.mendota$qval)

sig.mendota <- t(new_abun_mnorm[which(fdr.mendota$qval < 0.05),])
sig.mendota.key <- mendota_key[match(colnames(sig.mendota), mendota_key$Gene), ]
x <- table(as.character(sig.mendota.key[,3]))

# Do any taxa have more significantly cyclic genes than expected by chance?
# Expected value is 1584/10000 * 100 = 16%

subset <- c()
for(i in 1:length(x)){
  before <- length(grep(paste(names(x)[i], "$", sep = ""), mendota_key$Taxonomy))
  after <- as.numeric(x[i])
  subset[i] <- after/before *100
}

# For each gene, record the max time of expression and its taxonomic classification

sig.genes <- sig.mendota.key[,c(1, 3:4)]
sig.genes$Taxonomy <- gsub(";;;", "", sig.genes$Taxonomy)
sig.genes$Taxonomy <- gsub(";;", "", sig.genes$Taxonomy)
sig.genes$Taxonomy <- gsub(";$", "", sig.genes$Taxonomy)
spl <- strsplit(as.character(sig.genes$Taxonomy), ";")
sig.genes$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")

max.time <- c()
max.value <- c()
for(i in 1:dim(sig.genes)[1]){
  genedata <- new_abun_mnorm[which(rownames(new_abun_mnorm) == as.character(sig.genes$Gene[i])), ]
  max.time[i] <- colnames(genedata)[which(genedata == max(genedata))]
  max.value[i] <- max(genedata)
}

sig.genes$Timepoint <- max.time
sig.genes$Time <- max.time
sig.genes$Time[which(sig.genes$Time == "value.0" | sig.genes$Time == "value.24")] <- "5AM"
sig.genes$Time[which(sig.genes$Time == "value.4" | sig.genes$Time == "value.28")] <- "9AM"
sig.genes$Time[which(sig.genes$Time == "value.8" | sig.genes$Time == "value.32")] <- "1PM"
sig.genes$Time[which(sig.genes$Time == "value.12" | sig.genes$Time == "value.36")] <- "5PM"
sig.genes$Time[which(sig.genes$Time == "value.16" | sig.genes$Time == "value.40")] <- "9PM"
sig.genes$Time[which(sig.genes$Time == "value.20" | sig.genes$Time == "value.44")] <- "1AM"
sig.genes$Time <- factor(sig.genes$Time, levels = c("5AM", "9AM", "1PM", "5PM", "9PM", "1AM"))
sig.genes$Value <- max.value

sig.taxa <- aggregate(Value ~ Time + ShortTax, sig.genes, sum)
sig.taxa <- sig.taxa[grep("NO CLASSIFICATION", sig.taxa$ShortTax, invert = T), ]
sig.taxa <- sig.taxa[which(sig.taxa$Value > 30000000), ]

# Must have at least 10 diel genes
sig.taxa <- sig.taxa[which(sig.taxa$ShortTax == "Actinobacteria" | sig.taxa$ShortTax == "Bacteroidetes" | sig.taxa$ShortTax == "Bdellovibrionales" | sig.taxa$ShortTax == "Comamonadaceae" | sig.taxa$ShortTax == "Flavobacteriales" | sig.taxa$ShortTax == "Gemmatimonas"),]

ggplot(sig.taxa, aes(x = Time, y = ShortTax, color = ShortTax, size = log(Value))) + geom_point()


# Skip back to DESeq
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
sig.res <- res[which(res$padj < 0.05), ]
sig.res.day <- sig.res[which(sig.res$log2FoldChange < 0), ]
sig.res.night <- sig.res[which(sig.res$log2FoldChange > 0), ]
sig.res.day.key <- mendota_key[match(rownames(sig.res.day), mendota_key$Gene),]
sig.res.night.key <- mendota_key[match(rownames(sig.res.night), mendota_key$Gene),]

# Make plots of this aggreagted by short taxa and groupd products

sig.res.day.key$Taxonomy <- gsub(";;;", "", sig.res.day.key$Taxonomy)
sig.res.day.key$Taxonomy <- gsub(";;", "", sig.res.day.key$Taxonomy)
sig.res.day.key$Taxonomy <- gsub(";$", "", sig.res.day.key$Taxonomy)
spl <- strsplit(as.character(sig.res.day.key$Taxonomy), ";")
sig.res.day.key$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")

sig.res.night.key$Taxonomy <- gsub(";;;", "", sig.res.night.key$Taxonomy)
sig.res.night.key$Taxonomy <- gsub(";;", "", sig.res.night.key$Taxonomy)
sig.res.night.key$Taxonomy <- gsub(";$", "", sig.res.night.key$Taxonomy)
spl <- strsplit(as.character(sig.res.night.key$Taxonomy), ";")
sig.res.night.key$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")

# Make column by certain key words 
sig.res.day.key$Category <- "None"
sig.res.day.key$Category[grep("photo|Photo", sig.res.day.key$Product)] <- "Photosynthesis"
sig.res.day.key$Category[grep("rhodopsin|Rhodopsin", sig.res.day.key$Product)] <- "Rhodopsin"
sig.res.day.key$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate", sig.res.day.key$Product)] <- "Sugar degradation"
sig.res.day.key$Category[grep("rbcL|ribulose-bisphosphate carboxylase", sig.res.day.key$Product)] <- "RubisCO"
sig.res.day.key$Category[grep("putrescine|Putrescine|spermidine|Spermidine", sig.res.day.key$Product)] <- "Polyamines"

sig.res.night.key$Category <- "None"
sig.res.night.key$Category[grep("photo|Photo", sig.res.night.key$Product)] <- "Photosynthesis"
sig.res.night.key$Category[grep("rhodopsin|Rhodopsin", sig.res.night.key$Product)] <- "Rhodopsin"
sig.res.night.key$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate", sig.res.night.key$Product)] <- "Sugar degradation"
sig.res.night.key$Category[grep("rbcL|ribulose-bisphosphate carboxylase", sig.res.night.key$Product)] <- "RubisCO"
sig.res.night.key$Category[grep("putrescine|Putrescine|spermidine|Spermidine", sig.res.night.key$Product)] <- "Polyamines"

sig.res.day.key$Condition <- "day"
sig.res.night.key$Condition <- "night"

sig.res.key <- rbind(sig.res.day.key, sig.res.night.key)

ggplot(data = sig.res.key, aes(x = Category, y = totals, group = ShortTax)) + geom_bar(stat = "identity")

ggplot(data = sig.res.key[which(sig.res.key$Category != "None"), ], aes(x = Category, y = totals, fill = Condition)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "Mendota")
ggplot(data = sig.res.key[which(sig.res.key$Category != "None" & sig.res.key$Category != "Photosynthesis"), ], aes(x = Category, y = totals, fill = Condition)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "Mendota")

# Plug old marker gene analysis in here and repeat by lake