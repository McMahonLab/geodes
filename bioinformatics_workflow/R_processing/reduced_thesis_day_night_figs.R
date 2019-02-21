library(ggplot2)
library(cowplot)
library(reshape2)
library(DESeq2)
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
# mendota_key <- mendota_key[match(rownames(abun_mnorm), mendota_key$Gene), ]
# abun_mnorm <- abun_mnorm[grep("Cyano", mendota_key$Taxonomy, invert = T), ]
# mendota_key <- mendota_key[match(rownames(abun_mnorm), mendota_key$Gene), ]
# abun_mnorm <- abun_mnorm[grep("photo|Photo|chlorophyll|rhodopsin", mendota_key$Product, invert = T), ]
abun_mnorm <- abun_mnorm[1:20000,]
mendota_key <- mendota_key[match(rownames(abun_mnorm), mendota_key$Gene), ]

abun_tnorm <- tnorm[order(rowSums(tnorm), decreasing = T), ]
abun_tnorm <- abun_tnorm[1:20000,]
trout_key <- trout_key[match(rownames(abun_tnorm), trout_key$Gene), ]

abun_snorm <- snorm[order(rowSums(snorm), decreasing = T), ]
abun_snorm <- abun_snorm[1:20000,]
spark_key <- spark_key[match(rownames(abun_snorm), spark_key$Gene), ]

# abun_mnorm$Genes <- rownames(abun_mnorm)
# abun_mnorm <- melt(abun_mnorm)
# abun_mnorm$variable <- gsub(".nonrRNA", "", abun_mnorm$variable)
# abun_mnorm$Timepoint <- metadata$Timepoint[match(abun_mnorm$variable, metadata$Sample)]
# agg_abun_mnorm <- aggregate(value ~ Genes + Timepoint, data = abun_mnorm, mean)
# new_abun_mnorm <- reshape(agg_abun_mnorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
# rownames(new_abun_mnorm) <- new_abun_mnorm[, 1]
# new_abun_mnorm <- new_abun_mnorm[, 2:dim(new_abun_mnorm)[2]]


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

# fdr.mendota <- fdrtool(fisher.g.test(t(new_abun_mnorm)), statistic = "pvalue")
# length(which(fdr.mendota$pval < 0.05))/length(fdr.mendota$pval)
# quantile(fdr.mendota$qval)
# 
# sig.mendota <- t(new_abun_mnorm[which(fdr.mendota$qval < 0.05),])
# sig.mendota.key <- mendota_key[match(colnames(sig.mendota), mendota_key$Gene), ]
# x <- table(as.character(sig.mendota.key[,3]))
# 
# # Do any taxa have more significantly cyclic genes than expected by chance?
# # Expected value is 1584/10000 * 100 = 16%
# 
# subset <- c()
# for(i in 1:length(x)){
#   before <- length(grep(paste(names(x)[i], "$", sep = ""), mendota_key$Taxonomy))
#   after <- as.numeric(x[i])
#   subset[i] <- after/before *100
# }
# 
# # For each gene, record the max time of expression and its taxonomic classification
# 
# sig.genes <- sig.mendota.key[,c(1, 3:4)]
# sig.genes$Taxonomy <- gsub(";;;", "", sig.genes$Taxonomy)
# sig.genes$Taxonomy <- gsub(";;", "", sig.genes$Taxonomy)
# sig.genes$Taxonomy <- gsub(";$", "", sig.genes$Taxonomy)
# spl <- strsplit(as.character(sig.genes$Taxonomy), ";")
# sig.genes$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")
# 
# max.time <- c()
# max.value <- c()
# for(i in 1:dim(sig.genes)[1]){
#   genedata <- new_abun_mnorm[which(rownames(new_abun_mnorm) == as.character(sig.genes$Gene[i])), ]
#   max.time[i] <- colnames(genedata)[which(genedata == max(genedata))]
#   max.value[i] <- max(genedata)
# }
# 
# sig.genes$Timepoint <- max.time
# sig.genes$Time <- max.time
# sig.genes$Time[which(sig.genes$Time == "value.0" | sig.genes$Time == "value.24")] <- "5AM"
# sig.genes$Time[which(sig.genes$Time == "value.4" | sig.genes$Time == "value.28")] <- "9AM"
# sig.genes$Time[which(sig.genes$Time == "value.8" | sig.genes$Time == "value.32")] <- "1PM"
# sig.genes$Time[which(sig.genes$Time == "value.12" | sig.genes$Time == "value.36")] <- "5PM"
# sig.genes$Time[which(sig.genes$Time == "value.16" | sig.genes$Time == "value.40")] <- "9PM"
# sig.genes$Time[which(sig.genes$Time == "value.20" | sig.genes$Time == "value.44")] <- "1AM"
# sig.genes$Time <- factor(sig.genes$Time, levels = c("5AM", "9AM", "1PM", "5PM", "9PM", "1AM"))
# sig.genes$Value <- max.value
# 
# sig.taxa <- aggregate(Value ~ Time + ShortTax, sig.genes, sum)
# sig.taxa <- sig.taxa[grep("NO CLASSIFICATION", sig.taxa$ShortTax, invert = T), ]
# sig.taxa <- sig.taxa[which(sig.taxa$Value > 30000000), ]
# 
# # Must have at least 10 diel genes
# sig.taxa <- sig.taxa[which(sig.taxa$ShortTax == "Actinobacteria" | sig.taxa$ShortTax == "Bacteroidetes" | sig.taxa$ShortTax == "Bdellovibrionales" | sig.taxa$ShortTax == "Comamonadaceae" | sig.taxa$ShortTax == "Flavobacteriales" | sig.taxa$ShortTax == "Gemmatimonas"),]
# 
# ggplot(sig.taxa, aes(x = Time, y = ShortTax, color = ShortTax, size = log(Value))) + geom_point()
# 
# 
# # Skip back to DESeq
# abun_mnorm <- mnorm[order(rowSums(mnorm), decreasing = T), ]
# abun_mnorm <- abun_mnorm[1:20000,]
# 
# colnames(abun_mnorm) <- gsub(".nonrRNA", "", colnames(abun_mnorm))
# 
# input <- as.matrix(abun_mnorm)
# input <- input/100
# input <- round(input, digits = 0)
# 
# 
# conditions <- metadata$Time[match(colnames(abun_mnorm), metadata$Sample)]
# conditions[which(conditions == 9 | conditions == 13 | conditions == 17)] <- "day"
# conditions[which(conditions == 5 | conditions == 21 | conditions == 1)] <- "night"
# coldata <- data.frame(samples = colnames(abun_mnorm), conditions)
# 
# cds <- DESeqDataSetFromMatrix(countData = input,
#                               colData = coldata,
#                               design = ~ conditions)
# 
# cds <- estimateSizeFactors(cds) 
# cds <- estimateDispersions(cds) 
# 
# 
# dds <- DESeq(cds)
# res <- results(dds)
# sig.res <- res[which(res$padj < 0.05), ]
# sig.res.day <- sig.res[which(sig.res$log2FoldChange < 0), ]
# sig.res.night <- sig.res[which(sig.res$log2FoldChange > 0), ]
# sig.res.day.key <- mendota_key[match(rownames(sig.res.day), mendota_key$Gene),]
# sig.res.night.key <- mendota_key[match(rownames(sig.res.night), mendota_key$Gene),]
# 
# # Make plots of this aggreagted by short taxa and groupd products
# 
# sig.res.day.key$Taxonomy <- gsub(";;;", "", sig.res.day.key$Taxonomy)
# sig.res.day.key$Taxonomy <- gsub(";;", "", sig.res.day.key$Taxonomy)
# sig.res.day.key$Taxonomy <- gsub(";$", "", sig.res.day.key$Taxonomy)
# spl <- strsplit(as.character(sig.res.day.key$Taxonomy), ";")
# sig.res.day.key$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")
# 
# sig.res.night.key$Taxonomy <- gsub(";;;", "", sig.res.night.key$Taxonomy)
# sig.res.night.key$Taxonomy <- gsub(";;", "", sig.res.night.key$Taxonomy)
# sig.res.night.key$Taxonomy <- gsub(";$", "", sig.res.night.key$Taxonomy)
# spl <- strsplit(as.character(sig.res.night.key$Taxonomy), ";")
# sig.res.night.key$ShortTax <- sapply(lapply(spl, tail, 1), paste, collapse=";")
# 
# # Make column by certain key words 
# sig.res.day.key$Category <- "None"
# sig.res.day.key$Category[grep("photo|Photo", sig.res.day.key$Product)] <- "Photosynthesis"
# sig.res.day.key$Category[grep("rhodopsin|Rhodopsin", sig.res.day.key$Product)] <- "Rhodopsin"
# sig.res.day.key$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate", sig.res.day.key$Product)] <- "Sugar degradation"
# sig.res.day.key$Category[grep("rbcL|ribulose-bisphosphate carboxylase", sig.res.day.key$Product)] <- "RubisCO"
# sig.res.day.key$Category[grep("putrescine|Putrescine|spermidine|Spermidine", sig.res.day.key$Product)] <- "Polyamines"
# 
# sig.res.night.key$Category <- "None"
# sig.res.night.key$Category[grep("photo|Photo", sig.res.night.key$Product)] <- "Photosynthesis"
# sig.res.night.key$Category[grep("rhodopsin|Rhodopsin", sig.res.night.key$Product)] <- "Rhodopsin"
# sig.res.night.key$Category[grep("sugar|Sugar|ribose|Ribose|hexose|maltose|carbohydrate|Carbohydrate", sig.res.night.key$Product)] <- "Sugar degradation"
# sig.res.night.key$Category[grep("rbcL|ribulose-bisphosphate carboxylase", sig.res.night.key$Product)] <- "RubisCO"
# sig.res.night.key$Category[grep("putrescine|Putrescine|spermidine|Spermidine", sig.res.night.key$Product)] <- "Polyamines"
# 
# sig.res.day.key$Condition <- "day"
# sig.res.night.key$Condition <- "night"
# 
# sig.res.key <- rbind(sig.res.day.key, sig.res.night.key)
# 
# ggplot(data = sig.res.key, aes(x = Category, y = totals, group = ShortTax)) + geom_bar(stat = "identity")
# 
# ggplot(data = sig.res.key[which(sig.res.key$Category != "None"), ], aes(x = Category, y = totals, fill = Condition)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "Mendota")
# ggplot(data = sig.res.key[which(sig.res.key$Category != "None" & sig.res.key$Category != "Photosynthesis"), ], aes(x = Category, y = totals, fill = Condition)) + geom_bar(stat = "identity", position = "dodge") + labs(title = "Mendota")

# Aggregate by marker gene

reduced_searchterm <- "ammonia monooxygenase|nitrogenase|Nitrogenase|Nif|Chitobiase|chitobiase|chitinase|Chitinase|glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase|carbon monoxide dehydrogenase|Carbon monoxide dehydrogenase|carbon-monoxide dehydrogenase|Carbon-monoxide dehydrogenase|sulfite reductase|Sulfite reductase|glucosaminidase|hexosaminidase|glyoxal oxidase|Glyoxal oxidase|galactose oxidase|Galactose oxidase|laccase|Glutathione S-transferase|glutathione S-transferase|ligninase|lignin peroxidase|Lignin peroxidase|manganese peroxidase|Mn peroxidase|Manganese peroxidase|coenzyme M reductase|methane monooxygenase|Methane monooxygenase|nitrate reductase|Nitrate reductase|nitrite reductase|Nitrite reductase|nitric oxide reductase|Nitric oxide reductase|nitrous oxide reductase|Nitrous oxide reductase|nitrous-oxide reductase|Nitrous-oxide reductase|Bacillolysin|bacillolysin|Thermolysin|thermolysin|Fungalysin|fungalysin|nitrate oxidoreductase|Nitrate oxidoreductase|nitrite oxidoreductase|Nitrite oxidoreductase|lactaldehyde dehydrogenase|aldehyde dehydrogenase|phenoloxidase|Phenoloxidase|sox|Sox|sulfate thiol esterase|sulfur oxidation|Sulfur oxidation|sulfur-oxidizing|Sulfur-oxidizing|serine protease|urease|xylose isomerase|Xylose isomerase|citrate lyase|Citrate lyase|sulfate reductase|Sulfate reductase|sulfate kinase|adenylyltransferase|formaldehyde|Formaldehyde|sulfide dehydrogenase|Sulfide dehydrogenase|formate dehydrogenase|Formate dehydrogenase|formyltransferase|hydrazine|Hydrazine|methylamine|Methylamine|methanol|Methanol|DMSO reductase|Rubisco|RuBisCO|rubisco|bisphosphate carboxylase|sulfur dioxygenase|Sulfur dioxygenase|sulfide quinione reductase|Sulfide quinone reductase|opsin|rhamnulose-1-phosphate|Rhamnulose-1-phosphate|fuculose-phosphate|Fuculose phosphate|ribulokinase|mannose-6-phosphate|putrescine|spermidine|polymamine|Putrescine|Spermidine|Polyamine|transport|urea carboxylase|Urea carboxylase|sulfide-quinone reductase|Sulfide-quinione reductase|photo|chlorophyll|Photo|Chlorophyll|alkaline phosphatase|Alkaline phosphatase|cytochrome c|ATP synthase|Cytochrome c|Starch-binding|starch phosphorylase|glycogen/starch/alpha-glucan phosphorylases|Starch synthase|starch synthase|carboxylate transport"

keep <- grep(reduced_searchterm, mendota_key$Product)
abun_mnorm <- abun_mnorm[keep,]
mendota_key <- mendota_key[keep,]

keep <- grep(reduced_searchterm, trout_key$Product)
abun_tnorm <- abun_tnorm[keep,]
trout_key <- trout_key[keep,]

keep <- grep(reduced_searchterm, spark_key$Product)
abun_snorm <- abun_snorm[keep,]
spark_key <- spark_key[keep,]

abun_mnorm$Genes <- rownames(abun_mnorm)
abun_mnorm <- melt(abun_mnorm)
abun_mnorm$variable <- gsub(".nonrRNA", "", abun_mnorm$variable)
abun_mnorm$Timepoint <- metadata$Time[match(abun_mnorm$variable, metadata$Sample)]
abun_mnorm$Lake <- "Lake Mendota"
abun_mnorm$Product <- mendota_key$Product[match(abun_mnorm$Genes, mendota_key$Gene)]
abun_mnorm$Taxonomy <- mendota_key$Taxonomy[match(abun_mnorm$Genes, mendota_key$Gene)]

abun_tnorm$Genes <- rownames(abun_tnorm)
abun_tnorm <- melt(abun_tnorm)
abun_tnorm$variable <- gsub(".nonrRNA", "", abun_tnorm$variable)
abun_tnorm$Timepoint <- metadata$Time[match(abun_tnorm$variable, metadata$Sample)]
abun_tnorm$Lake <- "Trout Bog"
abun_tnorm$Product <- trout_key$Product[match(abun_tnorm$Genes, trout_key$Gene)]
abun_tnorm$Taxonomy <- trout_key$Taxonomy[match(abun_tnorm$Genes, trout_key$Gene)]

abun_snorm$Genes <- rownames(abun_snorm)
abun_snorm <- melt(abun_snorm)
abun_snorm$variable <- gsub(".nonrRNA", "", abun_snorm$variable)
abun_snorm$Timepoint <- metadata$Time[match(abun_snorm$variable, metadata$Sample)]
abun_snorm$Lake <- "Sparkling Lake"
abun_snorm$Product <- spark_key$Product[match(abun_snorm$Genes, spark_key$Gene)]
abun_snorm$Taxonomy <- spark_key$Taxonomy[match(abun_snorm$Genes, spark_key$Gene)]


#Combine all lakes
all_lakes <- rbind(abun_mnorm, abun_tnorm, abun_snorm)

all_lakes$Category <- "unknown"
all_lakes$Category[grep("ribulose-bisphosphate|Rubisco", all_lakes$Product)] <- "Carbon fixation"
all_lakes$Category[grep("sugar transport|carbohydrate ABC transporter|monosaccharide ABC transporter|mannose ABC transporter ATP-binding protein/fructose ABC transporter ATP-binding protein/ribose|Ribose/xylose/arabinose/galactoside ABC-type transport|Sugar (and other) transport|galactose transport", all_lakes$Product)] <- "Sugar transport"
all_lakes$Category[grep("rhodopsin|opsin", all_lakes$Product)] <- "Rhodopsin"
all_lakes$Category[grep("rhamnose transport|L-rhamnose-H+ transport", all_lakes$Product)] <- "Sugar transport"
all_lakes$Category[grep("spermidine|putrescine|Spermidine|Putrescine", all_lakes$Product)] <- "Polyamine metabolism"
all_lakes$Category[grep("ribose transport", all_lakes$Product)] <- "Sugar transport"
all_lakes$Category[grep("raffinose/stachyose/melibiose transport", all_lakes$Product)] <- "Sugar transport"
all_lakes$Category[grep("glucose/mannose transport", all_lakes$Product)] <- "Sugar transport"
all_lakes$Category[grep("xylose transport", all_lakes$Product)] <- "Sugar transport"
all_lakes$Category[grep("glucosaminidase", all_lakes$Product)] <- "Complex carbon degradation"
all_lakes$Category[grep("glucoside", all_lakes$Product)] <- "Complex carbon degradation"
all_lakes$Category[grep("Biopolymer transport|biopolymer transport", all_lakes$Product)] <- "Complex carbon degradation"
all_lakes$Category[grep("maltooligosaccharide transport|moltooligosaccharide transport", all_lakes$Product)] <- "Complex carbon degradation"
all_lakes$Category[grep("carotenoid", all_lakes$Product)] <- "Rhodopsin"
all_lakes$Category[grep("fructose transport", all_lakes$Product)] <- "Sugar transport"
all_lakes$Category[grep("C4-dicarboxylate transport", all_lakes$Product)] <- "Carboxylate transport"
all_lakes$Category[grep("chitobiose transport", all_lakes$Product)] <- "Complex carbon degradation"
all_lakes$Category[grep("nitrogenase|nitrogen fixation", all_lakes$Product)] <- "Nitrogen fixation"
all_lakes$Category[grep("chitinase|Chitinase", all_lakes$Product)] <- "Complex carbon degradation"
all_lakes$Category[grep("motility|gliding", all_lakes$Product)] <- "Motility"
all_lakes$Category[grep("nitrite reductase", all_lakes$Product)] <- "Nitrite/nitrate reduction"
all_lakes$Category[grep("nitrate reductase", all_lakes$Product)] <- "Nitrite/nitrate reduction"
all_lakes$Category[grep("tricarboxylate transport", all_lakes$Product)] <- "Carboxylate transport"
all_lakes$Category[grep("lactose/L-arabinose transport|arabinose ABC transport", all_lakes$Product)] <- "Sugar transport"
all_lakes$Category[grep("methylamine dehydrogenase|trimethylamine monooxygenase", all_lakes$Product)] <- "C1 metabolism"
all_lakes$Category[grep("nitrate/nitrite transport|Nitrate/nitrite transport|nitrite transport", all_lakes$Product)] <- "Nitrite/nitrate reduction"
all_lakes$Category[grep("methanol dehydrogenase|PQQ-dependent dehydrogenase", all_lakes$Product)] <- "C1 metabolism"
all_lakes$Category[grep("methane/ammonia monooxygenase|methane monooxygenase|ammonia monoxygenase", all_lakes$Product)] <- "C1 metabolism"
all_lakes$Category[grep("hexosaminidase", all_lakes$Product)] <- "Complex carbon degradation"
all_lakes$Category[grep("glycosyl hydrolase", all_lakes$Product)] <- "Complex carbon degradation"
all_lakes$Category[grep("citryl-CoA lyase", all_lakes$Product)] <- "Carbon fixation"
all_lakes$Category[grep("cellobiose transport", all_lakes$Product)] <- "Complex carbon degradation"
all_lakes$Category[grep("sulfur-oxidizing", all_lakes$Product)] <- "Sulfur oxidation"
all_lakes$Category[grep("photosystem II|Photosystem II", all_lakes$Product)] <- "Photosynthesis"
all_lakes$Category[grep("photosystem I |Photosystem I ", all_lakes$Product)] <- "Photosynthesis"
all_lakes$Category[grep("photosynthetic reaction centre|photosynthetic reaction center|Photosynthetic reaction center", all_lakes$Product)] <- "Photosynthesis"
all_lakes$Category[grep("chlorophyllide a reductase|chlorophyllide reductase", all_lakes$Product)] <- "Chlorophyll synthesis"
all_lakes$Category[grep("light-independent protochlorophyllide|ferredoxin protochlorophyllide|protochlorophyllide oxidoreductase", all_lakes$Product)] <- "Chlorophyll synthesis"
all_lakes$Category[grep("chlorophyll synthase", all_lakes$Product)] <- "Chlorophyll synthesis"
all_lakes$Category[grep("alkaline phosphatase|Alkaline phosphatase", all_lakes$Product)] <- "Alkaline phosphatase"
all_lakes$Category[grep("cytochrome c|ATP synthase|Cytochrome c", all_lakes$Product)] <- "Oxidative phosphorylation"
all_lakes$Category[grep("carboxylate transport", all_lakes$Product)] <- "Carboxylate transport"

all_lakes <- all_lakes[which(all_lakes$Category != "unknown"), ]
all_lakes$Condition <- "NA"
all_lakes$Condition[which(all_lakes$Timepoint == 9 | all_lakes$Timepoint == 13 | all_lakes$Timepoint == 17)] <- "day"
all_lakes$Condition[which(all_lakes$Timepoint == 5 | all_lakes$Timepoint == 21 | all_lakes$Timepoint == 1)] <- "night"


# Start testing! First by lake. Rerun for each pair of lakes

lake_comparison <- aggregate(value ~ variable + Category, all_lakes[which(all_lakes$Lake != "Sparkling Lake"),], sum)
lake_comparison <- reshape(lake_comparison, idvar = "Category", timevar = "variable", direction = "wide")
lake_comparison <- lake_comparison[which(is.na(lake_comparison[,2]) == F), ]
rownames(lake_comparison) <- lake_comparison[,1]
lake_comparison <- lake_comparison[,2:dim(lake_comparison)[2]]
colnames(lake_comparison) <- gsub("value.", "", colnames(lake_comparison))

input <- as.matrix(lake_comparison)
input <- input/500
input <- round(input, digits = 0)
input[which(is.na(input) == T)] <- 0


conditions <- metadata$Lake[match(colnames(lake_comparison), metadata$Sample)]
coldata <- data.frame(samples = colnames(lake_comparison), conditions)

cds <- DESeqDataSetFromMatrix(countData = input,
                              colData = coldata,
                              design = ~ conditions)

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds) 


dds <- DESeq(cds)
res <- results(dds)
sig.res <- res[which(res$padj < 0.05), ]
comp1 <- rownames(sig.res) #TB vs SP
comp2 <- rownames(sig.res) #ME vs SP
comp3 <- rownames(sig.res) #ME vs TB



plot.lake_comparison <- all_lakes[which(all_lakes$Category %in% unique(c(comp1, comp2, comp3))), ]
plot.lake_comparison <- aggregate(value ~ Category + Lake, plot.lake_comparison, mean)

# Photosystem I, photosystem II, and Rubisco all strongly and highly expressed in TB. Photosystem II most of all - more expressed than photosystem I in all lakes.

plot.lake_comparison$value[which(plot.lake_comparison$value > 100000000)] = 100000000
ggplot(data = plot.lake_comparison, aes(x = Category, y = value, fill = Lake)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = c("#5ab4ac", "#91bfdb", "#d8b365")) + labs(x = NULL, y = "Mean Transcripts/L") + coord_flip()


# Now in each lake, test day vs night.
metadata$Condition <- "NA"
metadata$Condition[which(metadata$Time == 9 | metadata$Time == 13 | metadata$Time == 17)] <- "day"
metadata$Condition[which(metadata$Time == 5 | metadata$Time == 21 | metadata$Time == 1)] <- "night"

time_comparison <- aggregate(value ~ variable + Category, all_lakes[which(all_lakes$Lake == "Lake Mendota"),], sum)
time_comparison <- reshape(time_comparison, idvar = "Category", timevar = "variable", direction = "wide")
time_comparison <- time_comparison[which(is.na(time_comparison[,2]) == F), ]
rownames(time_comparison) <- time_comparison[,1]
time_comparison <- time_comparison[,2:dim(time_comparison)[2]]
colnames(time_comparison) <- gsub("value.", "", colnames(time_comparison))

input <- as.matrix(time_comparison)
input <- input/500
input <- round(input, digits = 0)
input[which(is.na(input) == T)] <- 0


conditions <- metadata$Condition[match(colnames(time_comparison), metadata$Sample)]
coldata <- data.frame(samples = colnames(time_comparison), conditions)

cds <- DESeqDataSetFromMatrix(countData = input,
                              colData = coldata,
                              design = ~ conditions)

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds) 


dds <- DESeq(cds)
res <- results(dds)
sig.res <- res[which(res$padj < 0.05), ]

plot.time <- all_lakes[which(all_lakes$Category %in% c(rownames(sig.res))),]
plot.time <- plot.time[which(plot.time$Lake == "Lake Mendota"),]
plot.time <- aggregate(value ~ Category + Condition, plot.time, mean)

ggplot(data = plot.time, aes(x = Category, y = value, fill = Condition)) + geom_bar(stat = "identity", position = "dodge")  + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = c("darkgoldenrod2", "royalblue4")) + labs(x = NULL, y = "Mean Transcripts/L", title = "Lake Mendota") + coord_flip()


x <- aggregate(value ~ Taxonomy + Condition, all_lakes[which(all_lakes$Category == "Rhodopsin" & all_lakes$Lake == "Sparkling Lake"), ], mean)
ggplot(data = x, aes(x = Taxonomy, y = value, fill = Condition)) + geom_bar(stat = "identity", position = "dodge") + coord_flip()
