# Make a network just of metabolic marker genes

# Because this consumes a lot of RAM, delete objects you no longer need as soon as you no longer need them
### Load packages
library(ggplot2)
library(cowplot)
library(reshape2)
library(WGCNA)
allowWGCNAThreads(nThreads = 4)

### Load data (start with only one to save RAM and comment the rest out)
# Normalized read tables
#snorm <- read.csv("E:/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
#tnorm <- read.csv("E:/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mnorm <- read.csv("E:/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

# Gene keys
mendota_key <- read.csv("E:/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#spark_key <- read.csv("E:/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#trout_key <- read.csv("E:/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# Pull out genes encoding functional markers - annotations pulled from fungenes biogeochemical cycling section

#searchterm <- "alkaline phosphatase|Alkaline phosphatase|ammonia monooxygenase|nitrogenase|Nitrogenase|Nif|phytase|Phytase|acetyl-CoA hydrolase|Acetyl-CoA hydrolase|4-hydroxybutyrate CoA-transferase|Chitobiase|chitobiase|chitinase|Chitinase|glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase|carbon monoxide dehydrogenase|Carbon monoxide dehydrogenase|carbon-monoxide dehydrogenase|Carbon-monoxide dehydrogenase|sulfite reductase|Sulfite reductase|glucosaminidase|hexosaminidase|cytochrome|Cytochrome|peptidase|Peptidase|glyoxal oxidase|Glyoxal oxidase|galactose oxidase|Galactose oxidase|Histidine phosphatase|histidine phosphatase|iron hydrogenase|Iron hydrogenase|ferredoxin|Ferredoxin|hydrogenase, Fe-only|Fe-only hydrogenase|laccase|Glutathione S-transferase|glutathione S-transferase|ligninase|lignin peroxidase|Lignin peroxidase|manganese peroxidase|Mn peroxidase|Manganese peroxidase|coenzyme M reductase|methane monooxygenase|Methane monooxygenase|nitrate reductase|Nitrate reductase|nitrite reductase|Nitrite reductase|nitric oxide reductase|Nitric oxide reductase|quinol oxidase|Quinol oxidase|nitrous oxide reductase|Nitrous oxide reductase|nitrous-oxide reductase|Nitrous-oxide reductase|Bacillolysin|bacillolysin|Thermolysin|thermolysin|Fungalysin|fungalysin|protease|nitrate oxidoreductase|Nitrate oxidoreductase|nitrite oxidoreductase|Nitrite oxidoreductase|serine/threonine protein phosphatase|Serine/threonine protein phosphatase|Ser/Thr protein phosphatase|phosphonoacetaldehyde dehydrogenase|phosphonoacetate hydrolase|glyceraldehyde-3-phosphate dehydrogenase|lactaldehyde dehydrogenase|aldehyde dehydrogenase|Phosphonoacetaldehyde dehydrogenase|Phosphonoacetate hydrolase|Glyceraldehyde-3-phosphate dehydrogenase|Lactaldehyde dehydrogenase|Aldehyde dehydrogenase|phosphodiesterase|Phosphodiesterase|phenoloxidase|Phenoloxidase|sox|Sox|sulfate thiol esterase|sulfur oxidation|Sulfur oxidation|sulfur-oxidizing|Sulfur-oxidizing|serine protease|urease|xylose isomerase|Xylose isomerase|citrate lyase|Citrate lyase|sulfate reductase|Sulfate reductase|sulfate kinase|adenylyltransferase|formaldehyde|Formaldehyde|sulfide dehydrogenase|Sulfide dehydrogenase|formate dehydrogenase|Formate dehydrogenase|formyltransferase|hydrazine|Hydrazine|methylamine|Methylamine|methanol|Methanol|DMSO reductase|propionyl-CoA|Propionyl-CoA|Rubisco|RuBisCO|rubisco|bisphosphate carboxylase|sulfur dioxygenase|Sulfur dioxygenase|sulfide quinione reductase|Sulfide quinone reductase|hydroxybutryl-CoA|opsin|alcohol dehydrogenase|Alcohol dehydrogenase|rhamnulose-1-phosphate|Rhamnulose-1-phosphate|fuculose-phosphate|Fuculose phosphate|ribulokinase|epimerase|mannose-6-phosphate|putrescine|spermidine|polymamine|Putrescine|Spermidine|Polyamine|transport|urea carboxylase|Urea carboxylase|sulfide-quinone reductase|Sulfide-quinione reductase"

reduced_searchterm <- "ammonia monooxygenase|nitrogenase|Nitrogenase|Nif|Chitobiase|chitobiase|chitinase|Chitinase|glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase|carbon monoxide dehydrogenase|Carbon monoxide dehydrogenase|carbon-monoxide dehydrogenase|Carbon-monoxide dehydrogenase|sulfite reductase|Sulfite reductase|glucosaminidase|hexosaminidase|glyoxal oxidase|Glyoxal oxidase|galactose oxidase|Galactose oxidase|laccase|Glutathione S-transferase|glutathione S-transferase|ligninase|lignin peroxidase|Lignin peroxidase|manganese peroxidase|Mn peroxidase|Manganese peroxidase|coenzyme M reductase|methane monooxygenase|Methane monooxygenase|nitrate reductase|Nitrate reductase|nitrite reductase|Nitrite reductase|nitric oxide reductase|Nitric oxide reductase|nitrous oxide reductase|Nitrous oxide reductase|nitrous-oxide reductase|Nitrous-oxide reductase|Bacillolysin|bacillolysin|Thermolysin|thermolysin|Fungalysin|fungalysin|nitrate oxidoreductase|Nitrate oxidoreductase|nitrite oxidoreductase|Nitrite oxidoreductase|lactaldehyde dehydrogenase|aldehyde dehydrogenase|phenoloxidase|Phenoloxidase|sox|Sox|sulfate thiol esterase|sulfur oxidation|Sulfur oxidation|sulfur-oxidizing|Sulfur-oxidizing|serine protease|urease|xylose isomerase|Xylose isomerase|citrate lyase|Citrate lyase|sulfate reductase|Sulfate reductase|sulfate kinase|adenylyltransferase|formaldehyde|Formaldehyde|sulfide dehydrogenase|Sulfide dehydrogenase|formate dehydrogenase|Formate dehydrogenase|formyltransferase|hydrazine|Hydrazine|methylamine|Methylamine|methanol|Methanol|DMSO reductase|Rubisco|RuBisCO|rubisco|bisphosphate carboxylase|sulfur dioxygenase|Sulfur dioxygenase|sulfide quinione reductase|Sulfide quinone reductase|opsin|rhamnulose-1-phosphate|Rhamnulose-1-phosphate|fuculose-phosphate|Fuculose phosphate|ribulokinase|mannose-6-phosphate|putrescine|spermidine|polymamine|Putrescine|Spermidine|Polyamine|transport|urea carboxylase|Urea carboxylase|sulfide-quinone reductase|Sulfide-quinione reductase"

#marker_genes <- mendota_key[grep(searchterm, mendota_key$Product), ] #309595 genes
marker_genes <- mendota_key[grep(reduced_searchterm, mendota_key$Product), ] #181146 genes
#unique(as.character(marker_genes$Product))


# Light abundance filtering
abun_mnorm <- mnorm[match(marker_genes$Gene, rownames(mnorm)),]
abun_mnorm <- abun_mnorm[which(rowSums(abun_mnorm) > (dim(abun_mnorm)[2] * 100000)), ]
# Down to 143080 genes at 5000
# 117851 at 10000
# 92128 at 20000. That should run just fine on my work desktop.
# 60056 at 50000
# 40452 at 100000


#Aggregate by timepoint
rm(mnorm)

abun_mnorm$Genes <- rownames(abun_mnorm)
abun_mnorm <- melt(abun_mnorm)
abun_mnorm$variable <- gsub(".nonrRNA", "", abun_mnorm$variable)
abun_mnorm$Timepoint <- metadata$Timepoint[match(abun_mnorm$variable, metadata$Sample)]
agg_abun_mnorm <- aggregate(value ~ Genes + Timepoint, data = abun_mnorm, mean)
new_abun_mnorm <- reshape(agg_abun_mnorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(new_abun_mnorm) <- new_abun_mnorm[, 1]
new_abun_mnorm <- new_abun_mnorm[, 2:dim(new_abun_mnorm)[2]]
new_abun_mnorm <- t(new_abun_mnorm)

rm(agg_abun_mnorm)

# Check the soft threshold - only need to do this when running a subset of data for the 1st time

mpowers <- c()
for(i in 1:100){
  lottery <- sample(1:dim(new_abun_mnorm)[2], 2000, replace = F)
  sft = pickSoftThreshold(new_abun_mnorm[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
  mpowers[i] <- sft$powerEstimate
}
hist(mpowers)
mean(mpowers)
median(mpowers)
#Mendota Soft threshold = 6 in this case

mendota_net <- blockwiseModules(new_abun_mnorm, maxBlockSize = 1000, power = 6, loadTOM = F, saveTOMs = T, networkType = "signed", minModuleSize = 30, numericLabels = T, verbose = 3, threads = 4)

sig.mendota.key <- mendota_key[match(colnames(new_abun_mnorm), mendota_key$Gene), ]
sig.mendota.key$Cluster <- mendota_net$colors
sig.mendota.key$Totals <- colSums(new_abun_mnorm)[match(sig.mendota.key$Gene, colnames(new_abun_mnorm))]

write.csv(sig.mendota.key, "C:/Users/Goose and Gander/Documents/WGCNA_mendota_results_nondiel_2018-05-20.csv", row.names = F)
write.csv(mendota_net$MEs, "C:/Users/Goose and Gander/Documents/WGCNA_mendota_eigenvectors_nondiel_2018-05-20.csv", row.names = T)

rm(mendota_key)
rm(abun_mnorm)
rm(mendota_net)
rm(new_abun_mnorm)
rm(sig.mendota.key)
#Repeat for other lakes

### Load data (start with only one to save RAM and comment the rest out)
# Normalized read tables
snorm <- read.csv("E:/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
#tnorm <- read.csv("E:/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
#mnorm <- read.csv("E:/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

# Gene keys
#mendota_key <- read.csv("E:/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
spark_key <- read.csv("E:/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#trout_key <- read.csv("E:/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

marker_genes <- spark_key[grep(reduced_searchterm, spark_key$Product), ] #146155 genes
#unique(as.character(marker_genes$Product))


# Light abundance filtering
abun_snorm <- snorm[match(marker_genes$Gene, rownames(snorm)),]
abun_snorm <- abun_snorm[which(rowSums(abun_snorm) > (dim(abun_snorm)[2] * 50000)), ]
# 28864 genes at 100000
# 44066 genes at 50000


#Aggregate by timepoint
rm(snorm)

abun_snorm$Genes <- rownames(abun_snorm)
abun_snorm <- melt(abun_snorm)
abun_snorm$variable <- gsub(".nonrRNA", "", abun_snorm$variable)
abun_snorm$Timepoint <- metadata$Timepoint[match(abun_snorm$variable, metadata$Sample)]
agg_abun_snorm <- aggregate(value ~ Genes + Timepoint, data = abun_snorm, mean)
new_abun_snorm <- reshape(agg_abun_snorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(new_abun_snorm) <- new_abun_snorm[, 1]
new_abun_snorm <- new_abun_snorm[, 2:dim(new_abun_snorm)[2]]
new_abun_snorm <- t(new_abun_snorm)

rm(agg_abun_snorm)

# Check the soft threshold - only need to do this when running a subset of data for the 1st time

spowers <- c()
for(i in 1:100){
  lottery <- sample(1:dim(new_abun_snorm)[2], 2000, replace = F)
  sft = pickSoftThreshold(new_abun_snorm[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
  spowers[i] <- sft$powerEstimate
}
hist(spowers)
mean(spowers)
median(spowers)
#spark Soft threshold = 4 in this case

spark_net <- blockwiseModules(new_abun_snorm, maxBlockSize = 1000, power = 4, loadTOM = F, saveTOMs = F, networkType = "signed", minModuleSize = 30, numericLabels = T, verbose = 3, threads = 4)

sig.spark.key <- spark_key[match(colnames(new_abun_snorm), spark_key$Gene), ]
sig.spark.key$Cluster <- spark_net$colors
sig.spark.key$Totals <- colSums(new_abun_snorm)[match(sig.spark.key$Gene, colnames(new_abun_snorm))]

write.csv(sig.spark.key, "C:/Users/Goose and Gander/Documents/WGCNA_spark_results_nondiel_2018-05-20.csv", row.names = F)
write.csv(spark_net$MEs, "C:/Users/Goose and Gander/Documents/WGCNA_spark_eigenvectors_nondiel_2018-05-20.csv", row.names = T)

rm(spark_key)
rm(abun_snorm)
rm(spark_net)
rm(new_abun_snorm)
rm(sig.spark.key)

### Lastly, Trout
### Load data (start with only one to save RAM and comment the rest out)
# Normalized read tables
#snorm <- read.csv("E:/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
tnorm <- read.csv("E:/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
#mnorm <- read.csv("E:/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

# Gene keys
#mendota_key <- read.csv("E:/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#spark_key <- read.csv("E:/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
trout_key <- read.csv("E:/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

marker_genes <- trout_key[grep(reduced_searchterm, trout_key$Product), ] #53212 genes
#unique(as.character(marker_genes$Product))


# Light abundance filtering
abun_tnorm <- tnorm[match(marker_genes$Gene, rownames(tnorm)),]
abun_tnorm <- abun_tnorm[which(rowSums(abun_tnorm) > (dim(abun_tnorm)[2] * 10000)), ]
# 47697 at 100000


#Aggregate by timepoint
rm(tnorm)

abun_tnorm$Genes <- rownames(abun_tnorm)
abun_tnorm <- melt(abun_tnorm)
abun_tnorm$variable <- gsub(".nonrRNA", "", abun_tnorm$variable)
abun_tnorm$Timepoint <- metadata$Timepoint[match(abun_tnorm$variable, metadata$Sample)]
agg_abun_tnorm <- aggregate(value ~ Genes + Timepoint, data = abun_tnorm, mean)
new_abun_tnorm <- reshape(agg_abun_tnorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(new_abun_tnorm) <- new_abun_tnorm[, 1]
new_abun_tnorm <- new_abun_tnorm[, 2:dim(new_abun_tnorm)[2]]
new_abun_tnorm <- t(new_abun_tnorm)

rm(agg_abun_tnorm)

# Check the soft threshold - only need to do this when running a subset of data for the 1st time

tpowers <- c()
for(i in 1:100){
  lottery <- sample(1:dim(new_abun_tnorm)[2], 2000, replace = F)
  sft = pickSoftThreshold(new_abun_tnorm[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
  tpowers[i] <- sft$powerEstimate
}
hist(tpowers)
mean(tpowers)
median(tpowers)
#trout Soft threshold = 9 in this case

trout_net <- blockwiseModules(new_abun_tnorm, maxBlockSize = 1000, power = 9, loadTOM = F, saveTOMs = F, networkType = "signed", minModuleSize = 30, numericLabels = T, verbose = 3, threads = 4)

sig.trout.key <- trout_key[match(colnames(new_abun_tnorm), trout_key$Gene), ]
sig.trout.key$Cluster <- trout_net$colors
sig.trout.key$Totals <- colSums(new_abun_tnorm)[match(sig.trout.key$Gene, colnames(new_abun_tnorm))]

write.csv(sig.trout.key, "C:/Users/Goose and Gander/Documents/WGCNA_trout_results_nondiel_2018-05-20.csv", row.names = F)
write.csv(trout_net$MEs, "C:/Users/Goose and Gander/Documents/WGCNA_trout_eigenvectors_nondiel_2018-05-20.csv", row.names = T)

rm(trout_key)
rm(abun_tnorm)
rm(trout_net)
rm(new_abun_tnorm)
rm(sig.trout.key)
