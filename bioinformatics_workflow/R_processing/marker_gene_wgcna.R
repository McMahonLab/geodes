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
#snorm <- read.csv("D:/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
#tnorm <- read.csv("D:/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mnorm <- read.csv("D:/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

# Gene keys
mendota_key <- read.csv("D:/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#spark_key <- read.csv("D:/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#trout_key <- read.csv("D:/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Alex/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# Pull out genes encoding functional markers - annotations pulled from fungenes biogeochemical cycling section

searchterm <- "alkaline phosphatase|Alkaline phosphatase|ammonia monooxygenase|nitrogenase|Nitrogenase|Nif|phytase|Phytase|acetyl-CoA hydrolase|Acetyl-CoA hydrolase|4-hydroxybutyrate CoA-transferase|Chitobiase|chitobiase|chitinase|Chitinase|glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase|carbon monoxide dehydrogenase|Carbon monoxide dehydrogenase|carbon-monoxide dehydrogenase|Carbon-monoxide dehydrogenase|sulfite reductase|Sulfite reductase|glucosaminidase|hexosaminidase|cytochrome|Cytochrome|peptidase|Peptidase|glyoxal oxidase|Glyoxal oxidase|galactose oxidase|Galactose oxidase|Histidine phosphatase|histidine phosphatase|iron hydrogenase|Iron hydrogenase|ferredoxin|Ferredoxin|hydrogenase, Fe-only|Fe-only hydrogenase|laccase|Glutathione S-transferase|glutathione S-transferase|ligninase|lignin peroxidase|Lignin peroxidase|manganese peroxidase|Mn peroxidase|Manganese peroxidase|coenzyme M reductase|methane monooxygenase|Methane monooxygenase|nitrate reductase|Nitrate reductase|nitrite reductase|Nitrite reductase|nitric oxide reductase|Nitric oxide reductase|quinol oxidase|Quinol oxidase|nitrous oxide reductase|Nitrous oxide reductase|nitrous-oxide reductase|Nitrous-oxide reductase|Bacillolysin|bacillolysin|Thermolysin|thermolysin|Fungalysin|fungalysin|protease|nitrate oxidoreductase|Nitrate oxidoreductase|nitrite oxidoreductase|Nitrite oxidoreductase|serine/threonine protein phosphatase|Serine/threonine protein phosphatase|Ser/Thr protein phosphatase|phosphonoacetaldehyde dehydrogenase|phosphonoacetate hydrolase|glyceraldehyde-3-phosphate dehydrogenase|lactaldehyde dehydrogenase|aldehyde dehydrogenase|Phosphonoacetaldehyde dehydrogenase|Phosphonoacetate hydrolase|Glyceraldehyde-3-phosphate dehydrogenase|Lactaldehyde dehydrogenase|Aldehyde dehydrogenase|phosphodiesterase|Phosphodiesterase|phenoloxidase|Phenoloxidase|sox|Sox|sulfate thiol esterase|sulfur oxidation|Sulfur oxidation|sulfur-oxidizing|Sulfur-oxidizing|serine protease|urease|xylose isomerase|Xylose isomerase|citrate lyase|Citrate lyase|sulfate reductase|Sulfate reductase|sulfate kinase|adenylyltransferase|formaldehyde|Formaldehyde|sulfide dehydrogenase|Sulfide dehydrogenase|formate dehydrogenase|Formate dehydrogenase|formyltransferase|hydrazine|Hydrazine|methylamine|Methylamine|methanol|Methanol|DMSO reductase|propionyl-CoA|Propionyl-CoA|Rubisco|RuBisCO|rubisco|bisphosphate carboxylase|sulfur dioxygenase|Sulfur dioxygenase|sulfide quinione reductase|Sulfide quinone reductase|hydroxybutryl-CoA|opsin|alcohol dehydrogenase|Alcohol dehydrogenase|rhamnulose-1-phosphate|Rhamnulose-1-phosphate|fuculose-phosphate|Fuculose phosphate|ribulokinase|epimerase|mannose-6-phosphate|putrescine|spermidine|polymamine|Putrescine|Spermidine|Polyamine|transport|urea carboxylase|Urea carboxylase|sulfide-quinone reductase|Sulfide-quinione reductase"


marker_genes <- mendota_key[grep(searchterm, mendota_key$Product), ]
#unique(as.character(marker_genes$Product))


#Cuts 3 million plus genes to 300,000
# Light abundance filtering
abun_mnorm <- mnorm[match(marker_genes$Gene, rownames(mnorm)),]
abun_mnorm <- abun_mnorm[which(rowSums(abun_mnorm) > (dim(abun_mnorm)[2] * 20000)), ]

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
#Mendota Soft threshold = 7 in this case

mendota_net <- blockwiseModules(new_abun_mnorm, maxBlockSize = 1000, power = 7, loadTOM = F, saveTOMs = T, networkType = "signed", minModuleSize = 30, numericLabels = T, verbose = 3, threads = 4)

