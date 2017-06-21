#Networks!

library(WGCNA)
allowWGCNAThreads()

snorm <- read.csv("D:/GEODES_mapping_summaries/Sparkling_normalized_counts_2017-06-20.csv", header = T, row.names = 1)
tnorm <- read.csv("D:/GEODES_mapping_summaries/TroutBog_normalized_counts_2017-06-20.csv", header = T, row.names = 1)
mnorm <- read.csv("D:/GEODES_mapping_summaries/Mendota_normalized_counts_2017-06-20.csv", header = T, row.names = 1)

# metadata <- read.csv(file = "C:/Users/Alex/Desktop/geodes/analyses/04R_calculations/sample_metadata.csv", header = T, row.names = 1)
# metadata$condition <- paste(metadata$Lake, metadata$Timepoint, sep = ";")
# trout_meta <- metadata[which(metadata$Lake == "Trout"), ]
# 
# keep <- match(colnames(tnorm), rownames(trout_meta))
# keep <- keep[which(is.na(keep) == F)]
# tnorm <- tnorm[, keep]
# tnorm <- tnorm[which(rowSums(tnorm) > 0), ]

# Keep only reasonably abundant genes - I'm defining this as the top 25% of genes. This is roughly 25 normalized transcript counts for Sparkling, 17 for Trout, and 12 for Mendota. Sidenote, this seems to indicate that evenness is greatest in Mendota and lowest in Sparkling.

snorm <- snorm[which(rowSums(snorm) > quantile(rowSums(snorm), .75)), ]
tnorm <- tnorm[which(rowSums(tnorm) > quantile(rowSums(tnorm), .75)), ]
mnorm <- mnorm[which(rowSums(mnorm) > quantile(rowSums(mnorm), .75)), ]



# Pre-processing
# Cluster samples to identify weird ones

#Sparkling
sampletree <- hclust(dist(t(snorm)), method = "average")
plot(sampletree, main = "Sparkling", sub = "", xlab = "")
# No outliers

#Trout
sampletree <- hclust(dist(t(tnorm)), method = "average")
plot(sampletree, main = "Trout", sub = "", xlab = "")
# GEODES065 is a troublemaker. Remove.
clust <- cutreeStatic(sampletree, cutHeight = 1500000, minSize = 10)
table(clust)
keepSamples <- (clust == 1)
tnorm <- tnorm[, keepSamples]

#Mendota
sampletree <- hclust(dist(t(mnorm)), method = "average")
plot(sampletree, main = "Mendota", sub = "", xlab = "")
# GEODES158 is a little funky, but not off on its own. It can stay for now because it would be complicated to remove.

tnorm <- t(tnorm)
mnorm <- t(mnorm)
snorm <- t(snorm)

#Pick a soft threshold, using subsetting to avoid running the entire dataset. Start with Trout.
options(stringsAsFactors = FALSE)

tpowers <- c()
for(i in 1:100){
  lottery <- sample(1:dim(tnorm)[2], 2000, replace = F)
  sft = pickSoftThreshold(tnorm[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
  tpowers[i] <- sft$powerEstimate
}
hist(tpowers)
mean(tpowers)
median(tpowers)
# 3 is the winner!

#Sparkling
spowers <- c()
for(i in 1:100){
  lottery <- sample(1:dim(snorm)[2], 2000, replace = F)
  sft = pickSoftThreshold(snorm[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
  spowers[i] <- sft$powerEstimate
}
hist(spowers)
mean(spowers)
median(spowers)

# Also 3.

#Mendota
mpowers <- c()
for(i in 1:100){
  lottery <- sample(1:dim(mnorm)[2], 2000, replace = F)
  sft = pickSoftThreshold(mnorm[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
  mpowers[i] <- sft$powerEstimate
}
hist(mpowers)
mean(mpowers)
median(mpowers)

# and 3 again. Nice!

# Now make the network, starting with Trout. Complicated parameters so here's an explanation:
# maxBlockSize - the most genes that can be handled in one block. My desktop has 16GB RAM and should be able to handle 20000
# power is my soft threshold
# Leaving preclustering number at the default, based on number of genes (1613 for Trout). Can increase this to improve accuracy but slows down the analysis.
# TOM is a previously calcuated topological overlap matrix. I have not calculated this already so it is false. Nor do I want to save a TOM right now.
# leaving corType as Pearson for now - other option is bidweight midcorrelation
# networkType is signed, meaning use pos and neg correlations instead of absolute values
# minModuleSize - how many genes must be in a cluster for it to be considered a cluster? I'll say 10 arbitrarily
# Can switch merge cut height, but not sure what a good value is so I'll leave at 0.15
# Name clusters by number instead of color


trout_net <- blockwiseModules(tnorm, maxBlockSize = 27000, power = 3, loadTOM = F, saveTOMs = F, networkType = "signed", minModuleSize = 10, numericLabels = T, nThreads = 8, verbose = 3)



