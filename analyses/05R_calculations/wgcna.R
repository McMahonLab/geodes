#Networks!

library(WGCNA)
allowWGCNAThreads()

snorm <- read.csv("D:/GEODES_mapping_summaries/Sparkling_normalized_counts_2017-06-12.csv", header = T, row.names = 1)

# Pre-processing
# Cluster samples to identify weird ones

sampletree <- hclust(dist(t(snorm)), method = "average")
plot(sampletree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
# One outlier for Sparkling - GEODES065. Remove this sample.
# Overall replicates seem to cluster well in the tree, if not perfectly. Seems like there are some very similar timepoints.

clust <- cutreeStatic(sampletree, cutHeight = 1500000, minSize = 10)
table(clust)
keepSamples <- (clust == 1)
snorm_clean <- snorm[, keepSamples]
nGenes <- nrow(snorm_clean)
nSamples <- ncol(snorm_clean)
snorm_clean <- t(snorm_clean)
snorm_clean <- snorm_clean[, which(colSums(snorm_clean) > quantile(colSums(snorm_clean), 0.5))]

options(stringsAsFactors = FALSE)
sft = pickSoftThreshold(snorm_clean, verbose = 3, blockSize = 100)
