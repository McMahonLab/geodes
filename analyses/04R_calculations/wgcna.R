#Networks!

library(WGCNA)
allowWGCNAThreads()

snorm <- read.csv("D:/GEODES_mapping_summaries/Sparkling_normalized_counts_2017-06-12.csv", header = T, row.names = 1)

tnorm <- read.csv(file.choose(), header = T, row.names = 1)
metadata <- read.csv(file = "C:/Users/amlin/Desktop/geodes/analyses/04R_calculations/sample_metadata.csv", header = T, row.names = 1)
metadata$condition <- paste(metadata$Lake, metadata$Timepoint, sep = ";")
trout_meta <- metadata[which(metadata$Lake == "Trout"), ]

keep <- match(colnames(tnorm), rownames(trout_meta))
keep <- keep[which(is.na(keep) == F)]
tnorm <- tnorm[, keep]
tnorm <- tnorm[which(rowSums(tnorm) > 0), ]


# Pre-processing
# Cluster samples to identify weird ones

sampletree <- hclust(dist(t(tnorm)), method = "average")
plot(sampletree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
# One outlier for Sparkling - GEODES065. Remove this sample.
# Overall replicates seem to cluster well in the tree, if not perfectly. Seems like there are some very similar timepoints.

#clust <- cutreeStatic(sampletree, cutHeight = 1500000, minSize = 10)
#table(clust)
# keepSamples <- (clust == 1)
# snorm_clean <- snorm[, keepSamples]
nGenes <- nrow(tnorm)
nSamples <- ncol(tnorm)
tnorm_clean <- t(tnorm)
#tnorm_clean <- tnorm_clean[, which(colSums(snorm_clean) > quantile(colSums(snorm_clean), 0.5))]

options(stringsAsFactors = FALSE)
increments <- seq(1, dim(tnorm)[1], 5000)
tpowers <- c()
for(i in 2:length(increments)){
  sft = pickSoftThreshold(tnorm_clean[, increments[i-1]:increments[i]], verbose = 3, blockSize = 100, RsquaredCut = 0.7)
  tpowers[i-1] <- sft$powerEstimate
}

