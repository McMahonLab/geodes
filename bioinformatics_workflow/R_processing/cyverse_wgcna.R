# Cyverse network script

# Before starting R:
# mkdir geodes
# cd geodes
# iget Mendota_ID90_normalized_readcounts.csv
# iget Mendota_ID90_genekey_reclassified_2018-03-03.csv
# iget sample_metadata.csv

# Start R in sudo - password is the same as my Cyverse account

source("https://bioconductor.org/biocLite.R")
biocLite("GO.db")

chooseCRANmirror()
install.packages("cowplot")
install.packages("WGCNA")

library(ggplot2)
library(cowplot)
library(reshape2)
library(WGCNA)

mnorm <- read.csv("Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mendota_key <- read.csv("Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
metadata <- read.csv(file = "sample_metadata.csv", header = T)

abun_mnorm <- mnorm[which(rowSums(mnorm) > (dim(mnorm)[2] * 1000)), ]
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


####Calculate the mean of standard error for each timepoint for each gene

std_err <- function(x) (sd(x)/sqrt(length(x)))/mean(x) * 100

err_mnorm <- aggregate(value ~ Genes + Timepoint, data = abun_mnorm, std_err)
err_mnorm <- err_mnorm[which(is.na(err_mnorm$value) == F),]
err_mnorm <- aggregate(value ~ Genes, data = err_mnorm, mean)

keep <- err_mnorm$Genes[which(err_mnorm$value < 50)]

new_abun_mnorm <- new_abun_mnorm[, which(colnames(new_abun_mnorm) %in% keep)]

### WGCNA
# 1st, check the soft threshold in case removing samples changed it.

# mpowers <- c()
# for(i in 1:100){
#   lottery <- sample(1:dim(new_abun_mnorm)[2], 2000, replace = F)
#   sft = pickSoftThreshold(new_abun_mnorm[, lottery], verbose = 0, blockSize = 100, RsquaredCut = 0.7)
#   mpowers[i] <- sft$powerEstimate
# }
# hist(mpowers)
# mean(mpowers)
# median(mpowers)
# # Power of 7

allowWGCNAThreads(nThreads = 2)

mendota_net <- blockwiseModules(new_abun_mnorm, maxBlockSize = 1000, power = 7, loadTOM = F, saveTOMs = T, saveTOMFileBase = "mendota.network", networkType = "signed", minModuleSize = 30, numericLabels = T, verbose = 3, threads = 8)
