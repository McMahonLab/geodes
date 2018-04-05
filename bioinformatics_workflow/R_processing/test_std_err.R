# What threshold of error is needed to distinguish replicates from not replicates?
library(reshape2)

mnorm <- read.csv("D:/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)
metadata <- read.csv(file = "C:/Users/Alex/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)
abun_mnorm <- mnorm[which(rowSums(mnorm) > (dim(mnorm)[2] * 5000)), ]
rm(mnorm)

abun_mnorm$Genes <- rownames(abun_mnorm)
abun_mnorm <- melt(abun_mnorm)
abun_mnorm$variable <- gsub(".nonrRNA", "", abun_mnorm$variable)
abun_mnorm$Timepoint <- metadata$Timepoint[match(abun_mnorm$variable, metadata$Sample)]

std_err <- function(x) (sd(x)/sqrt(length(x)))/mean(x) * 100

genes2test <- unique(abun_mnorm$Genes)[1:1000]
withinreps <- c()
acrossreps <- c()
timepoints <- unique(metadata$Timepoint)
for(i in 1:length(genes2test)){
  set <- abun_mnorm[which(abun_mnorm$Genes == genes2test[i]),]
  means <- c()
  for(j in 1:length(timepoints)){
    reads <- set$value[which(set$Timepoint == timepoints[j])]
    means <- append(means, mean(reads), length(means))
    withinreps <- append(withinreps, std_err(reads), length(withinreps))
  }
  acrossreps <- append(acrossreps, std_err(means[which(is.na(means) == F)]), length(acrossreps))
}

#Conclusion: at least less than 40% error needed to distinguish replicates