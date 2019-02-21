
# Set up environment
path <- "/Users/Alex/Desktop/"
path2 <- "/Users/Alex/"

library(ggplot2)
library(cowplot)
library(reshape2)

zscore <- function(counts){
  z <- (counts - mean(counts)) / sd(counts)
  return(z)
}

# Sample data
metadata <- read.csv(file = paste(path2, "Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", sep = ""), header = T)
enviro_data <- read.csv(paste(path2, "Desktop/geodes/environmental_data/for_humans/compiled_field_data_for_R.csv", sep = ""), header = T)
enviro_data_no_sonde <- enviro_data[which(is.na(enviro_data$pH) == F), ]
chemdata <- read.csv(paste(path2, "Desktop/geodes/environmental_data/GEODES_TNTP.csv", sep = ""), header = T)
chemdata$Lake <- NA
chemdata$Lake[grep("ME", chemdata$Sample)] <- "Mendota"
chemdata$Lake[grep("SP", chemdata$Sample)] <- "Sparkling"
chemdata$Lake[grep("TB", chemdata$Sample)] <- "Trout"


mnorm <- read.csv(paste(path, "geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", sep = ""), header = T, row.names = 1)
mnorm <- mnorm[, which(colnames(mnorm) != "GEODES158.nonrRNA")]
colnames(mnorm) <- gsub(".nonrRNA", "", colnames(mnorm))
mendota_key <- read.csv(paste(path, "geodes_data_tables/Mendota_ID90_genekey_geneclassifications_2018-11-28.csv", sep = ""), header = T)

# Write a function to calculate CoV within replicates and average by gene
within_reps <- function(readcounts){
  readcounts2 <- readcounts[which(rowSums(readcounts) > 1000), ]
  colnames(readcounts2) <- gsub(".nonrRNA", "", colnames(readcounts2))
  timepoints <- metadata$Timepoint[match(colnames(readcounts2), metadata$Sample)]
  timelist <- unique(timepoints)
  time_table <- readcounts2[,which(timepoints == timelist[1])]
  results <- data.frame(timepoint1 = apply(time_table, 1, cv))
  
  for(i in 2:length(timelist)){
    if(length(which(timepoints == timelist[i])) > 1){
      time_table <- readcounts2[,which(timepoints == timelist[i])]
      next_timepoint <- data.frame(apply(time_table, 1, cv))
      colnames(next_timepoint) <- paste("timepoint", i, sep = "")
      results <- cbind(results, next_timepoint)
    }
  }
  return(results)
}

# mendota_within <- within_reps(mnorm)
# mendota_within <- apply(mendota_within, 1, mean, na.rm = T)
# mendota_across <- apply(mnorm, 1, cv, na.rm = T)
# cv_stat <- c(rep("within replicates", length(mendota_within)), rep("across samples", length(mendota_across)))
# mendota_cv_stats <- data.frame(type = cv_stat, stat = c(mendota_within, mendota_across))

# Save data to avoid rerunning the above code
# write.csv(mendota_cv_stats, file = paste(path2, "Desktop/intermediate_plotting_files/mendota_cv_stats.csv", sep = ""), quote = F, row.names = F)

# Read data back in
mendota_cv_stats <- read.csv(file = paste(path2, "Desktop/intermediate_plotting_files/mendota_cv_stats.csv", sep = ""), header = T)

p1 <- ggplot(mendota_cv_stats, aes(stat, fill = type)) + geom_density(alpha = 0.5) + scale_fill_manual(values = c("#d8b365", "#5ab4ac")) + labs(x = "Coefficent of Variance (%)", y = "Density",  title = "Lake Mendota") + theme(legend.title = element_blank())

# spark_within <- within_reps(snorm)
# spark_within <- apply(spark_within, 1, mean, na.rm = T)
# spark_across <- apply(snorm, 1, cv, na.rm = T)
# cv_stat <- c(rep("within replicates", length(spark_within)), rep("across samples", length(spark_across)))
# spark_cv_stats <- data.frame(type = cv_stat, stat = c(spark_within, spark_across))
# 
# p2 <- ggplot(spark_cv_stats, aes(stat, fill = type)) + geom_density(alpha = 0.5) + scale_fill_manual(values = c("#d8b365", "#5ab4ac")) + labs(x = "Coefficent of Variance (%)", title = "Sparkling Lake") + theme(legend.title = element_blank())
# 
# 
# trout_within <- within_reps(tnorm)
# trout_within <- apply(trout_within, 1, mean, na.rm = T)
# trout_across <- apply(tnorm, 1, cv, na.rm = T)
# cv_stat <- c(rep("within replicates", length(trout_within)), rep("across samples", length(trout_across)))
# trout_cv_stats <- data.frame(type = cv_stat, stat = c(trout_within, trout_across))
# 
# p3 <- ggplot(trout_cv_stats, aes(stat, fill = type)) + geom_density(alpha = 0.5) + scale_fill_manual(values = c("#d8b365", "#5ab4ac")) + labs(x = "Coefficent of Variance (%)", title = "Trout Bog") + theme(legend.title = element_blank())
# 
# plot_grid(p1, p2, p3, nrow = 3, labels = c("A", "B", "C"))

# Plot example trends
mnorm2 <- mnorm[1:1000,]
mnorm2$Genes <- rownames(mnorm2)
mnorm2 <- melt(mnorm2)
mnorm2$variable <- gsub(".nonrRNA", "", mnorm2$variable)
mnorm2$Time <- metadata$Timepoint[match(mnorm2$variable, metadata$Sample)]
agg_mnorm <- aggregate(value ~ Time + Genes, mnorm2, mean)
agg_mnorm$Time <- factor(agg_mnorm$Time, levels = c("0", "4", "8", "12", "16", "20", "24", "28", "32", "36", "40", "44"))

genes <- unique(agg_mnorm$Genes)
example <- data.frame("time" = agg_mnorm$Time[which(agg_mnorm$Genes == genes[1])], "eigenvector" = zscore(agg_mnorm$value[which(agg_mnorm$Genes == genes[32])]))

# Gene 6 is up down
# Gene 7 is on once only
# Gene 32 is day one


p2 <- ggplot(data = example, aes(x = time, y = eigenvector)) + geom_bar(stat = "identity") + labs(x = " ", y = " ")
example <- data.frame("time" = agg_mnorm$Time[which(agg_mnorm$Genes == genes[1])], "eigenvector" = zscore(agg_mnorm$value[which(agg_mnorm$Genes == genes[6])]))
p3 <- ggplot(data = example, aes(x = time, y = eigenvector)) + geom_bar(stat = "identity") + labs(x = " ", y = "Eigenvector of gene expression")
example <- data.frame("time" = agg_mnorm$Time[which(agg_mnorm$Genes == genes[1])], "eigenvector" = zscore(agg_mnorm$value[which(agg_mnorm$Genes == genes[7])]))
p4 <- ggplot(data = example, aes(x = time, y = eigenvector)) + geom_bar(stat = "identity") + labs(x = "Hours since start of time series", y = " ")

part2 <- plot_grid(p2, p3, p4, nrow = 3)
x <- plot_grid(p1, part2, nrow = 1, labels = c("A", "B"))

save_plot("/Users/Alex/Desktop/geodes/Manuscript/supplemental_materials/CoV.pdf", x, base_height = 6, base_aspect_ratio = 1.75)
