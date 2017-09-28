# Load packages
library(ggplot2)
library(reshape2)
library(cowplot)
library(OTUtable)

# Load data
metadata <- read.csv(file = "C:/Users/Alex/Desktop/geodes/analyses/04R_calculations/sample_metadata.csv", header = T, row.names = 1)
genekey <- read.table(file = "D:/GEODES_mapping_summaries/GEODES_genekey_2017-06-26.txt", header = F, sep = "", quote = "", fill = T, colClasses = c("character"))
MAG_data <- read.csv(file = "C:/Users/Alex/Desktop/geodes/analyses/03process_mapping_results/Readme.csv", header = T, row.names = 1)

snorm <- read.csv("D:/GEODES_mapping_summaries/Sparkling_normalized_counts_2017-06-20.csv", header = T, row.names = 1)
tnorm <- read.csv("D:/GEODES_mapping_summaries/TroutBog_normalized_counts_2017-06-20.csv", header = T, row.names = 1)
mnorm <- read.csv("D:/GEODES_mapping_summaries/Mendota_normalized_counts_2017-06-20.csv", header = T, row.names = 1)

spvals <- read.csv("D:/GEODES_mapping_summaries/Sparkling_pvalues_2017-06-20.csv", header = T, row.names = 1)
tpvals <- read.csv("D:/GEODES_mapping_summaries/TroutBog_pvalues_2017-06-20.csv", header = T, row.names = 1)
mpvals <- read.csv("D:/GEODES_mapping_summaries/Mendota_pvalues_2017-06-20.csv", header = T, row.names = 1)

sfold <- read.csv("D:/GEODES_mapping_summaries/Sparkling_foldchange_2017-06-20.csv", header = T, row.names = 1)
tfold <- read.csv("D:/GEODES_mapping_summaries/TroutBog_foldchange_2017-06-20.csv", header = T, row.names = 1)
mfold <- read.csv("D:/GEODES_mapping_summaries/Mendota_foldchange_2017-06-20.csv", header = T, row.names = 1)

# Identify clusters by significant changes between timepoints
# Using only the top 10%

spvals_abun <- spvals[which(rowSums(snorm) > quantile(rowSums(snorm), 0.9)), ]
spvals_bin <- spvals_abun < 0.05
spvals_bin <- spvals_bin[which(rowSums(spvals_bin) > 0), ]
spvals_bin <- spvals_bin*1
sparkling_kmeans <- kmeans(spvals_bin, centers = 20)

mpvals_abun <- mpvals[which(rowSums(mnorm) > quantile(rowSums(mnorm), 0.9)), ]
mpvals_bin <- mpvals_abun < 0.05
mpvals_bin <- mpvals_bin[which(rowSums(mpvals_bin) > 0), ]
mpvals_bin <- mpvals_bin*1
mendota_kmeans <- kmeans(mpvals_bin, centers = 20)

tpvals_abun <- tpvals[which(rowSums(tnorm) > quantile(rowSums(tnorm), 0.9)), ]
tpvals_bin <- tpvals_abun < 0.05
tpvals_bin <- tpvals_bin[which(rowSums(tpvals_bin) > 0), ]
tpvals_bin <- tpvals_bin*1
trout_kmeans <- kmeans(tpvals_bin, centers = 20)


# Plot information about clusters

# What is the trend of each cluster?

for(i in 1:20){
  # get read counts
  cluster_members <- names(sparkling_kmeans$cluster[which(sparkling_kmeans$cluster == i)])
  cluster_abundances <- snorm[match(cluster_members, rownames(snorm)), ]
  # z-score normalize
  z_abun <- zscore(cluster_abundances)
  # make long dataset for ggplot
  z_long <- melt(z_abun)
  #switch sample names to timpoints
  z_long$Timepoint <- metadata$Timepoint[match(z_long$Var2, rownames(metadata))]
  #Plot as line graph
  p <- ggplot(z_long, aes(x = Timepoint, y = value, group = Var1)) + geom_point(alpha = 0.2) + stat_summary(aes(group = NULL), fun.y = mean, geom = "line", color = "blue") + labs(y = "Normalized read count", title = paste("Sparkling Cluster", i))
  assign(paste("p", i, sep = ""), p)
}

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, nrow = 5)

for(i in 1:20){
  # get read counts
  cluster_members <- names(mendota_kmeans$cluster[which(mendota_kmeans$cluster == i)])
  cluster_abundances <- mnorm[match(cluster_members, rownames(mnorm)), ]
  # z-score normalize
  z_abun <- zscore(cluster_abundances)
  # make long dataset for ggplot
  z_long <- melt(z_abun)
  #switch sample names to timpoints
  z_long$Timepoint <- metadata$Timepoint[match(z_long$Var2, rownames(metadata))]
  #Plot as line graph
  p <- ggplot(z_long, aes(x = Timepoint, y = value, group = Var1)) + geom_point(alpha = 0.2) + stat_summary(aes(group = NULL), fun.y = mean, geom = "line", color = "blue") + labs(y = "Normalized read count", title = paste("Mendota Cluster", i))
  assign(paste("p", i, sep = ""), p)
}

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, nrow = 5)

for(i in 1:20){
  # get read counts
  cluster_members <- names(trout_kmeans$cluster[which(trout_kmeans$cluster == i)])
  cluster_abundances <- tnorm[match(cluster_members, rownames(tnorm)), ]
  # z-score normalize
  z_abun <- zscore(cluster_abundances)
  # make long dataset for ggplot
  z_long <- melt(z_abun)
  #switch sample names to timpoints
  z_long$Timepoint <- metadata$Timepoint[match(z_long$Var2, rownames(metadata))]
  #Plot as line graph
  p <- ggplot(z_long, aes(x = Timepoint, y = value, group = Var1)) + geom_point(alpha = 0.2) + stat_summary(aes(group = NULL), fun.y = mean, geom = "line", color = "blue") + labs(y = "Normalized read count", title = paste("Trout Cluster", i))
  assign(paste("p", i, sep = ""), p)
}

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, nrow = 5)

# What are the phylogenies and annotations of genes in each cluster?

whos_there <- function(lake, cluster_number){
  if(lake == "Sparkling"){
    results <- sparkling_kmeans
  }else if(lake == "Mendota"){
    results <- mendota_kmeans
  }else{
    results <- trout_kmeans
  }
  cluster_members <- names(results$cluster[which(results$cluster == cluster_number)])
  output <- genekey[match(cluster_members, genekey$V1), ]
  print(output)
}

# Plot expression levels of genomes over time
snorm_genomes <- genekey$V2[match(rownames(snorm), genekey$V1)]
uniq_snorm_genomes <- unique(snorm_genomes)

snorm_genome_table <- snorm[which(snorm_genomes == "pFN18A_DNA_transcript"),]
for(i in 1:length(uniq_snorm_genomes)){
  new_row <- rowSums(snorm[which(snorm_genomes == uniq_snorm_genomes[i]), ])
  snorm_genome_table <- rbind(snorm_genome_table, new_row)
}
snorm_genome_table <- snorm_genome_table[2:260,]
rownames(snorm_genome_table) <- uniq_snorm_genomes
snorm_genome_table <- snorm_genome_table[which(rownames(snorm_genome_table) != "pFN18A_DNA_transcript"),]
s_genomes <- snorm_genome_table[which(rowSums(snorm_genome_table) > quantile(rowSums(snorm_genome_table), 0.75)),]

#Keep only the maximum for each row
max_time1 <- c()
max_value1 <- c()
max_time2 <- c()
max_value2 <- c()
for(i in 1:dim(s_genomes)[1]){
  max1 <- which(s_genomes[i, ] == max(s_genomes[i, 1:17]))
  #print(length(max1))
  if(length(max1) == 1){
    max_time1[i] <- colnames(s_genomes)[max1] 
    max_value1[i] <- s_genomes[i, max1]
  }else{
    max1 <- max1[1]
    max_time1[i] <- colnames(s_genomes)[max1] 
    max_value1[i] <- s_genomes[i, max1]
  }
  max2 <- which(s_genomes[i, ] == max(s_genomes[i, 18:34]))
  #print(length(max2))
  if(length(max2) == 1){
    max_time2[i] <- colnames(s_genomes)[max2] 
    max_value2[i] <- s_genomes[i, max2]
  }else{
    max2 <- max2[1]
    max_time2[i] <- colnames(s_genomes)[max2] 
    max_value2[i] <- s_genomes[i, max2]
  }
}

max_time <- c(max_time1, max_time2)
max_value <- c(max_value1, max_value2)
max_genomes <- rep(rownames(s_genomes), 2)

sp_max <- data.frame(max_genomes, max_time, max_value)
colnames(sp_max) <- c("Genome", "Sample", "ReadCount")
sp_max$Timepoint <- metadata$Timepoint[match(sp_max$Sample, rownames(metadata))]
sp_max$Genome <- factor(sp_max$Genome, levels = sp_max$Genome[order(sp_max$Timepoint[1:65])])

sp_max$Class <- MAG_data$Class[match(sp_max$Genome, rownames(MAG_data))]
sp_max$Order <- MAG_data$Order[match(sp_max$Genome, rownames(MAG_data))]
sp_max$Phylum <- MAG_data$Phylum[match(sp_max$Genome, rownames(MAG_data))]

ggplot(sp_max, aes(x = Timepoint, y = Genome, size = log(ReadCount), color = Phylum)) + geom_point() + scale_color_brewer(palette = "Set2")

#Repeat with Mendota
mnorm_genomes <- genekey$V2[match(rownames(mnorm), genekey$V1)]
uniq_mnorm_genomes <- unique(mnorm_genomes)

mnorm_genome_table <- mnorm[which(mnorm_genomes == "pFN18A_DNA_transcript"),]
for(i in 1:length(uniq_mnorm_genomes)){
  new_row <- rowSums(mnorm[which(mnorm_genomes == uniq_mnorm_genomes[i]), ])
  mnorm_genome_table <- rbind(mnorm_genome_table, new_row)
}
mnorm_genome_table <- mnorm_genome_table[2:dim(mnorm_genome_table)[1],]
rownames(mnorm_genome_table) <- uniq_mnorm_genomes
mnorm_genome_table <- mnorm_genome_table[which(rownames(mnorm_genome_table) != "pFN18A_DNA_transcript"),]
m_genomes <- mnorm_genome_table[which(rowSums(mnorm_genome_table) > quantile(rowSums(mnorm_genome_table), 0.75)),]

#Keep only the maximum for each row
max_time1 <- c()
max_value1 <- c()
max_time2 <- c()
max_value2 <- c()
for(i in 1:dim(m_genomes)[1]){
  max1 <- which(m_genomes[i, ] == max(m_genomes[i, 1:17]))
  #print(length(max1))
  if(length(max1) == 1){
    max_time1[i] <- colnames(m_genomes)[max1] 
    max_value1[i] <- m_genomes[i, max1]
  }else{
    max1 <- max1[1]
    max_time1[i] <- colnames(m_genomes)[max1] 
    max_value1[i] <- m_genomes[i, max1]
  }
  max2 <- which(m_genomes[i, ] == max(m_genomes[i, 18:31]))
  #print(length(max2))
  if(length(max2) == 1){
    max_time2[i] <- colnames(m_genomes)[max2] 
    max_value2[i] <- m_genomes[i, max2]
  }else{
    max2 <- max2[1]
    max_time2[i] <- colnames(m_genomes)[max2] 
    max_value2[i] <- m_genomes[i, max2]
  }
}

max_time <- c(max_time1, max_time2)
max_value <- c(max_value1, max_value2)
max_genomes <- rep(rownames(m_genomes), 2)

me_max <- data.frame(max_genomes, max_time, max_value)
colnames(me_max) <- c("Genome", "Sample", "ReadCount")
me_max$Timepoint <- metadata$Timepoint[match(me_max$Sample, rownames(metadata))]
me_max$Genome <- factor(me_max$Genome, levels = me_max$Genome[order(me_max$Timepoint[1:65])])

me_max$Class <- MAG_data$Class[match(me_max$Genome, rownames(MAG_data))]
me_max$Order <- MAG_data$Order[match(me_max$Genome, rownames(MAG_data))]
me_max$Phylum <- MAG_data$Phylum[match(me_max$Genome, rownames(MAG_data))]

ggplot(me_max, aes(x = Timepoint, y = Genome, size = log(ReadCount), color = Phylum)) + geom_point() + scale_color_brewer(palette = "Set1")

#Repeat with Trout
tnorm_genomes <- genekey$V2[match(rownames(tnorm), genekey$V1)]
uniq_tnorm_genomes <- unique(tnorm_genomes)

tnorm_genome_table <- tnorm[which(tnorm_genomes == "pFN18A_DNA_transcript"),]
for(i in 1:length(uniq_tnorm_genomes)){
  new_row <- rowSums(tnorm[which(tnorm_genomes == uniq_tnorm_genomes[i]), ])
  tnorm_genome_table <- rbind(tnorm_genome_table, new_row)
}
tnorm_genome_table <- tnorm_genome_table[2:dim(tnorm_genome_table)[1],]
rownames(tnorm_genome_table) <- uniq_tnorm_genomes
tnorm_genome_table <- tnorm_genome_table[which(rownames(tnorm_genome_table) != "pFN18A_DNA_transcript"),]
t_genomes <- tnorm_genome_table[which(rowSums(tnorm_genome_table) > quantile(rowSums(tnorm_genome_table), 0.75)),]

#Keep only the maximum for each row
max_time1 <- c()
max_value1 <- c()
max_time2 <- c()
max_value2 <- c()
for(i in 1:dim(t_genomes)[1]){
  max1 <- which(t_genomes[i, ] == max(t_genomes[i, 1:20]))
  #print(length(max1))
  if(length(max1) == 1){
    max_time1[i] <- colnames(t_genomes)[max1] 
    max_value1[i] <- t_genomes[i, max1]
  }else{
    max1 <- max1[1]
    max_time1[i] <- colnames(t_genomes)[max1] 
    max_value1[i] <- t_genomes[i, max1]
  }
  max2 <- which(t_genomes[i, ] == max(t_genomes[i, 21:22]))
  #print(length(max2))
  if(length(max2) == 1){
    max_time2[i] <- colnames(t_genomes)[max2] 
    max_value2[i] <- t_genomes[i, max2]
  }else{
    max2 <- max2[1]
    max_time2[i] <- colnames(t_genomes)[max2] 
    max_value2[i] <- t_genomes[i, max2]
  }
}

max_time <- c(max_time1, max_time2)
max_value <- c(max_value1, max_value2)
max_genomes <- rep(rownames(t_genomes), 2)

tb_max <- data.frame(max_genomes, max_time, max_value)
colnames(tb_max) <- c("Genome", "Sample", "ReadCount")
tb_max$Timepoint <- metadata$Timepoint[match(tb_max$Sample, rownames(metadata))]
tb_max$Genome <- factor(tb_max$Genome, levels = tb_max$Genome[order(tb_max$Timepoint[1:56])])

tb_max$Class <- MAG_data$Class[match(tb_max$Genome, rownames(MAG_data))]
tb_max$Order <- MAG_data$Order[match(tb_max$Genome, rownames(MAG_data))]
tb_max$Phylum <- MAG_data$Phylum[match(tb_max$Genome, rownames(MAG_data))]

ggplot(tb_max, aes(x = Timepoint, y = Genome, size = log(ReadCount), color = Phylum)) + geom_point() + scale_color_brewer(palette = "Set1")

# Try with genes grouped by genome instead of summed genome
snorm_abun <- snorm[which(rowSums(snorm) > quantile(rowSums(snorm), 0.9)), ]
snorm_genomes <- genekey$V2[match(rownames(snorm_abun), genekey$V1)]

#Keep only the maximum for each row
max_time1 <- c()
max_value1 <- c()
max_time2 <- c()
max_value2 <- c()
for(i in 1:dim(snorm_abun)[1]){
  max1 <- which(snorm_abun[i, ] == max(snorm_abun[i, 1:17]))
  #print(length(max1))
  if(length(max1) == 1){
    max_time1[i] <- colnames(snorm_abun)[max1] 
    max_value1[i] <- snorm_abun[i, max1]
  }else{
    max1 <- max1[1]
    max_time1[i] <- colnames(snorm_abun)[max1] 
    max_value1[i] <- snorm_abun[i, max1]
  }
  max2 <- which(snorm_abun[i, ] == max(snorm_abun[i, 18:34]))
  #print(length(max2))
  if(length(max2) == 1){
    max_time2[i] <- colnames(snorm_abun)[max2] 
    max_value2[i] <- snorm_abun[i, max2]
  }else{
    max2 <- max2[1]
    max_time2[i] <- colnames(snorm_abun)[max2] 
    max_value2[i] <- snorm_abun[i, max2]
  }
  if(i %% 10000 == 0){
    print(i)
  }
}

max_time <- c(max_time1, max_time2)
max_value <- c(max_value1, max_value2)
max_genes <- rep(rownames(snorm_abun), 2)

sp_max <- data.frame(max_genes, max_time, max_value)
colnames(sp_max) <- c("Gene", "Sample", "ReadCount")
sp_max$Timepoint <- metadata$Timepoint[match(sp_max$Sample, rownames(metadata))]
sp_max$Genome <- rep(snorm_genomes, 2)
sp_max$Genome <- factor(sp_max$Genome, levels = sp_max$Genome[order(sp_max$Timepoint[1:length(snorm_genomes)])])

sp_max$Class <- MAG_data$Class[match(sp_max$Genome, rownames(MAG_data))]
sp_max$Order <- MAG_data$Order[match(sp_max$Genome, rownames(MAG_data))]
sp_max$Phylum <- MAG_data$Phylum[match(sp_max$Genome, rownames(MAG_data))]

ggplot(sp_max, aes(x = Timepoint, y = Genome, color = Phylum)) + geom_point(alpha = 0.25, size = 5) + scale_color_brewer(palette = "Set1")

#Repeat with other lakes
# Mendota
mnorm_abun <- mnorm[which(rowSums(mnorm) > quantile(rowSums(mnorm), 0.9)), ]
mnorm_genomes <- genekey$V2[match(rownames(mnorm_abun), genekey$V1)]

#Keep only the maximum for each row
max_time1 <- c()
max_value1 <- c()
max_time2 <- c()
max_value2 <- c()
for(i in 1:dim(mnorm_abun)[1]){
  max1 <- which(mnorm_abun[i, ] == max(mnorm_abun[i, 1:17]))
  #print(length(max1))
  if(length(max1) == 1){
    max_time1[i] <- colnames(mnorm_abun)[max1] 
    max_value1[i] <- mnorm_abun[i, max1]
  }else{
    max1 <- max1[1]
    max_time1[i] <- colnames(mnorm_abun)[max1] 
    max_value1[i] <- mnorm_abun[i, max1]
  }
  max2 <- which(mnorm_abun[i, ] == max(mnorm_abun[i, 18:31]))
  #print(length(max2))
  if(length(max2) == 1){
    max_time2[i] <- colnames(mnorm_abun)[max2] 
    max_value2[i] <- mnorm_abun[i, max2]
  }else{
    max2 <- max2[1]
    max_time2[i] <- colnames(mnorm_abun)[max2] 
    max_value2[i] <- mnorm_abun[i, max2]
  }
  if(i %% 10000 == 0){
    print(i)
  }
}

max_time <- c(max_time1, max_time2)
max_value <- c(max_value1, max_value2)
max_genes <- rep(rownames(mnorm_abun), 2)

me_max <- data.frame(max_genes, max_time, max_value)
colnames(me_max) <- c("Gene", "Sample", "ReadCount")
me_max$Timepoint <- metadata$Timepoint[match(me_max$Sample, rownames(metadata))]
me_max$Genome <- rep(mnorm_genomes, 2)
me_max$Genome <- factor(me_max$Genome, levels = me_max$Genome[order(me_max$Timepoint[1:length(mnorm_genomes)])])

me_max$Class <- MAG_data$Class[match(me_max$Genome, rownames(MAG_data))]
me_max$Order <- MAG_data$Order[match(me_max$Genome, rownames(MAG_data))]
me_max$Phylum <- MAG_data$Phylum[match(me_max$Genome, rownames(MAG_data))]

ggplot(me_max, aes(x = Timepoint, y = Genome, color = Phylum)) + geom_point(alpha = 0.25, size = 5) + scale_color_brewer(palette = "Set1")


#Trout
tnorm_abun <- tnorm[which(rowSums(tnorm) > quantile(rowSums(tnorm), 0.9)), ]
tnorm_genomes <- genekey$V2[match(rownames(tnorm_abun), genekey$V1)]

#Keep only the maximum for each row
max_time1 <- c()
max_value1 <- c()
max_time2 <- c()
max_value2 <- c()
for(i in 1:dim(tnorm_abun)[1]){
  max1 <- which(tnorm_abun[i, ] == max(tnorm_abun[i, 1:20]))
  #print(length(max1))
  if(length(max1) == 1){
    max_time1[i] <- colnames(tnorm_abun)[max1] 
    max_value1[i] <- tnorm_abun[i, max1]
  }else{
    max1 <- max1[1]
    max_time1[i] <- colnames(tnorm_abun)[max1] 
    max_value1[i] <- tnorm_abun[i, max1]
  }
  max2 <- which(tnorm_abun[i, ] == max(tnorm_abun[i, 21:22]))
  #print(length(max2))
  if(length(max2) == 1){
    max_time2[i] <- colnames(tnorm_abun)[max2] 
    max_value2[i] <- tnorm_abun[i, max2]
  }else{
    max2 <- max2[1]
    max_time2[i] <- colnames(tnorm_abun)[max2] 
    max_value2[i] <- tnorm_abun[i, max2]
  }
  if(i %% 10000 == 0){
    print(i)
  }
}

max_time <- c(max_time1, max_time2)
max_value <- c(max_value1, max_value2)
max_genes <- rep(rownames(tnorm_abun), 2)

tb_max <- data.frame(max_genes, max_time, max_value)
colnames(tb_max) <- c("Gene", "Sample", "ReadCount")
tb_max$Timepoint <- metadata$Timepoint[match(tb_max$Sample, rownames(metadata))]
tb_max$Genome <- rep(tnorm_genomes, 2)
tb_max$Genome <- factor(tb_max$Genome, levels = tb_max$Genome[order(tb_max$Timepoint[1:length(tnorm_genomes)])])

tb_max$Class <- MAG_data$Class[match(tb_max$Genome, rownames(MAG_data))]
tb_max$Order <- MAG_data$Order[match(tb_max$Genome, rownames(MAG_data))]
tb_max$Phylum <- MAG_data$Phylum[match(tb_max$Genome, rownames(MAG_data))]

ggplot(tb_max, aes(x = Timepoint, y = Genome, color = Phylum)) + geom_point(alpha = 0.25, size = 5) + scale_color_brewer(palette = "Set1")
