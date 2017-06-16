#Kraken analysis results
library(ggplot2)
library(cowplot)

metadata_file <- "C:/Users/Alex/Desktop/geodes/analyses/05R_calculations/sample_metadata.csv"
path2results <- "D:/GEODES_kraken/"

#From the metadata file, make a list of files to upload
samplenames <- read.csv(metadata_file, header = T)
samples <- as.character(samplenames$Sample)
samples[which(samplenames$site == "ME_epi")] <- paste(samples[which(samplenames$site == "ME_epi")], ".len150", sep = "")


for (i in 1:length(samples)){
  results <- read.table(paste(path2results, samples[i], "_nonrRNA.report", sep = ""), header = F, sep = "\t")
  colnames(results) <- c("Percent_reads", "Num_reads_root", "Num_reads_direct", "Rank_code", "NCBI_taxonomy", "Scientific_name")
  assign(samples[i], results)
}
order_data <- GEODES001[which(GEODES001$Rank_code == "O" & GEODES001$Percent_reads > 0), c(1, 6)]
order_data$Sample <- rep(samples[1], dim(order_data)[1])

for(j in 2:6){
  datapart <- get(samples[j])
  sub_order_data <- datapart[which(datapart$Rank_code == "O" & datapart$Percent_reads > 0), c(1, 6)]
  sub_order_data$Sample <- rep(samples[j], dim(sub_order_data)[1])
  order_data <- rbind(order_data, sub_order_data)
}

order_data$Lake <- samplenames$Lake[match(order_data$Sample, samples)]
agg_class <- aggregate(Percent_reads ~ Lake + Scientific_name, data = order_data, FUN = sum)

#num_ME <- length(which(samplenames$site == "ME_epi"))
#num_TE <- length(which(samplenames$site == "TB_epi"))
num_SP <- length(which(samplenames$Lake == "Sparkling"))

#agg_class$Percent_reads[which(agg_class$Lake == "ME_epi")] <- agg_class$Percent_reads[which(agg_class$Lake == "ME_epi")]/num_ME
#agg_class$Percent_reads[which(agg_class$Lake == "TB_epi")] <- agg_class$Percent_reads[which(agg_class$Lake == "TB_epi")]/num_TE
agg_class$Percent_reads[which(agg_class$Lake == "Sparkling")] <- agg_class$Percent_reads[which(agg_class$Lake == "Sparkling")]/num_SP

ggplot(data = agg_class, aes(x = Lake, y = Percent_reads, fill = Scientific_name)) + geom_bar(stat = "identity")

