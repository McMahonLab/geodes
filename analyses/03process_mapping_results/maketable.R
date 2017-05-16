# combine each featureCounts result file into a single table

# Need list of files to read

filelist <- read.table("D:/GEODES_mapping_summaries/bamfiles.txt", colClasses = c("character"))

# Read in all count files
output <- read.table(paste("D:/GEODES_mapping_summaries/", filelist$V1[1], ".CDS.txt", sep = ""), header = T, colClasses = c("character", "character", "numeric", "numeric", "character", "numeric", "numeric"))
locus_tags <- c(output$Geneid)
locus_table <- data.frame(output[, 7])
rownames(locus_table) <- locus_tags
for(i in 63:length(filelist$V1)){
  name <- filelist$V1[i]
  output <- read.table(paste("D:/GEODES_mapping_summaries/", name, ".CDS.txt", sep = ""), header = T, colClasses = c("character", "character", "numeric", "numeric", "character", "numeric", "numeric"))
  locus_table <- cbind(locus_table, output[, 7])
}

locus_tags <- unique(locus_tags)
index_col1 <- match(locus_tags, GEODES001_nonrRNA$Geneid)
col1 <- GEODES001_nonrRNA[index_col1, 7]
col1[which(is.na(col1) == T)] <- 0

locus_table <- data.frame(col1)
for(i in 2:length(filelist$V1)){
  dataset <- get(filelist$V1[i])
  column_index <- match(locus_tags, dataset$Geneid)
  column <- dataset[column_index, 7]
  column[which(is.na(column) == T)] <- 0
  locus_table <- cbind(locus_table, column)
}
colnames(locus_table) <- filelist$V1
rownames(locus_table) <- locus_tags
write.csv(locus_table, file = "C:/Users/Alex/Desktop/geodes/metaT_data/locus_tag_table_2017-04-25.csv")
