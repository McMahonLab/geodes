# combine each featureCounts result file into a single table

# Need list of files to read

filelist <- read.table("D:/GEODES_mapping_summaries/bamfiles.txt", colClasses = c("character"))

# Read in all count files
locus_tags <- c()
for(i in 1:length(filelist$V1)){
  name <- filelist$V1[i]
  output <- read.table(paste("D:/GEODES_mapping_summaries/", name, ".CDS.mapped.txt", sep = ""), header = T, colClasses = c("character", "character", "numeric", "numeric", "character", "numeric", "numeric"))
  assign(paste(name), output)
  locus_tags <- append(locus_tags, output$Geneid, length(locus_tags))
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
