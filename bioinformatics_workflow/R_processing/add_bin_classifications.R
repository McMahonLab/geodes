# Update gene key files with classifications from binning

# Gene keys
mendota_key <- read.csv("D:/geodes_data_tables/Mendota_ID90_genekey.csv", header = T)
spark_key <- read.csv("D:/geodes_data_tables/Sparkling_ID90_genekey.csv", header = T)
trout_key <- read.csv("D:/geodes_data_tables/Trout_ID90_genekey.csv", header = T)
# Binning results
bins <- read.csv("D:/geodes_data_tables/GEODES_bin_data.csv", header = T)
contigs <- read.table("D:/geodes_data_tables/GEODES_binned_contigs.txt")

# Keep only good bins

bins <- bins[which(bins$completeness > 30 & bins$contamination < 10), ]
keep <- match(contigs$V2, bins$bin)
contigs <- contigs[which(is.na(keep) == F), ]

# Replace in lake keys
mendota_key$Taxonomy <- as.character(mendota_key$Taxonomy)
search <- match(mendota_key$Genome, as.character(contigs$V1))
where_in_key <- which(is.na(search) == F)
where_in_contigs <- search[where_in_key]
matching_bins <- contigs$V2[where_in_contigs]
taxonomy_add <- as.character(bins$phylodist_taxonomy[match(matching_bins, bins$bin)])
mendota_key$Taxonomy[where_in_key] <- taxonomy_add

write.csv(mendota_key, "D:/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", row.names = F)

spark_key$Taxonomy <- as.character(spark_key$Taxonomy)
search <- match(spark_key$Genome, as.character(contigs$V1))
where_in_key <- which(is.na(search) == F)
where_in_contigs <- search[where_in_key]
matching_bins <- contigs$V2[where_in_contigs]
taxonomy_add <- as.character(bins$phylodist_taxonomy[match(matching_bins, bins$bin)])
spark_key$Taxonomy[where_in_key] <- taxonomy_add

write.csv(spark_key, "D:/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", row.names = F)

trout_key$Taxonomy <- as.character(trout_key$Taxonomy)
search <- match(trout_key$Genome, as.character(contigs$V1))
where_in_key <- which(is.na(search) == F)
where_in_contigs <- search[where_in_key]
matching_bins <- contigs$V2[where_in_contigs]
taxonomy_add <- as.character(bins$phylodist_taxonomy[match(matching_bins, bins$bin)])
trout_key$Taxonomy[where_in_key] <- taxonomy_add

write.csv(trout_key, "D:/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", row.names = F)
