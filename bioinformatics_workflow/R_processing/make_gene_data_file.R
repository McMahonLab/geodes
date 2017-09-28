# Make a metdata file once and for all

sp_norm <- read.csv("E:/Lake_data/GEODES_datafiles/Sparkling_normalized_counts_2017-06-20.csv", header = T, row.names = 1)
find_genes <- rownames(sp_norm)
rm(sp_norm)

me_norm <- read.csv("E:/Lake_data/GEODES_datafiles/Mendota_normalized_counts_2017-06-20.csv", header = T, row.names = 1)
find_genes <- append(find_genes, rownames(me_norm), length(find_genes))
rm(me_norm)

tb_norm <- read.csv("E:/Lake_data/GEODES_datafiles/TroutBog_normalized_counts_2017-06-20.csv", header = T, row.names = 1)
find_genes <- append(find_genes, rownames(tb_norm), length(find_genes))
rm(tb_norm)

find_genes <- unique(find_genes)
#Extract contig data
contig_data <- read.table("C:/Users/amlin/Desktop/geodes/analyses/03process_mapping_results/contig_metadata.txt", colClasses = c("character"))
type1 <- sapply(strsplit(find_genes,"_"), `[`, 1)
search1 <- pmatch(type1, contig_data, duplicates.ok = T)

write.table("E:/Lake_data/GEODES_datafiles/unique_genes_2017-06-22.txt", data.frame(find_genes))

#Link to genome data where possible

#Get product names

# For MAGs, get them from the gff file

# For MGs, get them from the product names file

# Get phylogenies

# For MAGS get them from the README

# For MGS, get them from the phylodist files

#Combine into a final dataframe and save

