##### Find specific products in nutrient cycling

library(ggplot2)
library(cowplot)
library(reshape2)
library(GeneCycle)

snorm <- read.csv("E:/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
tnorm <- read.csv("E:/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mnorm <- read.csv("E:/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

# Gene keys
mendota_key <- read.csv("E:/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
spark_key <- read.csv("E:/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
trout_key <- read.csv("E:/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# RuBisCo

rubisco_genes <- mendota_key$Gene[grep("RuBisCO|ribulose 1,5-bisphosphate carboxylase|ribulose-bisphosphate carboxylase", mendota_key$Product)]
rubisco_reads <- mnorm[match(rubisco_genes, rownames(mnorm)),]
rubisco_reads$Genes <- rownames(rubisco_reads)
rubisco_reads <- melt(rubisco_reads)
rubisco_reads$variable <- gsub(".nonrRNA", "", rubisco_reads$variable)
rubisco_reads$Timepoint <- metadata$Timepoint[match(rubisco_reads$variable, metadata$Sample)]
averaged_rubisco <- aggregate(value ~ Genes + Timepoint, data = rubisco_reads, mean)
averaged_rubisco$Taxonomy <- mendota_key$Taxonomy[match(averaged_rubisco$Genes, mendota_key$Gene)]
averaged_rubisco$Taxonomy <- gsub("unclassified Actinobacteria;unclassified Actinobacteria;", "", averaged_rubisco$Taxonomy) 
averaged_rubisco$Taxonomy <- gsub("Cyanobium;", "", averaged_rubisco$Taxonomy) 
averaged_rubisco$Taxonomy <- gsub("Synechococcus;", "", averaged_rubisco$Taxonomy) 
taxon_rubisco <- aggregate(value ~ Taxonomy + Timepoint, data = averaged_rubisco, sum)

ggplot(data = taxon_rubisco, aes(x = Timepoint, fill = Taxonomy, y = value)) + geom_line()

#Keep only order level groups
less <- c("Microcystaceae", "Synechococcaceae", "Thalassiosira")
keep <- c()
taxon <- c()
for(i in 1:length(less)){
  rows <- grep(less[i], taxon_rubisco$Taxonomy)
  keep <- append(keep, rows, length(keep))
  names <- rep(less[i], length(rows))
  taxon <- append(taxon, names, length(taxon))
}
less_table <- taxon_rubisco[keep,]
less_table$Group <- taxon
ggplot(data = less_table, aes(x = Timepoint, y = value, color = Group)) + geom_line()

# multiple sugar transport
sugar_genes <- mendota_key$Gene[grep("multiple sugar transport", mendota_key$Product)]
sugar_reads <- mnorm[match(sugar_genes, rownames(mnorm)),]
sugar_reads$Genes <- rownames(sugar_reads)
sugar_reads <- melt(sugar_reads)
sugar_reads$variable <- gsub(".nonrRNA", "", sugar_reads$variable)
sugar_reads$Timepoint <- metadata$Timepoint[match(sugar_reads$variable, metadata$Sample)]
averaged_sugar <- aggregate(value ~ Genes + Timepoint, data = sugar_reads, mean)
averaged_sugar$Taxonomy <- mendota_key$Taxonomy[match(averaged_sugar$Genes, mendota_key$Gene)]
averaged_sugar$Taxonomy <- gsub("unclassified Actinobacteria;", "", averaged_sugar$Taxonomy) 
averaged_sugar$Taxonomy <- gsub("Micrococcales;Microbacteriaceae;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Micrococcales;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Actinobacteria;Actinobacteria;;", "Bacteria;Actinobacteria;Actinobacteria;", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Corynebacteriales;Mycobacteriaceae;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Corynebacteriales;Mycobacteriaceae;Mycobacterium;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Chitinophagaceae;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Crocinitomicaceae;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Roseomonas;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Burkholderiaceae;Cupriavidus;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Burkholderiaceae;Polynucleobacter;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Burkholderiaceae;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Comamonadaceae;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Comamonadaceae;Pelomonas", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Oxalobacteraceae;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("unclassified Burkholderiales;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Cyanobium;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Synechococcus;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub("Mycobacterium;", "", averaged_sugar$Taxonomy)
averaged_sugar$Taxonomy <- gsub(";;", "", averaged_sugar$Taxonomy)
taxon_sugar <- aggregate(value ~ Timepoint + Taxonomy, data = averaged_sugar, sum)

less <- c("Actinobacteria;Actinobacteria")
keep <- c()
taxon <- c()
for(i in 1:length(less)){
  rows <- grep(less[i], taxon_sugar$Taxonomy)
  keep <- append(keep, rows, length(keep))
  names <- rep(less[i], length(rows))
  taxon <- append(taxon, names, length(taxon))
}
less_table <- taxon_sugar[keep,]
less_table$Group <- taxon
ggplot(data = less_table, aes(x = Timepoint, y = value, color = Group)) + geom_line()

# Reductive TCA
rTCA_genes <- mendota_key$Gene[grep("citrate lyase", mendota_key$Product)]
rTCA_reads <- mnorm[match(rTCA_genes, rownames(mnorm)),]
rTCA_reads$Genes <- rownames(rTCA_reads)
rTCA_reads <- melt(rTCA_reads)
rTCA_reads$variable <- gsub(".nonrRNA", "", rTCA_reads$variable)
rTCA_reads$Timepoint <- metadata$Timepoint[match(rTCA_reads$variable, metadata$Sample)]
averaged_rTCA <- aggregate(value ~ Genes + Timepoint, data = rTCA_reads, mean)
averaged_rTCA$Taxonomy <- mendota_key$Taxonomy[match(averaged_rTCA$Genes, mendota_key$Gene)]
averaged_rTCA$Taxonomy <- gsub("Micrococcales;", "", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub("Porphyrobacter;", "", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub("unclassified Actinobacteria;", "", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub("\\[Blank\\];;;", "", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub(";;;;", ";", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub("Actinomycetales;;;", ";", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub("Actinomycetales;acI;acI-A;acI-A6", "", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub("Actinomycetales;acI;acI-B;acI-B1", "", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub("Corynebacteriales;Mycobacteriaceae;Mycobacterium;", "", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub("Polaromonas;", "", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub("Massilia;", "", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub(";;", ";", averaged_rTCA$Taxonomy)
averaged_rTCA$Taxonomy <- gsub("Bacteria;Actinobacteria;Actinobacteria;", "Actinobacteria;Actinobacteria;", averaged_rTCA$Taxonomy)
taxon_rTCA <- aggregate(value ~ Taxonomy, data = averaged_rTCA, sum)
head(taxon_rTCA[order(taxon_rTCA$value, decreasing = T),],15)
taxon_rTCA <- aggregate(value ~ Taxonomy + Timepoint, data = averaged_rTCA, sum)


#Keep only order level groups
less <- c("Actinobacteria;Actinobacteria", "Xanthomonadales", "Comamonadaceae", "Oxalobacteraceae", "Erythrobacteraceae")
keep <- c()
taxon <- c()
for(i in 1:length(less)){
  rows <- grep(less[i], taxon_rTCA$Taxonomy)
  keep <- append(keep, rows, length(keep))
  names <- rep(less[i], length(rows))
  taxon <- append(taxon, names, length(taxon))
}
less_table <- taxon_rTCA[keep,]
less_table$Group <- taxon
ggplot(data = less_table, aes(x = Timepoint, y = value, color = Group)) + geom_line()

# rhamnose
rhamnose_genes <- mendota_key$Gene[grep("rhamnulose-1-phosphate", mendota_key$Product)]
rhamnose_reads <- mnorm[match(rhamnose_genes, rownames(mnorm)),]
rhamnose_reads$Genes <- rownames(rhamnose_reads)
rhamnose_reads <- melt(rhamnose_reads)
rhamnose_reads$variable <- gsub(".nonrRNA", "", rhamnose_reads$variable)
rhamnose_reads$Timepoint <- metadata$Timepoint[match(rhamnose_reads$variable, metadata$Sample)]
averaged_rhamnose <- aggregate(value ~ Genes + Timepoint, data = rhamnose_reads, mean)
averaged_rhamnose$Taxonomy <- mendota_key$Taxonomy[match(averaged_rhamnose$Genes, mendota_key$Gene)]
averaged_rhamnose$Taxonomy <- gsub("\\[Blank\\];;;", "", averaged_rhamnose$Taxonomy)
averaged_rhamnose$Taxonomy <- gsub("Actinomycetales;acI;acI-B;acI-B1", "", averaged_rhamnose$Taxonomy)
averaged_rhamnose$Taxonomy <- gsub("Bacteria;Actinobacteria;Actinobacteria", "Actinobacteria;Actinobacteria", averaged_rhamnose$Taxonomy)
averaged_rhamnose$Taxonomy <- gsub("Puniceicoccales;;;", "", averaged_rhamnose$Taxonomy)
averaged_rhamnose$Taxonomy <- gsub("Microbacteriaceae;", "", averaged_rhamnose$Taxonomy)
taxon_rhamnose <- aggregate(value ~ Taxonomy, data = averaged_rhamnose, sum)
head(taxon_rhamnose[order(taxon_rhamnose$value, decreasing = T),],15)
taxon_rhamnose <- aggregate(value ~ Taxonomy + Timepoint, data = averaged_rhamnose, sum)


#Keep only order level groups
less <- c("Spartobacteriales", "Opitutae", "Sphingobacteriales", "Actinobacteria")
keep <- c()
taxon <- c()
for(i in 1:length(less)){
  rows <- grep(less[i], taxon_rhamnose$Taxonomy)
  keep <- append(keep, rows, length(keep))
  names <- rep(less[i], length(rows))
  taxon <- append(taxon, names, length(taxon))
}
less_table <- taxon_rhamnose[keep,]
less_table$Group <- taxon
ggplot(data = less_table, aes(x = Timepoint, y = value, color = Group)) + geom_line()

#fucose
fucose_genes <- mendota_key$Gene[grep("fuculose-phosphate", mendota_key$Product)]
fucose_reads <- mnorm[match(fucose_genes, rownames(mnorm)),]
fucose_reads$Genes <- rownames(fucose_reads)
fucose_reads <- melt(fucose_reads)
fucose_reads$variable <- gsub(".nonrRNA", "", fucose_reads$variable)
fucose_reads$Timepoint <- metadata$Timepoint[match(fucose_reads$variable, metadata$Sample)]
averaged_fucose <- aggregate(value ~ Genes + Timepoint, data = fucose_reads, mean)
averaged_fucose$Taxonomy <- mendota_key$Taxonomy[match(averaged_fucose$Genes, mendota_key$Gene)]
averaged_fucose$Taxonomy <- gsub("\\[Blank\\];;;", "", averaged_fucose$Taxonomy)
averaged_fucose$Taxonomy <- gsub("Planctomycetaceae;Pirellula;", "", averaged_fucose$Taxonomy)
averaged_fucose$Taxonomy <- gsub("Planctomycetaceae;", "", averaged_fucose$Taxonomy)
averaged_fucose$Taxonomy <- gsub("unclassified Actinobacteria;", "", averaged_fucose$Taxonomy)
averaged_fucose$Taxonomy <- gsub("Bacteria;Actinobacteria;", "Actinobacteria;", averaged_fucose$Taxonomy)
averaged_fucose$Taxonomy <- gsub("Synechococcus;", "", averaged_fucose$Taxonomy)
averaged_fucose$Taxonomy <- gsub("Opitutales;", "", averaged_fucose$Taxonomy)
averaged_fucose$Taxonomy <- gsub("Rhodopirellula;", "", averaged_fucose$Taxonomy)
averaged_fucose$Taxonomy <- gsub("Micrococcales;Microbacteriaceae;", "", averaged_fucose$Taxonomy)
taxon_fucose <- aggregate(value ~ Taxonomy, data = averaged_fucose, sum)
head(taxon_fucose[order(taxon_fucose$value, decreasing = T),],15)
taxon_fucose <- aggregate(value ~ Taxonomy + Timepoint, data = averaged_fucose, sum)


#Keep only order level groups
less <- c("Planctomycetales", "Actinobacteria;Actinobacteria", "Synechococcaceae", "Opitutae", "Acetobacteraceae")
keep <- c()
taxon <- c()
for(i in 1:length(less)){
  rows <- grep(less[i], taxon_fucose$Taxonomy)
  keep <- append(keep, rows, length(keep))
  names <- rep(less[i], length(rows))
  taxon <- append(taxon, names, length(taxon))
}
less_table <- taxon_fucose[keep,]
less_table$Group <- taxon
ggplot(data = less_table, aes(x = Timepoint, y = value, color = Group)) + geom_line()

#arabinose
arabinose_genes <- mendota_key$Gene[grep("ribulokinase", mendota_key$Product)]
arabinose_reads <- mnorm[match(arabinose_genes, rownames(mnorm)),]
arabinose_reads$Genes <- rownames(arabinose_reads)
arabinose_reads <- melt(arabinose_reads)
arabinose_reads$variable <- gsub(".nonrRNA", "", arabinose_reads$variable)
arabinose_reads$Timepoint <- metadata$Timepoint[match(arabinose_reads$variable, metadata$Sample)]
averaged_arabinose <- aggregate(value ~ Genes + Timepoint, data = arabinose_reads, mean)
averaged_arabinose$Taxonomy <- mendota_key$Taxonomy[match(averaged_arabinose$Genes, mendota_key$Gene)]
averaged_arabinose$Taxonomy <- gsub("\\[Blank\\];;;", "", averaged_arabinose$Taxonomy)
averaged_arabinose$Taxonomy <- gsub("\\[Blank\\];", "", averaged_arabinose$Taxonomy)
averaged_arabinose$Taxonomy <- gsub("Microcystis;", "", averaged_arabinose$Taxonomy)
averaged_arabinose$Taxonomy <- gsub("Synechococcus;", "", averaged_arabinose$Taxonomy)
averaged_arabinose$Taxonomy <- gsub("Cyanobium;", "", averaged_arabinose$Taxonomy)
averaged_arabinose$Taxonomy <- gsub("Planctomycetaceae;Pirellula;Pirellula sp. SH-Sr6A;", "", averaged_arabinose$Taxonomy)
taxon_arabinose <- aggregate(value ~ Taxonomy, data = averaged_arabinose, sum)
head(taxon_arabinose[order(taxon_arabinose$value, decreasing = T),],15)
taxon_arabinose <- aggregate(value ~ Taxonomy + Timepoint, data = averaged_arabinose, sum)


#Keep only order level groups
less <- c("Synechococcaceae", "Microcystaceae", "Comamonadaceae", "Planctomycetales", "Chitinophagales")
keep <- c()
taxon <- c()
for(i in 1:length(less)){
  rows <- grep(less[i], taxon_arabinose$Taxonomy)
  keep <- append(keep, rows, length(keep))
  names <- rep(less[i], length(rows))
  taxon <- append(taxon, names, length(taxon))
}
less_table <- taxon_arabinose[keep,]
less_table$Group <- taxon
ggplot(data = less_table, aes(x = Timepoint, y = value, color = Group)) + geom_line()

#galactose
galactose_genes <- mendota_key$Gene[grep("aldose 1-epimerase", mendota_key$Product)]
galactose_reads <- mnorm[match(galactose_genes, rownames(mnorm)),]
galactose_reads$Genes <- rownames(galactose_reads)
galactose_reads <- melt(galactose_reads)
galactose_reads$variable <- gsub(".nonrRNA", "", galactose_reads$variable)
galactose_reads$Timepoint <- metadata$Timepoint[match(galactose_reads$variable, metadata$Sample)]
averaged_galactose <- aggregate(value ~ Genes + Timepoint, data = galactose_reads, mean)
averaged_galactose$Taxonomy <- mendota_key$Taxonomy[match(averaged_galactose$Genes, mendota_key$Gene)]
averaged_galactose$Taxonomy <- gsub("\\[Blank\\];;;", "", averaged_galactose$Taxonomy)
averaged_galactose$Taxonomy <- gsub("\\[Blank\\];", "", averaged_galactose$Taxonomy)
averaged_galactose$Taxonomy <- gsub("Microcystis;", "", averaged_galactose$Taxonomy)
averaged_galactose$Taxonomy <- gsub("Synechococcus;", "", averaged_galactose$Taxonomy)
averaged_galactose$Taxonomy <- gsub("Cyanobium;", "", averaged_galactose$Taxonomy)

taxon_galactose <- aggregate(value ~ Taxonomy, data = averaged_galactose, sum)
head(taxon_galactose[order(taxon_galactose$value, decreasing = T),],15)
taxon_galactose <- aggregate(value ~ Taxonomy + Timepoint, data = averaged_galactose, sum)


# #Keep only order level groups
# less <- c("Synechococcaceae", "Microcystaceae", "Comamonadaceae", "Planctomycetales", "Chitinophagales")
# keep <- c()
# taxon <- c()
# for(i in 1:length(less)){
#   rows <- grep(less[i], taxon_galactose$Taxonomy)
#   keep <- append(keep, rows, length(keep))
#   names <- rep(less[i], length(rows))
#   taxon <- append(taxon, names, length(taxon))
# }
# less_table <- taxon_galactose[keep,]
# less_table$Group <- taxon
# ggplot(data = less_table, aes(x = Timepoint, y = value, color = Group)) + geom_line()

#xylose
xylose_genes <- mendota_key$Gene[grep("xylose isomerase", mendota_key$Product)]
xylose_reads <- mnorm[match(xylose_genes, rownames(mnorm)),]
xylose_reads$Genes <- rownames(xylose_reads)
xylose_reads <- melt(xylose_reads)
xylose_reads$variable <- gsub(".nonrRNA", "", xylose_reads$variable)
xylose_reads$Timepoint <- metadata$Timepoint[match(xylose_reads$variable, metadata$Sample)]
averaged_xylose <- aggregate(value ~ Genes + Timepoint, data = xylose_reads, mean)
averaged_xylose$Taxonomy <- mendota_key$Taxonomy[match(averaged_xylose$Genes, mendota_key$Gene)]
averaged_xylose$Taxonomy <- gsub("\\[Blank\\];;;", "", averaged_xylose$Taxonomy)
averaged_xylose$Taxonomy <- gsub("\\[Blank\\];", "", averaged_xylose$Taxonomy)
averaged_xylose$Taxonomy <- gsub("Microcystis;", "", averaged_xylose$Taxonomy)
averaged_xylose$Taxonomy <- gsub("Synechococcus;", "", averaged_xylose$Taxonomy)
averaged_xylose$Taxonomy <- gsub("Cyanobium;", "", averaged_xylose$Taxonomy)

taxon_xylose <- aggregate(value ~ Taxonomy, data = averaged_xylose, sum)
head(taxon_xylose[order(taxon_xylose$value, decreasing = T),],15)
taxon_xylose <- aggregate(value ~ Taxonomy + Timepoint, data = averaged_xylose, sum)


# #Keep only order level groups
# less <- c("Synechococcaceae", "Microcystaceae", "Comamonadaceae", "Planctomycetales", "Chitinophagales")
# keep <- c()
# taxon <- c()
# for(i in 1:length(less)){
#   rows <- grep(less[i], taxon_xylose$Taxonomy)
#   keep <- append(keep, rows, length(keep))
#   names <- rep(less[i], length(rows))
#   taxon <- append(taxon, names, length(taxon))
# }
# less_table <- taxon_xylose[keep,]
# less_table$Group <- taxon
# ggplot(data = less_table, aes(x = Timepoint, y = value, color = Group)) + geom_line()

# Now do nitrogen metabolism
nif_genes <- mendota_key$Gene[grep("sulfate reductase", mendota_key$Product)]
nif_reads <- mnorm[match(nif_genes, rownames(mnorm)),]
nif_reads$Genes <- rownames(nif_reads)
nif_reads <- melt(nif_reads)
nif_reads$variable <- gsub(".nonrRNA", "", nif_reads$variable)
nif_reads$Timepoint <- metadata$Timepoint[match(nif_reads$variable, metadata$Sample)]
averaged_nif <- aggregate(value ~ Genes + Timepoint, data = nif_reads, mean)
averaged_nif$Taxonomy <- mendota_key$Taxonomy[match(averaged_nif$Genes, mendota_key$Gene)]


taxon_nif <- aggregate(value ~ Taxonomy, data = averaged_nif, sum)
head(taxon_nif[order(taxon_nif$value, decreasing = T),],15)
taxon_nif <- aggregate(value ~ Timepoint, data = averaged_nif, sum)

ggplot(data = taxon_nif, aes(x = Timepoint, y = value)) + geom_line()


###Make a more streamlined function for this, including other lake datasets

genesearch <- function(lake, gene){
  if(lake == "Mendota"){
    reads <- mnorm
    key <- mendota_key
  }else if(lake == "Sparkling"){
    reads <- snorm
    key <- spark_key
  }else if(lake == "Trout"){
    reads <- tnorm
    key <- trout_key
  }
  search_genes <- key$Gene[grep(gene, key$Product)]
  search_reads <- reads[match(search_genes, rownames(reads)),]
  search_reads$Genes <- rownames(search_reads)
  search_reads <- melt(search_reads)
  search_reads$variable <- gsub(".nonrRNA", "", search_reads$variable)
  search_reads$Timepoint <- metadata$Timepoint[match(search_reads$variable, metadata$Sample)]
  averaged_search <- aggregate(value ~ Genes + Timepoint, data = search_reads, mean)
  averaged_search$Taxonomy <- key$Taxonomy[match(averaged_search$Genes, key$Gene)]
  
  
  taxon_search <- aggregate(value ~ Taxonomy, data = averaged_search, sum)
  head(taxon_search[order(taxon_search$value, decreasing = T),],5)
  
}