# Marker gene analysis - look mom, no clustering!

library(ggplot2)
library(cowplot)
library(reshape2)

### Load data (start with only one to save RAM and comment the rest out)
# Normalized read tables
#snorm <- read.csv("E:/geodes_data_tables/Sparkling_ID90_normalized_readcounts.csv", header = T, row.names = 1)
#tnorm <- read.csv("E:/geodes_data_tables/Trout_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mnorm <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)

# Gene keys
mendota_key <- read.csv("C:/Users/Goose and Gander/Documents/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#spark_key <- read.csv("E:/geodes_data_tables/Sparkling_ID90_genekey_reclassified_2018-03-03.csv", header = T)
#trout_key <- read.csv("E:/geodes_data_tables/Trout_ID90_genekey_reclassified_2018-03-03.csv", header = T)

# Sample data
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

# Choose your search term:
#searchterm <- c("Rubisco|RuBisCO|rubisco|bisphosphate carboxylase")
#searchterm <- c("nitrogenase|Nitrogenase|Nif")
#searchterm <- c("Chitobiase|chitobiase|chitinase|Chitinase")
#searchterm <- c("glycoside hydrolase|Glycoside hydrolase|glycosyl hydrolase")
#searchterm <- c("sulfite reductase|Sulfite reductase")
#searchterm <- c("glucosaminidase|hexosaminidase")
#searchterm <- c("putrescine|spermidine|polymamine|Putrescine|Spermidine|Polyamine")
#searchterm <- c("cobalamin|thiamin|pyridoxal")
#searchterm <- c("rhamnulose-1-phosphate|Rhamnulose-1-phosphate|fuculose-phosphate|Fuculose phosphate|ribulokinase|mannose-6-phosphate")
#searchterm <- c("aminopeptidase|Aminopeptidase")
#searchterm <- c("alkaline phosphatase|Alkaline phosphatase")
#searchterm <- c("sugar transport|ose transport")
#searchterm <- c("rhodopsin|Rhodopsin")
#searchterm <- c("phytoene|lycopene|carotene|Phytoene|Lycopene|Carotene")
#searchterm <- c("photosynth|Photosynth")
#searchterm <- c("cellulase|cellulose")
#searchterm <- c("sulfatase")
#searchterm <- c("sulfite reductase|Sulfite reductase")
#searchterm <- c("peroxidase|peroxide|catalase")
#searchterm <- c("ammonia monooxygenase|hydroxylamine oxidoreductase")
searchterm <- c("nitrite reductase|Nitrite reductase|nitrite oxidoreductase|Nitrite oxidoreductase")
#searchterm <- c("nitrate reductase|Nitrate reductase|nitrate oxidoreductase|Nitrate oxidoreductase")
#searchterm <- c("urease")

marker_genes <- mendota_key[grep(searchterm, mendota_key$Product), ]
marker_genes_mnorm <- mnorm[match(marker_genes$Gene, rownames(mnorm)), ]

#Aggregate by timepoint
marker_genes_mnorm$Genes <- rownames(marker_genes_mnorm)
marker_genes_mnorm <- melt(marker_genes_mnorm)
marker_genes_mnorm$variable <- gsub(".nonrRNA", "", marker_genes_mnorm$variable)
marker_genes_mnorm$Timepoint <- metadata$Timepoint[match(marker_genes_mnorm$variable, metadata$Sample)]
agg_marker_genes_mnorm <- aggregate(value ~ Genes + Timepoint, data = marker_genes_mnorm, mean)
new_marker_genes_mnorm <- reshape(agg_marker_genes_mnorm, idvar = "Genes", timevar = "Timepoint", direction = "wide")
rownames(new_marker_genes_mnorm) <- new_marker_genes_mnorm[, 1]
new_marker_genes_mnorm <- new_marker_genes_mnorm[, 2:dim(new_marker_genes_mnorm)[2]]

# Add taxonomy information and sum by unique taxonomy
new_marker_genes_mnorm$Genes <- rownames(new_marker_genes_mnorm)
new_marker_genes_mnorm$Taxonomy <- marker_genes$Taxonomy[match(new_marker_genes_mnorm$Genes, marker_genes$Gene)]
new_marker_genes_mnorm <- melt(new_marker_genes_mnorm)
new_marker_genes_mnorm$Taxonomy <- gsub("Bacteria;", "", new_marker_genes_mnorm$Taxonomy)
new_marker_genes_mnorm$Taxonomy <- gsub("Eukaryota;", "", new_marker_genes_mnorm$Taxonomy)
new_marker_genes_mnorm$Taxonomy <- gsub("Proteobacteria;", "", new_marker_genes_mnorm$Taxonomy)
new_marker_genes_mnorm$Phylum <- sapply(strsplit(as.character(new_marker_genes_mnorm$Taxonomy),";"), `[`, 1)
new_marker_genes_mnorm$Phylum <- gsub("Heterokonta,Coscinodiscophyceae,Thalassiosirales,Thalassiosiraceae,Thalassiosira,pseudonana", "Heterokonta", new_marker_genes_mnorm$Phylum)
#summed_markers <- aggregate(value ~ variable + Taxonomy, new_marker_genes_mnorm, sum)
summed_markers <- aggregate(value ~ variable + Phylum, new_marker_genes_mnorm, sum)
summed_markers <- summed_markers[grep("unclassified|NO CLASS|None|Blank", summed_markers$Phylum, invert = T), ]

# Plot total sums by phyla
summed_phyla <- aggregate(value ~ Phylum, new_marker_genes_mnorm, sum)
summed_phyla <- summed_phyla[grep("unclassified|NO CLASS|None|Blank", summed_phyla$Phylum, invert = T), ]
summed_phyla <- summed_phyla[which(summed_phyla$value > 1000000),]
p2 <- ggplot(data = summed_phyla, aes(x = Phylum, y = value)) + geom_bar(stat = "identity") + coord_flip() + labs(x = NULL, y = "Summed Read Counts", title = "Marker Gene Reads Assigned")

summed_markers <- summed_markers[which(summed_markers$Phylum %in% summed_phyla$Phylum),]

#Change to fold change
taxa <- unique(summed_markers$Phylum)
for(i in 1:length(taxa)){
  subset <- summed_markers[which(summed_markers$Phylum == taxa[i]), ]
  ave <- mean(subset$value)
  subset$value <- subset$value - ave
  maxrange <- max(abs(subset$value))
  subset$value <- subset$value * (10/maxrange)
  summed_markers[which(summed_markers$Phylum == taxa[i]), ] <- subset
}

p1 <- ggplot(data = summed_markers, aes(x = variable, y = Phylum, fill = value)) + geom_tile()+ scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "yellow") + scale_x_discrete(labels=c("value.0" = "5AM", "value.4" = "9AM", "value.8" = "1PM", "value.12" = "5PM", "value.16" = "9PM", "value.20" = "1AM", "value.24" = "5AM", "value.28" = "9AM", "value.32" = "1PM", "value.36" = "5PM", "value.40" = "9PM", "value.44" = "1AM")) + theme(legend.position = "none") + labs(x = NULL, y = NULL, title = paste(searchterm))

full_plot <- plot_grid(p1, p2)
full_plot

# correlation matrix
phyla <- c("Cyanobacteria", "Actinobacteria", "Bacteroidetes")
cor.matrix <- matrix(nrow = length(phyla), ncol = length(phyla), 0)
rownames(cor.matrix) <- colnames(cor.matrix) <- phyla
for(i in 1:length(phyla)){
  thing1 <- summed_markers[which(summed_markers$Phylum == phyla[i]), ]
    for(j in 1:length(phyla)){
      thing2 <- summed_markers[which(summed_markers$Phylum == phyla[j]), ]
      cor.matrix[i, j] <- cor(thing1$value, thing2$value)
    }
}
cor.matrix
