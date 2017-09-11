geodes057 <- read.table("E:/geodes_MG/phylodist_classifications/GEODES057.contig.classification.perc70.minhit3.txt", header = T, sep = "\t", colClasses = c('character'))

geodes057_results <- table(geodes057[,2])
unclassified <- geodes057_results[which(names(geodes057_results) == "NO CLASSIFICATION MH")]/sum(geodes057_results)
print(unclassified)
geodes057_results <- geodes057_results[which(geodes057_results < 100000 & geodes057_results > 100)]
geodes057_results <- sort(geodes057_results)

library(ggplot2)
library(cowplot)
plot057 <- data.frame(geodes057_results)
p1 <- ggplot(plot057,aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) + coord_flip() + background_grid(major = "xy", minor = "none") + labs(title = "GEODES057 Trout Bog", x = "", y = "# Contigs")
save_plot("C:/Users/amlin/Desktop/geodes/GEODES057_contig_classifications.pdf", plot = p1, base_height = 10, base_aspect_ratio = 2)

######
geodes058 <- read.table("E:/geodes_MG/phylodist_classifications/GEODES058.contig.classification.perc70.minhit3.txt", header = T, sep = "\t", colClasses = c('character'))

geodes058_results <- table(geodes058[,2])
unclassified <- geodes058_results[which(names(geodes058_results) == "NO CLASSIFICATION MH")]/sum(geodes058_results)
print(unclassified)
geodes058_results <- geodes058_results[which(geodes058_results < 100000 & geodes058_results > 100)]
geodes058_results <- sort(geodes058_results)

plot058 <- data.frame(geodes058_results)
p2 <- ggplot(plot058,aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) + coord_flip() + background_grid(major = "xy", minor = "none") + labs(title = "GEODES058 Trout Bog", x = "", y = "# Contigs")
save_plot("C:/Users/amlin/Desktop/geodes/GEODES058_contig_classifications.pdf", plot = p2, base_height = 10, base_aspect_ratio = 2)

#####
geodes005 <- read.table("E:/geodes_MG/phylodist_classifications/GEODES005.contig.classification.perc70.minhit3.txt", header = T, sep = "\t", colClasses = c('character'))

geodes005_results <- table(geodes005[,2])
unclassified <- geodes005_results[which(names(geodes005_results) == "NO CLASSIFICATION MH")]/sum(geodes005_results)
print(unclassified)
geodes005_results <- geodes005_results[which(geodes005_results < 100000 & geodes005_results > 100)]
geodes005_results <- sort(geodes005_results)

plot005 <- data.frame(geodes005_results)
p3 <- ggplot(plot005,aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) + coord_flip() + background_grid(major = "xy", minor = "none") + labs(title = "GEODES005 Sparkling Lake", x = "", y = "# Contigs")
save_plot("C:/Users/amlin/Desktop/geodes/GEODES005_contig_classifications.pdf", plot = p3, base_height = 10, base_aspect_ratio = 2)

#####
geodes006 <- read.table("E:/geodes_MG/phylodist_classifications/GEODES006.contig.classification.perc70.minhit3.txt", header = T, sep = "\t", colClasses = c('character'))

geodes006_results <- table(geodes006[,2])
unclassified <- geodes006_results[which(names(geodes006_results) == "NO CLASSIFICATION MH")]/sum(geodes006_results)
print(unclassified)
geodes006_results <- geodes006_results[which(geodes006_results < 100000 & geodes006_results > 100)]
geodes006_results <- sort(geodes006_results)

plot006 <- data.frame(geodes006_results)
p4 <- ggplot(plot006,aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) + coord_flip() + background_grid(major = "xy", minor = "none") + labs(title = "GEODES006 Sparkling Lake", x = "", y = "# Contigs")
save_plot("C:/Users/amlin/Desktop/geodes/GEODES006_contig_classifications.pdf", plot = p4, base_height = 10, base_aspect_ratio = 2)

#####
geodes117 <- read.table("E:/geodes_MG/phylodist_classifications/GEODES117.contig.classification.perc70.minhit3.txt", header = T, sep = "\t", colClasses = c('character'))

geodes117_results <- table(geodes117[,2])
unclassified <- geodes117_results[which(names(geodes117_results) == "NO CLASSIFICATION MH")]/sum(geodes117_results)
print(unclassified)
geodes117_results <- geodes117_results[which(geodes117_results < 100000 & geodes117_results > 100)]
geodes117_results <- sort(geodes117_results)

plot117 <- data.frame(geodes117_results)
p5 <- ggplot(plot117,aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) + coord_flip() + background_grid(major = "xy", minor = "none") + labs(title = "GEODES117 Lake Mendota", x = "", y = "# Contigs")
save_plot("C:/Users/amlin/Desktop/geodes/GEODES117_contig_classifications.pdf", plot = p5, base_height = 10, base_aspect_ratio = 2)

#####
geodes118 <- read.table("E:/geodes_MG/phylodist_classifications/GEODES118.contig.classification.perc70.minhit3.txt", header = T, sep = "\t", colClasses = c('character'))

geodes118_results <- table(geodes118[,2])
unclassified <- geodes118_results[which(names(geodes118_results) == "NO CLASSIFICATION MH")]/sum(geodes118_results)
print(unclassified)
geodes118_results <- geodes118_results[which(geodes118_results < 100000 & geodes118_results > 100)]
geodes118_results <- sort(geodes118_results)

plot118 <- data.frame(geodes118_results)
p6 <- ggplot(plot118,aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) + coord_flip() + background_grid(major = "xy", minor = "none") + labs(title = "GEODES118 Lake Mendota", x = "", y = "# Contigs")
save_plot("C:/Users/amlin/Desktop/geodes/GEODES118_contig_classifications.pdf", plot = p6, base_height = 10, base_aspect_ratio = 2)
