geodes057 <- read.table("D:/GEODES057.contig.classification.perc70.minhit3.txt", header = T, sep = "\t", colClasses = c('character'))

geodes057_results <- table(geodes057[,1])
geodes057_results <- geodes057_results[which(geodes057_results < 100000 & geodes057_results > 10)]
geodes057_results <- sort(geodes057_results)

library(ggplot2)
library(cowplot)
plot057 <- data.frame(geodes057_results)
ggplot(plot057,aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) + coord_flip() + background_grid(major = "xy", minor = "none") + labs(title = "GEODES057 Trout Bog", x = "", y = "# Contigs")

geodes058 <- read.table("D:/GEODES058.contig.classification.perc70.minhit3.txt", header = T, sep = "\t", colClasses = c('character'))

geodes058_results <- table(geodes058[,1])
geodes058_results <- geodes058_results[which(geodes058_results < 100000 & geodes058_results > 10)]
geodes058_results <- sort(geodes058_results)

library(ggplot2)
library(cowplot)
plot058 <- data.frame(geodes058_results)
ggplot(plot058,aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) + coord_flip() + background_grid(major = "xy", minor = "none") + labs(title = "GEODES058 Trout Bog", x = "", y = "# Contigs")