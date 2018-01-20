# Fourier transformations!

library(ggplot2)
library(cowplot)
library(GeneCycle)

#snorm <- read.csv("D:/geodes_data_tables/Sparkling_normalized.csv", header = T, row.names = 1)
#tnorm <- read.csv("D:/geodes_data_tables/Trout_normalized.csv", header = T, row.names = 1)
mnorm <- read.csv("D:/geodes_data_tables/Mendota_normalized.csv", header = T, row.names = 1)
metadata <- read.csv(file = "C:/Users/Alex/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
  plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))
  
  # TODO: why this scaling is necessary?
  plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
  
  plot(plot.data, t="h", lwd=2, main="", 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
  return(plot.data)
}

# Mendota has the best coverage, so I'll start with that one.
# For each gene:
# - average by timepoint
# - eyeball it
# - apply fast Fourier Transformation
# - plot periodogram/spectrogram

timepoints <- metadata$Timepoint[match(colnames(mnorm), metadata$Sample)]
index = 1
input_data <- as.numeric(mnorm[index,])
to.agg <- data.frame(timepoints, input_data)
ggplot(data = to.agg, aes(x = timepoints, y = input_data)) + geom_point()
averaged <- aggregate(input_data ~ timepoints, data = to.agg, mean)
ggplot(data = averaged, aes(x = timepoints, y = input_data)) + geom_point() + geom_line()
fourier <- fft(averaged$input_data)
plot.frequency.spectrum(fourier)
gstat <- max(Mod(fourier))/sum(Mod(fourier)[1:5])

fdr.out <- fdrtool(fisher.g.test(t(mnorm[1:10000,])), statistic = "pvalue")
length(which(fdr.out$qval < 0.05))
length(which(fdr.out$pval < 0.05))

sig <- mnorm[which(fdr.out$pval < 0.05),]
index = 25
input_data <- as.numeric(sig[index,])
to.agg <- data.frame(timepoints, input_data)
ggplot(data = to.agg, aes(x = timepoints, y = input_data)) + geom_point()
averaged <- aggregate(input_data ~ timepoints, data = to.agg, mean)
ggplot(data = averaged, aes(x = timepoints, y = input_data)) + geom_point() + geom_line()
fourier <- fft(averaged$input_data)
plot.frequency.spectrum(fourier)

