# source("https://bioconductor.org/biocLite.R")
# biocLite("rain")
# bioLite("DESeq2")

# Load Necessary Libraries
library(rain)
library(DESeq2)

################################################################################################
###################### Sample RAIN Diel Analysis for Mock Dataset ##############################
################################################################################################
#load data
#table1 <- read.table(file='sample_input_RAIN.txt', header=TRUE, row.names=1, sep="\t", quote="")
#Mendota
mnorm <- read.csv("E:/geodes_data_tables/Mendota_ID90_normalized_readcounts.csv", header = T, row.names = 1)
mendota_key <- read.csv("E:/geodes_data_tables/Mendota_ID90_genekey_reclassified_2018-03-03.csv", header = T)
metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

colnames(mnorm) <- gsub(".nonrRNA", "", colnames(mnorm))
timepoints <- metadata$Timepoint[match(colnames(mnorm), metadata$Sample)]
#DESeq only accepts integers - round numbers in read counts
mnorm <- round(mnorm/11)
mnorm <- mnorm[which(rowSums(mnorm) > 15000),]

# use DESeq to variance stabilize data before analysis
colData <- data.frame(condition=factor(timepoints))
dds <- DESeqDataSetFromMatrix(as.matrix(mnorm), colData, formula(~ condition))
var <- varianceStabilizingTransformation(dds, blind=T)
vsd <- assay(var)

# setup time vector for RAIN
# The following code gives a vector where the times given in "t" are represented by a 1 in the vector "ft".
# The time series in question lasted 46 hours and sampling times were approximately every 4 hours, so the 
# vector "ft" gives a 1 for each hour representing a sampling time. This would be easier if sampling was
# strictly every 4 hours, because then ft could be a vector of 12 1's and deltat=4 in the RAIN call below. 
t <- c(0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44)
ft <- t
for (i in 1:length(ft)) {
	if (ft[i] %in% t) {
		ft[i] = length(which(timepoints == ft[i]))
	}
	else {
		ft[i] = 0
	}
}


rainresult <- rain(t(vsd), period=24, measure.sequence = ft, deltat = 4)
fdr <- p.adjust(as.numeric(rainresult$pVal), method="fdr") # correct p-values for multiple testing

# How many genes are significantly periodic?
sum(fdr < 0.05)

# Get estimated peak times from smoothed profiles
peak_times = list()
for( i in 1:dim(rainresult)[1]) {
	name <- row.names(rainresult)[i]
	plot.new()
	d = xspline(c(1:12), vsd[name,], shape= -0.5 , draw=F)
	maxtime <- d$x[d$y == max(d$y)]
	maxtime_c <- (maxtime * 4) - 2
		if (maxtime_c > 24) {
			maxtime_f <- maxtime_c - 24
		}

		else {
			maxtime_f <- maxtime_c
		}
	peak_times[i] <- maxtime_f
}
peak_times <- as.numeric(peak_times)
diel <- fdr <= 0.05
table3 <- cbind(rainresult[,1], fdr, rainresult[,2:4], peak_times, diel)
colnames(table3) <- c("pVal", "fdr", "phase", "peak.shape", "period", "peak_times", "diel")
full_table = table3[order(table3$pVal),]

# Make new table of compiled gene info, rain results, and counts. 
write.table(full_table, file="Diel_Analysis.results", sep="\t", quote=F)

###############################
############ End ##############
###############################
