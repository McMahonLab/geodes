# Use the gff file to make a metadata file to interpret rownames in the mapping counts output

# Read in gff file
gff <- read.csv("C:/Users/Alex/Desktop/geodes/analyses/02mapping/mapping_database.csv", header = T, fill = NA)

# Read in data table to interpret

metaT_table <- read.table("C:/Users/Alex/Desktop/geodes/analyses/03process_mapping_results/GEODES_genes_2017-02-22.txt", row.names=1, header = T)
metaT_table <- metaT_table[which(rowSums(metaT_table) > 0), ]

# The rownames in the table are the locus_tag field. I want the contig the gene was on (so I can figure out what genome) and gene product it encodes.
# Ideally, I'd like a table where V1 = row name, V2 = genome, V3 = gene product

# Can I write a bash script to get contig names out of the gff or fasta files for each genome?
# For each gff file, 1. get column 1 2. remove duplicates 3. add second column with basename of file 4. concatenate all the outputs together


#I'm jst going to go into interactive mode for this

#cp /mnt/gluster/amlinz/ref_genomes/gff_files/*.gff .
#for file in *.gff;do awk '{print $1}' $file > temp.txt; sort -u temp.txt > $file.txt;awk '{print FILENAME $0}' $file.txt > temp2.txt;mv temp2.txt $file.txt;done
#cat *.gff.txt > contig_metadata.txt

#Downloading output file to 03processing
contig_file <- read.table("C:/Users/Alex/Desktop/geodes/analyses/03process_mapping_results/contig_metadata.txt")
#split up the columns - first part should be the same length, which makes this easier. Doesn't work for 
contigs <- data.frame(substr(contig_file$V1, start = 1, stop = 10), substr(contig_file$V1, start = 19, stop = 100), stringsAsFactors = F)
contigs[41209, ] <- c("standard", "pFN18A_DNA_transcript")
colnames(contigs) <- c("Genome", "Contig")

rowkey <- data.frame(rownames(metaT_table), stringsAsFactors = F)
colnames(rowkey) <- c("Hits")
x <- match(rowkey$Hits, gff$Locus_tag)
rowkey$Contig <- gff$Contig[x]
rowkey$Product <- gff$Product[x]
y <- match(rowkey$Contig, contigs$Contig)
rowkey$Genome <- contigs$Genome[y]

#Last thing I want is the phylogenetic assignment of each genome - get this from the ref_MAGs_SAGs readme

genome_data <- read.csv("C:/Users/Alex/Desktop/geodes/analyses/03process_mapping_results/README.csv", header = T, colClasses = c("character"))
phylogeny <- c()
for(i in 1:dim(genome_data)[1]){
  phylogeny[i]<- paste(genome_data[i, 2:6], collapse = ";")
}

z <- match(rowkey$Genome, genome_data$IMG.OID)
rowkey$Phylogeny <- phylogeny[z]
