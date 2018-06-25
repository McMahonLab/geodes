# Test FW fave genomes

metadata <- read.csv(file = "C:/Users/Goose and Gander/Desktop/geodes/bioinformatics_workflow/R_processing/sample_metadata.csv", header = T)

read_counts <- read.csv("C:/Users/Goose and Gander/Desktop/GEODES_refMAGsSAGs_mapping/GEODES2refMAGsSAGs_normalized_readcounts.csv", header = T, row.names = 1)

genekey <- read.csv("C:/Users/Goose and Gander/Desktop/GEODES_refMAGsSAGs_mapping/GEODES2refMAGsSAGS_geneinfo.txt", header = T, row.names = 1)
colnames(read_counts) <- gsub(".nonrRNA", "", colnames(read_counts))
read_counts$GEODES014 <- NULL
read_counts$GEODES033 <- NULL
read_counts$GEODES065 <- NULL
read_counts$GEODES158 <- NULL

genome <- "2545555807"
lake <- "Sparkling"
lakekey <- metadata$Lake[match(colnames(read_counts), metadata$Sample)]

#What genes are expressed?

genomekey <- genekey[which(as.character(genekey$X1) == genome), ]
genomekey$Totals <- rowSums(read_counts[match(rownames(genomekey), rownames(read_counts)), which(lakekey == lake)])
genomekey <- genomekey[order(genomekey$Totals, decreasing = T), ]
genomekey[1:100, ]
genomekey[grep("transport", genomekey$X3), ]
barplot(zscore(colSums(read_counts[match(rownames(genomekey), rownames(read_counts)), which(lakekey == lake)])), las = 2)

#Mendota
#acI-A6, 2602042079 - monosaccharide, sugar, dipeptide, carbohydrate, amino acid, iron, heavy metal cations, potassium, ammonium, nucleoside, branched amino, daunorubicin resistance, cobalt, multidrug, maltose, lipid, sperm/put, tricarboxylate, xylose, nicotinamide mononucleotide. Also, rhodopsin
#LD28, 2524023178 - phosphate, amino acid, multidrug, magnesium/cobalt, peptide/nickel, potassium. Also, rhodopsin, PQQ-dehydrogenase, Sox
#acI-B1, 2236661002 - sugar, amino acid, phosphate, iron, branched amino, sperm/put, peptide/nickel, multidrug, daunorubicin, polysaccharide/polyol, cobalt, maltose. Also - bacteriorhodopsin
#LD12, 2236347069 - mannitol, amino acid, branched amino, potassium, phosphate, C4-dicarboxylate, phosphate, polysaccharide-polyol, antibiotic, antimicrobial peptide, ammonium. Also - rhodopsin, sulfite oxidase, malate synthase

#Trout Bog
#acI-B1, 2293339186 - xylose, monosaccharide, peptide/nickel, ammonium, potassium, iron, benzoate, antibiotic, hemin, molybdenum, heavy metal cations, branched amino, sperm/put, molybdate, multidrug. Also, CO dehydrogenase, g3p dehydrogenase. NO rhodopsin
#LD28, 2524023178 - not much expressed here actually.
#PnecC - 2519103005 - also not much expressed.

#Sparkling
#acI-B1, 2236661009 - sugar, branched amino, nicotinamide mononucleotide, phosphate, ribose/xylose/arabinose/galactoside, sperm/put, multidrug, ribose, molybdenum, magnesium/cobalt. also - rhodopsin, g3p dehydrogenase, rhamnulose-1-phosphate aldolase, xylulose-5-phosphate
# acI-A7, 22642465190 - iron, branched amino, xylose, amino acid, sugar, oligopeptide, dipeptide, ammonium, phosphate, cations, sperm/put, magnesium/cobalt, multidrug, potassium, peptide/nickel, nicotinamide mononucleotide, polysaccharide/polyol. Also, rhodopsin, IG-like domain
# 
# acI-A7, 2236661001 - sugar, branched amino, oligopeptide, ammonium, multidrug, peptide/nickel, iron, amino acid, proline/glycine betaine, sperm/put, phosphate, cations. Also, NO rhodopsin
# 
# acI-A7, 2236661005 - sugar, branched amino, xylose, oligopeptide, multidrug, amino acid, drug resistance, sperm/put, melibiose, ammonium, phosphate, peptide/nickel. Also, rhodopsin, g3p dehydrogenase, beta-galactosidase, rhodopsin
#acI-A6, 2602042079 - sugar, monosaccharide, amino acid, carbohydrate, dipeptide, sperm/put, branched amino, amino acid, nucleoside, iron, xylose, nicotinamide mononucleotide, cobalt, lipid, peptide/nickel, daunorubicin, potassium, carbohydrate, tricarboxylate, multidrug, ammonium, maltose, heavy metals. Also, fucose isomerase, rhodopsin.
#Lhab-A1, 2524023174 - branched amino, C4-dicarboxylate, amino acid, type I secretion, chromate, nitrate/sulfonate/bicarbonate, carbohydrate, cation, phosphate, peptide, biopolymer, peptide/nickel, nitrate/nitrite, iron, melibiose, sperm/put, magnesium/cobalt, ammonium. Also, cyanophycin synthase.
#acI-A5, 2606217190 - g3p, sulfur, carbohydrate, sugar, nitrate/sulfonate/bicarbonate, ammonium, branched amino, sperm/put, daunorubicin resistance, multidrug, cellobiose, phosphate, siderophore, nucleoside, cobalt, glutathione, peptide/nickel, molybdenum, molybdenate, nitrate/nitrite. Also, beta-galactosidase/glucosidase.
#LD12, 2236661000 - mannitol, amino acid, branched amino, potassium, C4-dicarboxylate, phosphate, multidrug, magnesium, cation, long-chain fatty acid. Also, rhodopsin, rhodanese
#LD12, 2236347069 - phosphate, amino acid, branched amino, antibiotic, magnesium, multidrug, antimicrobial peptide, cation, long-chain fatty acid, C4-dicarboxylate, ammonium. Also, rhodopsin
#bacI-A1, 2545555807 - ion, magnesium/cobalt, multidrug, malonate, amino acid/polyamine, iron, anitmicrobial peptide, phosphate. Also, rhodopsin.

# Are the Sparkling acI-A7s really different genomes? Yes, based on ANI.
# 2264265190 vs 2236661001 92% r2 = 0.99
# 2264265190 vs 2236661005 86% r2 = 0.99
# 2236661001 vs 2236661005 88% r2 = 0.98

# Are they correlated? report above

lake <- "Sparkling"
genome <- "2236347069"
genomekey <- genekey[which(as.character(genekey$X1) == genome), ]
group1 <- colSums(read_counts[match(rownames(genomekey), rownames(read_counts)), which(lakekey == lake)])
group2 <- colSums(read_counts[match(rownames(genomekey), rownames(read_counts)), which(lakekey == lake)])
group3 <- colSums(read_counts[match(rownames(genomekey), rownames(read_counts)), which(lakekey == lake)])

cor(group1, group2) 
cor(group1, group3)
cor(group2, group3)

#completeness?
# 2264265190 48%
# 2236661001 34%
# 2236661005 80% (Damariscotta)
