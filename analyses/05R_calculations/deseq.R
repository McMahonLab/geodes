library(ff)
library(DESeq)

# gene_table <- read.table.ffdf(file = "D:/GEODES_mapping_summaries/GEODES_genes_2017-05-23.txt.gz", header = T)
# rownames(gene_table) <- as.character(gene_table[, 1])
# gene_table <- gene_table[, 2:dim(gene_table)[2]]
# colnames(gene_table) <- gsub("_nonrRNA", "", colnames(gene_table))

metadata <- read.csv(file = "C:/Users/Alex/Desktop/geodes/analyses/05R_calculations/sample_metadata.csv", header = T, row.names = 1)
metadata$condition <- paste(metadata$Lake, metadata$Timepoint, sep = ";")
coldata <- data.frame(metadata[,4])

gene_table1 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table1.csv", header = T, row.names = 1)
gene_table2 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table2.csv", header = T, row.names = 1)
gene_table3 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table3.csv", header = T, row.names = 1)
gene_table4 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table4.csv", header = T, row.names = 1)
gene_table5 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table5.csv", header = T, row.names = 1)
gene_table6 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table6.csv", header = T, row.names = 1)
gene_table7 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table7.csv", header = T, row.names = 1)
gene_table8 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table8.csv", header = T, row.names = 1)
gene_table9 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table9.csv", header = T, row.names = 1)
gene_table10 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table10.csv", header = T, row.names = 1)
gene_table11 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table11.csv", header = T, row.names = 1)
gene_table12 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table12.csv", header = T, row.names = 1)
gene_table13 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table13.csv", header = T, row.names = 1)
gene_table17 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table17.csv", header = T, row.names = 1)
gene_table18 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table18.csv", header = T, row.names = 1)
gene_table19 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table19.csv", header = T, row.names = 1)
gene_table20 <- read.csv(file = "D:/GEODES_mapping_summaries/gene_table20.csv", header = T, row.names = 1)

sums1 <- colSums(gene_table1)
sums2 <- colSums(gene_table2)
sums3 <- colSums(gene_table3)
sums4 <- colSums(gene_table4)
sums5 <- colSums(gene_table5)
sums6 <- colSums(gene_table6)
sums7 <- colSums(gene_table7)
sums8 <- colSums(gene_table8)
sums9 <- colSums(gene_table9)
sums10 <- colSums(gene_table10)
sums11 <- colSums(gene_table11)
sums12 <- colSums(gene_table12)
sums13 <- colSums(gene_table13)
sums17 <- colSums(gene_table17)
sums18 <- colSums(gene_table18)
sums19 <- colSums(gene_table19)
sums20 <- colSums(gene_table20)

library_size_df <- rbind(sums1, sums2, sums3, sums4, sums5, sums6, sums7, sums8, sums9, sums10, sums11, sums12, sums13, sums17, sums18, sums19, sums20)
library_size <- colSums(library_size_df)
effect_size <- library_size/min(library_size)


datasets <- c("gene_table1", "gene_table2", "gene_table3", "gene_table4", "gene_table5", "gene_table6", "gene_table7", "gene_table8", "gene_table8", "gene_table9", "gene_table10", "gene_table11", "gene_table12", "gene_table13", "gene_table17", "gene_table18", "gene_table19", "gene_table20")

for(j in 17:length(datasets)){
  table <- get(datasets[j])
  colnames(table) <- gsub("_nonrRNA", "", colnames(table))
  cds <- newCountDataSet(countData = table, conditions = metadata$condition, sizeFactors = effect_size)
  cds <- estimateDispersions(cds)
  #Sparkling
  res1 <- nbinomTest(cds, "Sparkling;0", "Sparkling;4")
  res2 <- nbinomTest(cds, "Sparkling;4", "Sparkling;8")
  res3 <- nbinomTest(cds, "Sparkling;8", "Sparkling;12")
  res4 <- nbinomTest(cds, "Sparkling;12", "Sparkling;16")
  res5 <- nbinomTest(cds, "Sparkling;16", "Sparkling;20")
  res6 <- nbinomTest(cds, "Sparkling;20", "Sparkling;24")
  res7 <- nbinomTest(cds, "Sparkling;24", "Sparkling;28")
  res8 <- nbinomTest(cds, "Sparkling;28", "Sparkling;32")
  res9 <- nbinomTest(cds, "Sparkling;32", "Sparkling;36")
  res10 <- nbinomTest(cds, "Sparkling;36", "Sparkling;40")
  res11 <- nbinomTest(cds, "Sparkling;40", "Sparkling;44")
  
  sp_pvals <- data.frame(res1$pval, res2$pval, res3$pval, res4$pval, res5$pval, res6$pval, res7$pval, res8$pval, res9$pval, res10$pval, res11$pval)
  sp_fold <- data.frame(res1$foldChange, res2$foldChange, res3$foldChange, res4$foldChange, res5$foldChange, res6$foldChange, res7$foldChange, res8$foldChange, res9$foldChange, res10$foldChange, res11$foldChange)
  rownames(sp_fold) <- rownames(sp_pvals) <- res1$id
  colnames(sp_fold) <- colnames(sp_pvals) <- c("SP0:4", "SP4:8", "SP8:12", "SP12:16", "SP16:20", "SP20:24", "SP24:28", "SP28:32", "SP32:36", "SP36:40", "SP:40:44")
  
  howmanynas <- c()
  for(i in 1:dim(sp_pvals)[1]){
    howmanynas[i] <- length(which(is.na(sp_pvals[i, ]) == T))
  }
  
  sp_pvals <- sp_pvals[which(howmanynas < 11), ]
  sp_fold <- sp_fold[which(howmanynas < 11), ]
  sp_norm <- table[match(rownames(sp_pvals), rownames(table)), ]
  
  assign(paste("sp_pvals", j, sep = ""), sp_pvals)
  assign(paste("sp_fold", j, sep = ""), sp_fold)
  assign(paste("sp_norm", j, sep = ""), sp_norm)
  
  print(paste("Sparkling", j))
  #Trout
  res1 <- nbinomTest(cds, "Trout;0", "Trout;4")
  res2 <- nbinomTest(cds, "Trout;4", "Trout;8")
  res3 <- nbinomTest(cds, "Trout;8", "Trout;12")
  res4 <- nbinomTest(cds, "Trout;12", "Trout;16")
  res5 <- nbinomTest(cds, "Trout;16", "Trout;20")
  res6 <- nbinomTest(cds, "Trout;20", "Trout;24")
  res7 <- nbinomTest(cds, "Trout;24", "Trout;28")
  #missing timepoints
  res11 <- nbinomTest(cds, "Trout;40", "Trout;44")
  
  tb_pvals <- data.frame(res1$pval, res2$pval, res3$pval, res4$pval, res5$pval, res6$pval, res7$pval, res11$pval)
  tb_fold <- data.frame(res1$foldChange, res2$foldChange, res3$foldChange, res4$foldChange, res5$foldChange, res6$foldChange, res7$foldChange, res11$foldChange)
  rownames(tb_fold) <- rownames(tb_pvals) <- res1$id
  colnames(tb_fold) <- colnames(tb_pvals) <- c("TB0:4", "TB4:8", "TB8:12", "TB12:16", "TB16:20", "TB20:24", "TB24:28", "TB:40:44")
  
  howmanynas <- c()
  for(i in 1:dim(tb_pvals)[1]){
    howmanynas[i] <- length(which(is.na(tb_pvals[i, ]) == T))
  }
  
  tb_pvals <- tb_pvals[which(howmanynas < 8), ]
  tb_fold <- tb_fold[which(howmanynas < 8), ]
  tb_norm <- table[match(rownames(tb_pvals), rownames(table)), ]
  
  assign(paste("tb_pvals", j, sep = ""), tb_pvals)
  assign(paste("tb_fold", j, sep = ""), tb_fold)
  assign(paste("tb_norm", j, sep = ""), tb_norm)
  
  print(paste("Trout", j))
  
  #Mendota
  res1 <- nbinomTest(cds, "Mendota;0", "Mendota;4")
  res2 <- nbinomTest(cds, "Mendota;4", "Mendota;8")
  res3 <- nbinomTest(cds, "Mendota;8", "Mendota;12")
  res4 <- nbinomTest(cds, "Mendota;12", "Mendota;16")
  res5 <- nbinomTest(cds, "Mendota;16", "Mendota;20")
  res6 <- nbinomTest(cds, "Mendota;20", "Mendota;24")
  res7 <- nbinomTest(cds, "Mendota;24", "Mendota;28")
  res8 <- nbinomTest(cds, "Mendota;28", "Mendota;32")
  res9 <- nbinomTest(cds, "Mendota;32", "Mendota;36")
  res10 <- nbinomTest(cds, "Mendota;36", "Mendota;40")
  res11 <- nbinomTest(cds, "Mendota;40", "Mendota;44")
  
  me_pvals <- data.frame(res1$pval, res2$pval, res3$pval, res4$pval, res5$pval, res6$pval, res7$pval, res8$pval, res9$pval, res10$pval, res11$pval)
  me_fold <- data.frame(res1$foldChange, res2$foldChange, res3$foldChange, res4$foldChange, res5$foldChange, res6$foldChange, res7$foldChange, res8$foldChange, res9$foldChange, res10$foldChange, res11$foldChange)
  rownames(me_fold) <- rownames(me_pvals) <- res1$id
  colnames(me_fold) <- colnames(me_pvals) <- c("ME0:4", "ME4:8", "ME8:12", "ME12:16", "ME16:20", "ME20:24", "ME24:28", "ME28:32", "ME32:36", "ME36:40", "ME:40:44")
  
  howmanynas <- c()
  for(i in 1:dim(me_pvals)[1]){
    howmanynas[i] <- length(which(is.na(me_pvals[i, ]) == T))
  }
  
  me_pvals <- me_pvals[which(howmanynas < 11), ]
  me_fold <- me_fold[which(howmanynas < 11), ]
  me_norm <- table[match(rownames(me_pvals), rownames(table)), ]
  
  assign(paste("me_pvals", j, sep = ""), me_pvals)
  assign(paste("me_fold", j, sep = ""), me_fold)
  assign(paste("me_norm", j, sep = ""), me_norm)
  
  print(paste("Mendota", j))
}

sp_pvals <- rbind(sp_pvals1, sp_pvals2, sp_pvals3, sp_pvals4, sp_pvals5, sp_pvals6, sp_pvals7, sp_pvals8, sp_pvals9, sp_pvals10, sp_pvals11, sp_pvals12, sp_pvals13, sp_pvals14, sp_pvals15, sp_pvals16, sp_pvals17, sp_pvals18)
sp_fold <- rbind(sp_fold1, sp_fold2, sp_fold3, sp_fold4, sp_fold5, sp_fold6, sp_fold7, sp_fold8, sp_fold9, sp_fold10, sp_fold11, sp_fold12, sp_fold13, sp_fold14, sp_fold15, sp_fold16, sp_fold17, sp_fold18)
sp_norm <- rbind(sp_norm1, sp_norm2, sp_norm3, sp_norm4, sp_norm5, sp_norm6, sp_norm7, sp_norm8, sp_norm9, sp_norm10, sp_norm11, sp_norm12, sp_norm13, sp_norm14, sp_norm15, sp_norm16, sp_norm17, sp_norm18)

tb_pvals <- rbind(tb_pvals1, tb_pvals2, tb_pvals3, tb_pvals4, tb_pvals5, tb_pvals6, tb_pvals7, tb_pvals8, tb_pvals9, tb_pvals10, tb_pvals11, tb_pvals12, tb_pvals13, tb_pvals14, tb_pvals15, tb_pvals16, tb_pvals17, tb_pvals18)
tb_fold <- rbind(tb_fold1, tb_fold2, tb_fold3, tb_fold4, tb_fold5, tb_fold6, tb_fold7, tb_fold8, tb_fold9, tb_fold10, tb_fold11, tb_fold12, tb_fold13, tb_fold14, tb_fold15, tb_fold16, tb_fold17, tb_fold18)
tb_norm <- rbind(tb_norm1, tb_norm2, tb_norm3, tb_norm4, tb_norm5, tb_norm6, tb_norm7, tb_norm8, tb_norm9, tb_norm10, tb_norm11, tb_norm12, tb_norm13, tb_norm14, tb_norm15, tb_norm16, tb_norm17, tb_norm18)

me_pvals <- rbind(me_pvals1, me_pvals2, me_pvals3, me_pvals4, me_pvals5, me_pvals6, me_pvals7, me_pvals8, me_pvals9, me_pvals10, me_pvals11, me_pvals12, me_pvals13, me_pvals14, me_pvals15, me_pvals16, me_pvals17, me_pvals18)
me_fold <- rbind(me_fold1, me_fold2, me_fold3, me_fold4, me_fold5, me_fold6, me_fold7, me_fold8, me_fold9, me_fold10, me_fold11, me_fold12, me_fold13, me_fold14, me_fold15, me_fold16, me_fold17, me_fold18)
me_norm <- rbind(me_norm1, me_norm2, me_norm3, me_norm4, me_norm5, me_norm6, me_norm7, me_norm8, me_norm9, me_norm10, me_norm11, me_norm12, me_norm13, me_norm14, me_norm15, me_norm16, me_norm17, me_norm18)

for(i in 1:dim(tb_pvals)[2]){
  print(length(which(tb_pvals[,i] < 0.05)))
}

write.csv(sp_pvals, file = "D:/GEODES_mapping_summaries/Sparkling_pvalues_2017-06-08.csv")
write.csv(tb_pvals, file = "D:/GEODES_mapping_summaries/TroutBog_pvalues_2017-06-08.csv")
write.csv(me_pvals, file = "D:/GEODES_mapping_summaries/Mendota_pvalues_2017-06-08.csv")

write.csv(sp_fold, file = "D:/GEODES_mapping_summaries/Sparkling_foldchange_2017-06-08.csv")
write.csv(tb_fold, file = "D:/GEODES_mapping_summaries/TroutBog_foldchange_2017-06-08.csv")
write.csv(me_fold, file = "D:/GEODES_mapping_summaries/Mendota_foldchange_2017-06-08.csv")

write.csv(sp_norm, file = "D:/GEODES_mapping_summaries/Sparkling_normalized_counts_2017-06-08.csv")
write.csv(tb_norm, file = "D:/GEODES_mapping_summaries/TroutBog_normalized_counts_2017-06-08.csv")
write.csv(me_norm, file = "D:/GEODES_mapping_summaries/Mendota_normalized_counts_2017-06-08.csv")

