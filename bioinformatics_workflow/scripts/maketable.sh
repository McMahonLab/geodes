#!/bin/bash
awk '{print $1}' GEODES005-sampled.90.CDS.txt | tail -n +2 > rownames.txt
for file in *.CDS.txt;do awk '{print $7}' $file | tail -n +3 > temp.txt;sample=$(basename $file |cut -d'.' -f1);echo $sample | cat - temp.txt > temp2.txt && mv temp2.txt $file;done
                                                                    
files=$(echo *CDS.txt)
paste rownames.txt $files > GEODES_metaG_ID90_2018-03-10.txt

mkdir split_genes

#Remove rows that sum to zero
awk 'NR > 1{s=0; for (i=3;i<=NF;i++) s+=$i; if (s!=0)print}' GEODES_metaG_ID90_2018-03-10.txt > GEODES_metaG_ID90_2018-03-10.readcounts.txt
awk '{print $1}' GEODES_metaG_ID90_2018-03-10.readcounts.txt > genes.txt
split -l 1000 -a 4 -d genes.txt split_genes/genes

#Add column name
head -1 GEODES_metaG_ID90_2018-03-10.txt > colnames.txt
cat colnames.txt GEODES_metaG_ID90_2018-03-10.readcounts.txt > temp.txt
mv temp.txt GEODES_metaG_ID90_2018-03-10.readcounts.txt

#1500 should aim for just under 10000 jobs
# Make a list of files to run - only doing a couple to test
ls split_genes > split_genes.txt

