#!/bin/bash
#Combine htseq-count results into one table
cp /mnt/gluster/amlinz/GEODES_mapping_summaries/* .

awk '{print $1}' GEODES001_nonrRNA.CDS.out > rownames.txt
for file in *.out;do awk '{print $2}' $file > temp.txt;sample=$(basename $file |cut -d'.' -f1);echo $sample | cat - temp.txt > temp2.txt && mv temp2.txt $file;done
echo "\t" | cat - rownames.txt > temp3.txt && mv temp3.txt rownames.txt
files=$(echo *.out)
paste rownames.txt $files > GEODES_genes_2017-02-22.txt

rm *.out
rm rownames.txt
rm temp*.txt
