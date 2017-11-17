#!/bin/bash
awk '{print $1}' GEODES001-nonrRNA.MAGsSAGs.CDS.txt | tail -n +2 > rownames.txt
for file in *.CDS.txt;do awk '{print $7}' $file | tail -n +3 > temp.txt;sample=$(basename $file |cut -d'.' -f1);echo $sample | cat - temp.txt > temp2.txt && mv temp2.txt $file;done
                                                                    
files=$(echo *CDS.txt)
paste rownames.txt $files > GEODES_genes_2017-10-16.txt
