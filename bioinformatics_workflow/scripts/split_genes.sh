#!/bin/bash
mkdir split_genes

awk 'NR > 1{s=0; for (i=3;i<=NF;i++) s+=$i; if (s!=0)print}' GEODES_ID75_2017-12-27.txt > GEODES_ID75_2017-12-27.readcounts.txt
awk '{print $1}' GEODES_ID75_2017-12-27.readcounts.txt > genes.txt
split -l 1000 -a 4 -d genes.txt split_genes/genes

#1500 should aim for just under 10000 jobs
# Make a list of files to run - only doing a couple to test
ls split_genes > split_genes.txt

