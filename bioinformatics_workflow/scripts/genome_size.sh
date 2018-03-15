#!/bin/bash

cp /mnt/gluster/amlinz/GEODES_nonrRNA/$1-nonrRNA.fastq.gz .
gzip -d $1-nonrRNA.fastq.gz

sed -i '/^$/d' $1-nonrRNA.fastq  
awk 'NR%4 ==2 {print $1}' $1-nonrRNA.fastq > bases.txt
wc -m bases.txt

rm bases.txt
rm $1-nonrRNA.fastq


