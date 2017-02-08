#!/bin/bash
#Concatenate sortmerna output

cp /mnt/gluster/amlinz/GEODES_nonrRNA/$1??_nonrRNA.fastq ./
cat $1??_nonrRNA.fastq > $1_nonrRNA.fastq
nonrRNAcount=$(wc -l $1_nonrRNA.fastq)
gzip $1_nonrRNA.fastq
mv $1_nonrRNA.fastq.gz /mnt/gluster/amlinz/GEODES_nonrRNA_concat/
rm *_nonrRNA.fastq

cp /mnt/gluster/amlinz/GEODES_rRNA/$1*??_rRNA.fastq ./
cat $1??_rRNA.fastq > $1_rRNA.fastq
rRNAcount=$(wc -l $1_rRNA.fastq)
gzip $1_rRNA.fastq
mv $1_rRNA.fastq.gz /mnt/gluster/amlinz/GEODES_rRNA_concat/
rm *_rRNA.fastq

echo ${nonrRNAcount},${rRNAcount} > $1_rRNA_results.txt
