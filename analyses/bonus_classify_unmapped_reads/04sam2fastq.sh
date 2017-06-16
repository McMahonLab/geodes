#!/bin/bash
#Convert sam to fastq and split
tar xvf samtools.tar.gz

#Copy file from gluster
cp /mnt/gluster/amlinz/GEODES_mapping_concat/$1.all.sam .

#Extract fastq file of unmapped reads
samtools view -S -b -f 4 $1.all.sam > $1.unaligned.bam 
samtools bam2fq $1.unaligned.bam > $1.unaligned.fastq

mkdir $1

# Split by number of lines
split -l 100000 $1.unaligned.fastq $1/$1

# Rename files to include .fastq extension
for file in $1/*;do mv $file $file.unaligned.fastq;done

cp $1/* /mnt/gluster/amlinz/GEODES_kraken_split/
rm $1.unaligned.fastq
rm *sam
rm *bam
rm -r $1
rm -r samtools
rm samtools.tar.gz
