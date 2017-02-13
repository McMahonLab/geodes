#!/bin/bash

cp /mnt/gluster/amlinz/GEODES_metaT/$1.filter-MTF.fastq.gz ./
gzip -d $1.filter-MTF.fastq.gz

# make a folder for each metagenome
mkdir $1

# Split by number of lines
split -l 1000000 $1.filter-MTF.fastq $1/$1

# Rename files to include .fastq extension
for file in $1/*;do mv "$file" "$file.fastq";done

# Grace includes a read file in the folder with the split files, but I don't think that makes sense for using with SortmeRNA
# I'm also not going to zip up the file, because I'd just need to unzip back on gluster
mv $1 /mnt/gluster/amlinz/GEODES_metaT_split/
rm $1.filter-MTF.fastq
