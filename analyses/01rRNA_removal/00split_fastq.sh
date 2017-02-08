#!/bin/bash

cp /mnt/gluster/amlinz/GEODES_metaT/GEODES001.filter-MTF.fastq.gz ./
gzip -d GEODES001.filter-MTF.fastq.gz

# make a folder for each metagenome
code=`echo "GEODES001.filter-MTF.fastq" | cut -d'.' -f1`
mkdir ${code}

# Split by number of lines
split -l 1000000 GEODES001.filter-MTF.fastq ${code}/${code}

# Rename files to include .fastq extension
for file in ${code}/*;do mv "$file" "$file.fastq";done

# Grace includes a read file in the folder with the split files, but I don't think that makes sense for using with SortmeRNA
# I'm also not going to zip up the file, because I'd just need to unzip back on gluster
mv ${code} /mnt/gluster/amlinz/GEODES_metaT_split/
rm GEODES001.filter-MTF.fastq
