#!/bin/bash
#Map metatranscriptome reads to my database of reference genomes

#Transfer metaT from gluster
cp $1 .
name=$(basename $1 |cut -d'.' -f1)

#Unzip program and database
tar xvf bwa.tar.gz
gzip -d mapping_database.fna.gz

#Index the reference database
./bwa-0.7.12/bwa index mapping_database.fna

#Run the mapping step - using 5 processors
./bwa-0.7.12/bwa mem -t 3 mapping_database.fna ${name}.fastq > ${name}.sam
cp ${name}.sam /mnt/gluster/amlinz/GEODES_mapping_split/

#Clean up
rm -rf bwa
rm ${name}.fastq
rm *.sam
rm mapping_database*
rm *.tar.gz
