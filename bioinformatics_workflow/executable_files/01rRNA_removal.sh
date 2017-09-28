#!/bin/bash
#Sort metatranscriptomic reads into rRNA and non-rRNA files

#Transfer the fasta file from gluster
name=$(basename $1 |cut -d'.' -f1)

#Unzip files
tar -xvf sortmerna-2.1-linux-64.tar.gz
gzip -d $name.fastq.gz
cd sortmerna-2.1-linux-64

#Index the rRNA databases
./indexdb_rna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:\./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db

#Run the sorting program
./sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db --reads ../$name.fastq  --fastx --aligned ../$name-rRNA --other ../$name-nonrRNA --log -v -m 1 -a 1

cd ..
gzip *RNA.fastq
mkdir $1-rRNA
mkdir $1-nonrRNA
mv *-rRNA.fastq.gz $1-rRNA
mv *nonrRNA.fastq.gz $1-nonrRNA

#Remove files
rm sortmerna-2.1-linux-64.tar.gz
rm -r sortmerna-2.1-linux-64
