#!/bin/bash 
#Sort metatranscriptomic reads into rRNA and non-rRNA files

#Transfer the fasta file from gluster
cp $1 ./
name=$(basename $1 |cut -d'.' -f1)

#Unzip files
tar -xvf sortmerna-2.1-linux-64.tar.gz

cd sortmerna-2.1-linux-64

#Index the rRNA databases
./indexdb_rna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:\./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db

#Run the sorting program
./sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db --reads ../${name}.fastq  --fastx --aligned ${name}_rRNA --other ${name}_nonrRNA --log -v -m 1 -a 1

#Let's unpack the parameters I've chosen. --ref refers to the databases I've just indexed. --reads says use the fastq file provided. --fastx means I would like a fastq file as output. --aligned is the name for rRNA reads, --other is the name for non-rRNA reads. --log says compute statistics about the run. -v means verbose. I've left the alignment setting at the default, --best. There's a faster setting, but I'm hoping I won't need it. And finally, -m tells sortmerna it can use 3GB of RAM to hold each piece of the fastq file as it processes. Hopefully that will give it a fighting chance.
#Move the output files back to gluster
mv ${name}_rRNA.fastq /mnt/gluster/amlinz/GEODES_rRNA/
mv ${name}_nonrRNA.fastq /mnt/gluster/amlinz/GEODES_nonrRNA/

#Remove files
cd ..
rm ${name}.fastq
rm sortmerna-2.1-linux-64.tar.gz
rm -r sortmerna-2.1-linux-64
