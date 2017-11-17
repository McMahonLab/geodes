#!/bin/bash

#mkdir metaG_gffs
cd metaG_gffs
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES005.assembled.gff GEODES005.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES006.assembled.gff GEODES006.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES057.assembled.gff GEODES057.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES058.assembled.gff GEODES058.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES117.assembled.gff GEODES117.
split -500000 /mnt/gluster/amlinz/metagenome_assemblies/gff/GEODES118.assembled.gff GEODES118.
cd ..
# make list of files to run
ls metaG_gffs > metaG_gffs.txt

for file in ref_MAGs_SAGs/fastas/*;do name=$(basename $file | cut -d'.' -f1); echo $name; done > refMAGs_SAGs_list.txt
for file in refseq_algae/fastas/*;do name=$(basename $file .fna.gz); echo $name; done > algae_list.txt

