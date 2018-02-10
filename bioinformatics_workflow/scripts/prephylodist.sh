#!/bin/bash

mkdir contig_lists
cat metagenome.txt | while read line;
  do awk '{print substr($1,1,18)}' /mnt/gluster/amlinz/metagenome_assemblies/phylogeny/$line.assembled.phylodist | sort | uniq >   /home/amlinz/$line-contigs.txt;
  split -l 1500 -a 4 -d $line-contigs.txt contig_lists/$line-contigs;
done

# Make a list of files to run - only doing a couple to test
ls contig_lists > metaG_contigs.txt
