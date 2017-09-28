#!/bin/bash
mkdir metaG_gffs
cat metagenome.txt | while read line;
  do split -l 5000 -a 4 -d /mnt/gluster/amlinz/metagenome_assemblies/gff/$line.assembled.gff metaG_gffs/$line;
done
#2500 should aim for about 1000 jobs per metagenome assembly, so 6000 total
# Make a list of files to run - only doing a couple to test
ls metaG_gffs > metaG_gffs.txt
