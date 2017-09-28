#!/bin/bash
mv GEODES*gz /home/amlinz/GEODES_metaT_split/
cd /home/amlinz/GEODES_metaT_split
ls *.tar.gz |xargs -n1 tar -xvzf
rm *.tar.gz
cd /home/amlinz/
gzip /home/amlinz/GEODES_metaT_split/*/*fastq
mv /home/amlinz/GEODES_metaT_split/*/*fastq.gz /home/amlinz/GEODES_metaT_split/
rmdir /home/amlinz/GEODES_metaT_split/*
ls /home/amlinz/GEODES_metaT_split/ > path2splitfastqs.txt
