#!/bin/bash
for file in /mnt/gluster/amlinz/filtered/*; do sample=$(basename $file |cut -d'.' -f1); echo $sample;done > samplenames.txt
#head -3 samplenames.txt > test_samplenames.txt; mv test_samplenames.txt samplenames.txt
