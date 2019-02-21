#!/bin/bash

for RNA in /mnt/gluster/amlinz/GEODES_nonrRNA/*;
	do sample=$(basename $RNA .fastq.gz);
	gzip -d $RNA;
	length=$(grep "@HISEQ" /mnt/gluster/amlinz/GEODES_nonrRNA/$sample.fastq | wc -l) 
	echo $sample $length;
	gzip /mnt/gluster/amlinz/GEODES_nonrRNA/$sample.fastq;
	done > MT_size.txt
