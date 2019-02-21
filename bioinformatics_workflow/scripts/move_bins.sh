#!/bin/bash

#while read line; do tar xvzf /mnt/gluster/amlinz/$line-binning.tar.gz;
#	for file in $line-binning/*fasta; do
#		name=$(basename "$file" .fasta);
#		mv $file metaG_bins/$name.fna;
#		done;
#	rm -r $line-binning;
#	done < metaG_samples.txt

for file in metaG_bins/*; do
	name=$(basename "$file" .fna);
	echo $name;
	done > bins_to_classify.txt

head -3 bins_to_classify.txt > testbins_to_classify.txt

