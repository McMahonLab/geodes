#!/bin/bash
#Downsample metagenomes to 10 percent
#Transfer metaG from gluster
cp /mnt/gluster/amlinz/GEODES_metagenomes/$1-metaG.fastq.gz .

#Unzip program
tar -xvzf BBMap_36.99.tar.gz
gzip -d $1-metaG.fastq.gz

#BBmap to sample
bbmap/reformat.sh in=$1-metaG.fastq out=$1-sampled.fastq samplerate=0.1

gzip $1-sampled.fastq

#Copy file back to gluster
cp $1-sampled.fastq.gz /mnt/gluster/amlinz/downsampled_metagenomes/

#Clean up
rm -r bbmap
rm *.fastq
rm *.gz
