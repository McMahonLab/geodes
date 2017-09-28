#!/bin/bash
#Map metatranscriptome reads to my pre-indexed database of reference genomes
#Transfer metaT from gluster
#Not splitting the metaTs anymore
cp /mnt/gluster/amlinz/GEODES_nonrRNA/$1 .
cp /mnt/gluster/amlinz/ref.tar.gz .

#Unzip program and database
tar -xvzf BBMap_36.99.tar.gz
tar -xvf samtools.tar.gz
tar -xvzf ref.tar.gz
gzip -d $1
name=$(basename $1 | cut -d'.' -f1)
sed -i '/^$/d' $name.fastq

#Run the mapping step
bbmap/bbmap.sh in=$name.fastq out=$name.mapped.sam minid=0.8 trd=T sam=1.3 threads=1 build=1 usemodulo=T mappedonly=T -Xmx30g

# I want to store the output as bam. Use samtools to convert.
samtools view -b -S -o $name.mapped.bam $name.mapped.sam

#Copy bam file back to gluster
cp $name.mapped.bam /mnt/gluster/amlinz/GEODES_mapping_results/

#Clean up
rm -r bbmap
rm -r ref
rm *.bam
rm *.sam
rm *.fastq
rm *.gz
