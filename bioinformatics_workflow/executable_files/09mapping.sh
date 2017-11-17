#!/bin/bash
#Map metatranscriptome reads to my pre-indexed database of reference genomes
#Transfer metaT from gluster
#Not splitting the metaTs anymore
cp /mnt/gluster/amlinz/GEODES_nonrRNA/$1 .
cp /mnt/gluster/amlinz/ref_MAGs_SAGs.tar.gz .

#Unzip program and database
tar -xvzf BBMap_36.99.tar.gz
tar -xvf samtools.tar.gz
tar -xvzf ref_MAGs_SAGs.tar.gz
gzip -d $1
name=$(basename $1 | cut -d'.' -f1)
sed -i '/^$/d' $name.fastq

#Run the mapping step
bbmap/bbmap.sh in=$name.fastq out=$name.MAGsSAGs.mapped.sam minid=0.95 trd=T sam=1.3 threads=1 build=1 mappedonly=T -Xmx20g

# I want to store the output as bam. Use samtools to convert.
./samtools/bin/samtools view -b -S -o $name.MAGsSAGs.mapped.bam $name.MAGsSAGs.mapped.sam

#Copy bam file back to gluster
cp $name.MAGsSAGs.mapped.bam /mnt/gluster/amlinz/GEODES_mapping_results/

#Clean up
rm -r bbmap
rm -r ref
rm *.bam
rm *.sam
rm *.fastq
rm *.gz
