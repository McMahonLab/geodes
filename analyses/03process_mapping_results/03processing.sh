#!/bin/bash
#Make a count table of mapped reads

#Transfer data from gluster
cp /mnt/gluster/amlinz/GEODES_mapping_concat/$1.all.sam .
cp /mnt/gluster/amlinz/python.tar.gz .

#Unzip programs
tar xzf python.tar.gz
tar xvf samtools.tar.gz
gzip -d mapping_database.gff.gz

#Update the path variable
mkdir home
export PATH=$(pwd)/python/bin:$PATH
export HOME=$(pwd)/home

#Manipulate the output
samtools view -b -S $1.all.sam > $1.bam 
samtools sort -o $1.sorted $1.bam 
samtools index $1.sorted.bam

#Count reads
python ./python/bin/htseq-count -f bam -r pos -s no -a 0 -t CDS -i locus_tag -m intersection-strict -o $1.CDS.sam $1.all.sam mapping_database.gff > $1.CDS.out

#Save output to gluster. Using cp/rm instead of mv so it will overwrite old output in gluster.
cp $1.CDS.out /mnt/gluster/amlinz/GEODES_mapping_summaries/
rm $1.CDS.out

#Clean up
rm -rf python
rm -rf samtools
rm *.sam
#rm *.bam
#rm *.bai
rm mapping_database*
rm *.tar.gz
