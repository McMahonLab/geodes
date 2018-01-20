#!/bin/bash
#Count mapped reads
cp /mnt/gluster/amlinz/GEODES_mapping_results/$1.90.mapped.bam .
gzip -d nonredundant_database.gff.gz

#Unzip programs
tar -xvzf subreads.tar.gz
#Pair reads
#Count reads
./subread-1.5.2-source/bin/featureCounts -t CDS -g ID --fracOverlap 0.75 -M --fraction --donotsort -a nonredundant_database.gff -o $1.90.CDS.txt $1.90.mapped.bam

#Clean up - don't delete the text file, that is being sent back to the submit node
rm *.tar.gz
rm -r subread-1.5.2-source
rm *bam
rm *gff
rm *.txt.summary
