#!/bin/bash
#Count mapped reads

#Unzip programs
tar -xvzf subreads.tar.gz
#Pair reads
#./subread-1.5.2-source/bin/utilities/repair -i /mnt/gluster/amlinz/GEODES_mapping_results/$1.mapped.MAGsSAGs.bam -o $1.repaired.bam -T 1
#Count reads
./subread-1.5.2-source/bin/featureCounts -t CDS -g ID --fracOverlap 0.75 -M --fraction --donotsort -a /mnt/gluster/amlinz/ref_MAGs_SAGs_database.gff -o $1.MAGsSAGs.CDS.txt /mnt/gluster/amlinz/GEODES_mapping_results/$1.MAGsSAGs.mapped.bam

#Clean up - don't delete the text file, that is being sent back to the submit node
rm *.tar.gz
rm -r subread-1.5.2-source
