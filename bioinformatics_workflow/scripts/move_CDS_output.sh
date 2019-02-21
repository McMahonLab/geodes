#!/bin/bash
mkdir /home/amlinz/CDS_fastas
mv /home/amlinz/CDS*.fna /home/amlinz/CDS_fastas
cat /home/amlinz/CDS_fastas/* /home/amlinz/ref_MAGs_SAGs/fastas/pFN18A_DNA_transcript.fna > CDS_regions.fna
echo "pFN18A_DNA_transcript  internal standard internal standard internal standard" > internalstd.txt
cat /home/amlinz/*.table.txt internalstd.txt > CDS_regions_genekey.txt
gzip CDS_regions.fna
cp CDS_regions.fna.gz /mnt/gluster/amlinz/
cp CDS_regions_genekey.txt /mnt/gluster/amlinz/
