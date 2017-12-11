#!/bin/bash

mv /home/amlinz/CDS*.fna /home/amlinz/CDS_fastas
cat /home/amlinz/CDS_fastas/* /home/amlinz/ref_MAGs_SAGs/fastas/pFN18A_DNA_transcript.fna > /mnt/gluster/amlinz/CDS_regions.fna
cat /home/amlinz/*.table.txt "pFN18A_DNA_transcript\tinternal standard\tinternal standard\tinternal standard" > /mnt/gluster/amlinz/CDS_regions_genekey.txt
gzip /mnt/gluster/amlinz/CDS_regions.fna
