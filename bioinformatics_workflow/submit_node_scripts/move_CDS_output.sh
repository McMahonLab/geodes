#!/bin/bash

mv /home/amlinz/CDS*.fna /home/amlinz/CDS_fastas
cat /home/amlinz/CDS_fastas/* > /mnt/gluster/amlinz/CDS_regions.fna
gzip /mnt/gluster/amlinz/CDS_regions.fna

