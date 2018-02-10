#!/bin/bash
gzip -d CDS_regions_genekey.txt.gz
while read line; do grep $line CDS_regions_genekey.txt; done < $1 > $1.gene.info.txt

rm $1
rm CDS_regions_genekey.txt
