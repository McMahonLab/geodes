#!/bin/bash
#Cluster coding regions to get nonredundant genes and make a dummy gff file to go with it

tar xvzf cd-hit.tar.gz

cp /mnt/gluster/amlinz/CDS_regions.fna.gz .
gzip -d CDS_regions.fna.gz

./cd-hit-v4.6.8-2017-0621/cd-hit-est -i CDS_regions.fna -o nonredundant_database.fna -c 0.97 -M 36000 -T 12 -d 50

gzip nonredundant_database.fna
mv nonredundant_database.fna.gz /mnt/gluster/amlinz
mv nonredundant_database.fna.clstr /mnt/gluster/amlinz

rm *tar.gz
rm *fna
rm *clstr
rm *gff
rm *txt
rm -r cd-hit-v4.6.8-2017-0621
