#!/bin/bash
#Build a re-usable mapping index
cp /mnt/gluster/amlinz/nonredundant_database.fna.gz .
gzip -d nonredundant_database.fna.gz
tr -d ' ' < nonredundant_database.fna > nonredundant_database_nospace.fna

#Unzip bbmap and build the index
tar -xvzf BBMap_36.99.tar.gz
bbmap/bbmap.sh ref=nonredundant_database_nospace.fna usemodulo=T -Xmx48g

# Make ref/ a tarball and move to gluster
tar czvf ref.tar.gz ref/

cp ref.tar.gz /mnt/gluster/amlinz/
rm ref.tar.gz
rm nonredundant_database.fna
rm -rf bbmap/
