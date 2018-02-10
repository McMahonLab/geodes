#!/bin/bash
#Build a re-usable mapping index
cp /mnt/gluster/amlinz/ref_MAGs_SAGs_database.fna .

tar -xvzf BBMap_36.99.tar.gz
bbmap/bbmap.sh ref=ref_MAGs_SAGs_database.fna -Xmx10g

# Make ref/ a tarball and move to gluster
tar czvf ref_MAGs_SAGs.tar.gz ref/

cp ref_MAGs_SAGs.tar.gz /mnt/gluster/amlinz/
rm ref_MAGs_SAGs.tar.gz
rm *database.fna
rm -rf bbmap/
