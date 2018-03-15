#!/bin/bash
tar -xvzf BBMap_36.99.tar.gz

cp /mnt/gluster/amlinz/$1.datafiles2.tar.gz .
tar -xvzf $1.datafiles2.tar.gz
gzip -d $1.assembled.fna.gz

./bbmap/reformat.sh in=$1.assembled.fna out=$1-filtered.assembled.fna minlength=1000

gzip $1-filtered.assembled.fna

cp $1-filtered.assembled.fna.gz /mnt/gluster/amlinz/filtered_assemblies/

rm *gz
rm *fna
rm *txt
rm *product_names
rm -r bbmap
