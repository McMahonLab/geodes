#!/bin/bash

tar xvzf genometools.tar.gz
export PATH=$(pwd)/genometools/bin:$PATH

echo '##gff-version 3' | cat - /mnt/gluster/amlinz/semiprocessed-nrdb.gff > temp.gff
gt gff3 -sort yes -tidy -retainids -o nonredundant_database.gff temp.gff
gzip nonredundant_database.gff
mv nonredundant_database.gff.gz /mnt/gluster/amlinz

rm *tar.gz
rm -r genometools
rm *gff
