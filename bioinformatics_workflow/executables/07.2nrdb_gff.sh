#!/bin/bash

tar xvzf genometools.tar.gz
export PATH=$(pwd)/genometools/bin:$PATH

gzip -d unprocessed-nrdb.gff.gz
echo '##gff-version 3' | cat - unprocessed-nrdb.gff > temp.gff
gt gff3 -sort yes -tidy -retainids -o sorted_database.gff temp.gff

echo "pFN18A_DNA_transcript	>Cluster X	CDS	1	917	.	+	0	ID=pFN18A-DNA-transcript" > temp.gff
cat sorted_database.gff temp.gff > nonredundant_database.gff
sed -i 's/>//' nonredundant_database.gff

gzip nonredundant_database.gff
mv nonredundant_database.gff.gz /mnt/gluster/amlinz

rm *tar.gz
rm -r genometools
rm *gff
