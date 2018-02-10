#!/bin/bash
tar xvzf genometools.tar.gz
tar xvzf python.tar.gz

#Update the path variable
mkdir home
export PATH=$(pwd)/python/bin:$PATH
export HOME=$(pwd)/home
export PATH=$(pwd)/genometools/bin:$PATH

cp /mnt/gluster/amlinz/metagenome_assemblies/fastas/GEODES006.assembled.fna.gz .
gzip -d GEODES006.assembled.fna.gz
cp /mnt/gluster/amlinz/phylodist_results/GEODES006.contig.classification.perc70.minhit3.txt .
cp /mnt/gluster/amlinz/metagenome_assemblies/product_names/GEODES006.assembled.product_names .

mv $1 GEODES006.assembled.gff
chmod +x metaG_parsing.py

awk -F'\t' -vOFS='\t' '{gsub("-1", "-", $7); gsub("1", "+", $7); print}' GEODES006.assembled.gff > f1_GEODES006.assembled.gff
echo '##gff-version 3' | cat - f1_GEODES006.assembled.gff > temp && mv temp GEODES006.assembled.gff
#gt gff3 -sort yes -tidy -retainids -force -o GEODES006.assembled.gff f1_GEODES006.assembled.gff

python metaG_parsing.py "GEODES006"

mv GEODES006.parsed.gff $1.parsed.gff
mv GEODES006.table.txt $1.table.txt

gt gff3 -sort yes -tidy -retainids -o sorted.$1.gff $1.parsed.gff #clean up the the gff sorter

# Extract the coding regions from the fasta file
gt extractfeat -type CDS -seqid no -retainids yes -seqfile GEODES006.assembled.fna -matchdescstart sorted.$1.gff >  CDS.$1.fna

rm GEODES006*
rm *gff
rm metaG_parsing.py
rm *tar.gz
rm -r genometools
rm -r home
rm -rf python
