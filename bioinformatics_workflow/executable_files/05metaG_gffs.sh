#!/bin/bash
tar xvzf genometools.tar.gz
tar xvzf python.tar.gz

#Update the path variable
mkdir home
export PATH=$(pwd)/python/bin:$PATH
export HOME=$(pwd)/home
export PATH=$(pwd)/genometools/bin:$PATH

metaG=$(basename $1 | cut -d'.' -f1)

cp /mnt/gluster/amlinz/metagenome_assemblies/fastas/$metaG.assembled.fna.gz .
gzip -d $metaG.assembled.fna.gz
cp /mnt/gluster/amlinz/metagenome_assemblies/phylogeny/$metaG.assembled.phylodist .
cp /mnt/gluster/amlinz/metagenome_assemblies/product_names/$metaG.assembled.product_names .
cp /mnt/gluster/amlinz/phylodist_results/$metaG.contig.classification.perc70.minhit3.txt .

mv $1 $metaG.assembled.gff
chmod +x metaG_parsing.py

awk -F'\t' -vOFS='\t' '{gsub("-1", "-", $7); gsub("1", "+", $7); print}' $metaG.assembled.gff > f1_$metaG.assembled.gff
echo '##gff-version 3' | cat - f1_$metaG.assembled.gff > temp && mv temp $metaG.assembled.gff

python metaG_parsing.py $metaG

mv $metaG.parsed.gff $1.parsed.gff
mv $metaG.table.txt $1.table.txt

gt gff3 -sort yes -tidy -retainids -o sorted.$1.gff $1.parsed.gff #clean up the the gff sorter

# Extract the coding regions from the fasta file
gt extractfeat -type CDS -seqid no -retainids yes -seqfile $metaG.assembled.fna -matchdescstart sorted.$1.gff >  CDS.$1.fna


rm *minhit3.txt
rm *product_names
rm $metaG.assembled.fna
rm *gff
rm metaG_parsing.py
rm *tar.gz
rm -r genometools
rm -r home
rm -rf python
