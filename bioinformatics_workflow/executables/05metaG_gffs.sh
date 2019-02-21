#!/bin/bash
tar xvzf genometools.tar.gz
tar xvzf python.tar.gz

#Update the path variable
mkdir home
export PATH=$(pwd)/python/bin:$PATH
export HOME=$(pwd)/home
export PATH=$(pwd)/genometools/bin:$PATH

metaG=$(echo $1 | cut -c1-9)
tar xvzf $metaG.datafiles2.tar.gz
gzip -d $metaG.assembled.fna.gz

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
