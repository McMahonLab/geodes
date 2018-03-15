#!/bin/bash

tar xvzf python.tar.gz
#Update the path variable
mkdir home
export PATH=$(pwd)/python/bin:$(pwd)/samtools/bin:$(pwd)/bwa:$PATH
export HOME=$(pwd)/home

# Unzip assembly data files
tar xvzf for_bin_classification.tar.gz

# Make all scripts executable
chmod +x *.sh
chmod +x *.py

# Make list of contigs
grep '>' $1.fna > $1.contigs

# Make gff file
./makeBinGFF.sh $1.contigs for_bin_classification/

# Get bin name with no extension
name=$(basename $1 .fna)

# Make phylodist file
python makeBinPhylodist.py $name.gff for_bin_classification

# Make COG file
python makeBinCOGS.py $name.gff for_bin_classification

# Filter COGs to include only phylogeny marker genes
python filterPhyloCOGs.py $name.cog.txt

# Get the consensus phylogeny
python classifyWphylodistWcutoffs.py $name.phylodist.subphylocog.txt

# Simplify output name
mv *classonly.txt $name.perc70.minhit3.classonly.txt

# Clean up
rm *.tar.gz
rm *.py
rm *.sh
rm *.tsv
rm *.gff
rm *.fna
rm *.fna.contigs
rm *cog.txt
rm *phylodist.txt
rm *minhit3.txt
rm -r for_bin_classification/
rm -rf python/
rm -r home/
