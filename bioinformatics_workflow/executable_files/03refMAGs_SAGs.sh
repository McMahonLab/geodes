#!/bin/bash
#Install genome tools and python
tar xvzf genometools.tar.gz
tar xvzf python.tar.gz

#Update the path variable
mkdir home
export PATH=$(pwd)/python/bin:$PATH
export HOME=$(pwd)/home
export PATH=$(pwd)/genometools/bin:$PATH

#Remove CRISPRs
grep -v "CRISPR" $1.gff > temp.gff && mv temp.gff $1.gff
#Setup and run python script
chmod +x ref_MAGs_SAGs_parsing.py
python ref_MAGs_SAGs_parsing.py $1

#Remove all the duplicate gff-version lines
grep -v "##gff-version 3" $1.parsed.gff > int1.gff
# Move sequence region lines to top of file 
grep "##sequence-region" int1.gff > sequence_regions.gff
grep -v "##sequence-region" int1.gff > not_sequence_regions.gff

# Put it all back together and clean up
echo '##gff-version 3' | cat - sequence_regions.gff not_sequence_regions.gff > $1.fixed.gff 
gt gff3 -sort yes -tidy -retainids -o sorted.$1.gff $1.fixed.gff #clean up the the gff sorter

# Extract the coding regions from the fasta file
gt extractfeat -type CDS -seqid no -retainids yes -seqfile $1.fna -matchdescstart sorted.$1.gff >  CDS.$1.fna

rm *tar.gz
rm *csv
rm *py
rm -r genometools
rm -r home
rm -rf python
rm $1.fna
rm $1.fna.*
rm *.gff
