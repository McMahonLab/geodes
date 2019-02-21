#!/bin/bash
#Install genome tools and python
tar xvzf genometools.tar.gz
tar xvzf python.tar.gz

#Update the path variable
mkdir home
export PATH=$(pwd)/python/bin:$PATH
export HOME=$(pwd)/home
export PATH=$(pwd)/genometools/bin:$PATH

gzip -d $1.gff.gz
gzip -d $1.fna.gz

#Remove CRISPRs
sed -i '/CRISPR/d' $1.gff

#Split up gff file
grep "##\|#!" $1.gff > top_of_file.txt
grep "CDS" $1.gff > bottom_of_file.txt
awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8}' bottom_of_file.txt > part1.txt;
awk -F '\t' '{print $9}' bottom_of_file.txt > part2.txt;
#Get IDs
awk -F ";" '{print $1}' part2.txt > IDs.txt;
#Get start location
awk '{print $4}' bottom_of_file.txt > start.txt;

#look up phylogeny and make that column, too
rows=$(wc -l bottom_of_file.txt | awk '{print $1}')
assignment=$(grep -n $1 algae_phylogeny.txt | awk '{print $2}')
yes $assignment | head -n $rows > classification.txt
awk -F ',' '{print  $6}' classification.txt > species.txt

#Pull out product tags
awk -F 'product=' '{print $2}' part2.txt | awk -F ';' '{print $1}' > products.txt
#Replace blank lines with "None given"
sed -i -e 's/^$/None given/' products.txt
sed -i -e 's/^/product=/' products.txt

#Put the columns back together
paste -d "\t'" part1.txt IDs.txt > cat1.txt;
paste -d "_" cat1.txt start.txt > cat2.txt;
paste -d "_" cat2.txt species.txt > cat3.txt;
paste -d ";" cat3.txt products.txt > newbottom.txt

#Put the sequence region tags back on top
cat top_of_file.txt newbottom.txt > new.gff;

#Sort with genometools
gt gff3 -sort yes -tidy -retainids -o sorted.$1.gff new.gff;
# Extract the coding regions from the fasta file
gt extractfeat -type CDS -seqid no -retainids yes -seqfile $1.fna -matchdescstart sorted.$1.gff >  CDS.$1.fna

#Run the python script to build the key table

chmod +x ref_algae_parsing.py
python ref_algae_parsing.py $1

rm *tar.gz
rm algae_phylogeny.txt
rm *py
rm -r genometools
rm -r home
rm -rf python
rm $1.fna
rm $1.fna.*
rm *.gff
rm bottom_of_file.txt
rm top_of_file.txt
rm part1.txt
rm part2.txt
rm IDs.txt
rm classification.txt
rm products.txt
rm cat*txt
rm species.txt
rm start.txt
rm newbottom.txt
