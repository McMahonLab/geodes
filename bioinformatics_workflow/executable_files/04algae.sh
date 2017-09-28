#!/bin/bash
#Rename genes to include relevant info from the gff file, then extract a fasta of only coding regions
tar xvzf genometools.tar.gz
export PATH=$(pwd)/genometools/bin:$PATH

gzip -d *gz

grep "##\|#!" $1.gff > top_of_file.txt
grep -v "##\|#!" $1.gff > bottom_of_file.txt
grep "CDS" bottom_of_file.txt > temp.txt && temp.txt bottom_of_file.txt

awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8}' bottom_of_file.txt > part1.txt;
awk '{print $9}' bottom_of_file.txt > part2.txt;
awk -F ";" '{print $1}' part2.txt > IDs.txt;
awk -F "product=" '{print $2}' part2.txt > gene.txt;
awk '{print $1}' bottom_of_file.txt > genome.txt;
awk '{print $4}' bottom_of_file.txt > start.txt;
#look up phylogeny and make that column, too
rows=$(wc -l bottom_of_file.txt | awk '{print $1}')
assignment=$(grep -n $1 algae_phylogeny.txt | awk '{print $2}')
yes $assignment | head -n $rows > classification.txt
paste -d "\t'" part1.txt IDs.txt > cat1.txt;
paste -d "," cat1.txt start.txt > cat2.txt;
paste -d "," cat2.txt genome.txt > cat3.txt;
paste -d "," cat3.txt classification.txt > cat4.txt; #paste phylogeny to ID
paste -d "," cat4.txt gene.txt > newbottom.txt;
cat top_of_file.txt newbottom.txt > new.gff;
gt gff3 -sort yes -tidy -retainids -o sorted_$1.gff new.gff;

gt extractfeat -type CDS -seqid no -retainids yes -seqfile $1.fna -matchdescstart sorted_$1.gff >  CDS_$1.fna

rm *tar
rm -r genometools
rm *txt
rm *gff

