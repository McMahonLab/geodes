#!/bin/bash
#Rename genes to include relevant info from the gff file, then extract a fasta of only coding regions
tar xvzf genometools.tar.gz
export PATH=$(pwd)/genometools/bin:$PATH

grep "##\|#!" $1.gff > top_of_file.txt
grep -v "##\|#!" $1.gff > bottom_of_file.txt #take headers off the top of the file
grep "CDS" bottom_of_file.txt > temp.txt && mv temp.txt bottom_of_file.txt

awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8}' bottom_of_file.txt > part1.txt #split by first 8 columns
awk '{print $9}' bottom_of_file.txt > part2.txt #put the last of column of tags in its own file
awk -F ";" '{print $1}' part2.txt > IDs.txt #take out just the ID tag
awk -F "product=" '{print $2}' part2.txt > gene.txt #get the product name
rows=$(wc -l bottom_of_file.txt | awk '{print $1}')
yes $1 | head -n $rows > genome.txt
#look up phylogeny and make that column, too
assignment=$(grep -n $1 Readme.csv | awk -F "," '{print $3,$4,$5,$6,$7,$8}' OFS=",")
yes $assignment | head -n $rows > classification.txt
awk '{print $4}' bottom_of_file.txt > start.txt #get the gene start location
paste -d "\t'" part1.txt IDs.txt > cat1.txt #paste first 8 columns to ID tag
paste -d "," cat1.txt start.txt > cat2.txt #paste start location to ID
paste -d "," cat2.txt genome.txt > cat3.txt #paste genome to ID
paste -d "," cat3.txt classification.txt > cat4.txt #paste phylogeny to ID
paste -d "," cat4.txt gene.txt > newbottom.txt #paste product name to ID
cat top_of_file.txt newbottom.txt > new.gff #put the top of the file back on
gt gff3 -sort yes -tidy -retainids -o sorted_$1.gff new.gff #clean up the the gff sorter


gt extractfeat -type CDS -seqid no -retainids yes -seqfile $1.fna -matchdescstart sorted_$1.gff >  CDS_$1.fna

rm *tar.gz
rm -r genometools
rm *txt
rm *gff

