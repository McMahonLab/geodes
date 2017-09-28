#!/bin/bash

tar xvzf genometools.tar.gz
export PATH=$(pwd)/genometools/bin:$PATH

#no top of file characters
grep "CDS" $1 > bottom_of_file.txt
metaG=$(echo $1 | cut -c1-9)

awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8}' bottom_of_file.txt > part1.txt #split by first 8 columns - metagenome assembly has issue with whitespace in the fields.
awk -F'\t' '{print $9}' bottom_of_file.txt > part2.txt #put the last of column of tags in its own file

awk -F";" '{for(i=1;i<=NF;i++){if ($i ~ /locus_tag=*/){print $i}}}' part2.txt > locust.txt
#These files have the relevant info in locustag instead of id
awk '{print substr($1,11)}' locust.txt > IDS.txt
# from product name file
while read line;
  do grep -wF $line /mnt/gluster/amlinz/metagenome_assemblies/product_names/$metaG.assembled.product_names;
  done < IDS.txt > gene.txt
awk -F'\t' '{print $2}' gene.txt > product.txt

rows=$(wc -l bottom_of_file.txt | awk '{print $1}')
yes $metaG | head -n $rows > genome.txt
#look up phylogeny and make that column, too

while read line; do echo $line | cut -c1-18; done < IDS.txt > contigs.txt
sed -i -e 's/^$/NA/' contigs.txt
while read line;
  do line2=`grep -wF "$line" /mnt/gluster/amlinz/phylodist_results/$metaG.contig.classification.perc70.minhit3.txt`;
  [ ! -z "$line2" ] && echo $line2 || echo "NA   No_classification"
  done < contigs.txt > class.txt
  
awk '{print $2}' class.txt > classification.txt

#add the ID tag back in
sed -i -e 's/^/ID=/' IDS.txt

awk -F'\t' '{print $4}' bottom_of_file.txt > start.txt #get the gene start location
paste -d "\t" part1.txt IDS.txt > cat1.txt #paste first 8 columns to ID tag
paste -d "," cat1.txt start.txt > cat2.txt #paste start location to ID
paste -d "," cat2.txt genome.txt > cat3.txt #paste genome to ID
paste -d "," cat3.txt classification.txt > cat4.txt #paste phylogeny to ID
paste -d "," cat4.txt product.txt > new.gff #paste product name to ID

#Fix the strand issue in the assemblies
awk -F'\t' -vOFS='\t' '{gsub("-1", "-", $7); gsub("1", "+", $7); print}' new.gff > f1.gff

#Add first comment line to the assemblies
echo '##gff-version 3' | cat - f1.gff > temp && mv temp f1.gff

awk -F'\t' -vOFS='\t' '{gsub(";", ":", $9); print}' f1.gff > temp && mv temp f1.gff

gt gff3 -sort yes -tidy -retainids -o sorted_$1.gff f1.gff #clean up the the gff sorter

cp /mnt/gluster/amlinz/metagenome_assemblies/fastas/$metaG.assembled.fna.gz .
gzip -d $metaG.assembled.fna.gz
gt extractfeat -type CDS -seqid no -retainids yes -seqfile $metaG.assembled.fna -matchdescstart sorted_$1.gff >  CDS_$1.fna

rm *tar.gz
rm -r genometools
rm *txt
rm *gff
rm $1
rm $metaG.assembled.*
