#!/bin/bash
#Make a dummy gff file to go with the nonredundant database
#start=`expr $1 - 1499`
#sed -n "$start","$1"p /mnt/gluster/amlinz/fasta_headers_nospace.txt > newfile.txt
#echo $1 > vartest.txt
#mv $1 newfile.txt
#Link cluster 
#maxsize=$(sort -nrk1,1 /mnt/gluster/amlinz/nonredundant_database.fna.clstr  | head -1 | cut -f1)
#maxsize=$(($maxsize + 1))
maxsize=300000
touch endpoint.txt
cat $1 | while read line; do var=$(echo $line | cut -c1-19); grep -B $maxsize $var /mnt/gluster/amlinz/nonredundant_database.fna.clstr > splitfile.clstr; cluster=`grep "Cluster" splitfile.clstr | tail -1`; echo $cluster; length=`tail -1 splitfile.clstr | awk '{print $2}' | sed 's/[^0-9]*//g'`; echo $length >> endpoint.txt;done > clusters.txt

#already, now the good part. Make a couple of dummy columns and put the whole shebang together.
# I need:
# first line is version of gff
# genome name - something like "NR_database"
# source - cluster number
# type - CDS
# start - 1
# stop - my calculated endpoints
# strand and frame info - ., +, 0
# tags - use the fasta header as ID

rows=$(wc -l $1 | awk '{print $1}')
yes "NR_gene_database" | head -n $rows > genome.txt
yes "CDS" | head -n $rows > type.txt
yes "1" | head -n $rows > start.txt
yes "." | head -n $rows > info1.txt
yes "+" | head -n $rows > info2.txt
yes "0" | head -n $rows > info3.txt
yes "ID" | head -n $rows > ID.txt

paste -d "=" ID.txt $1 > tags.txt
paste genome.txt clusters.txt type.txt start.txt endpoint.txt info1.txt info2.txt info3.txt tags.txt > $1-nrdb.gff
 
rm *txt
rm $1
