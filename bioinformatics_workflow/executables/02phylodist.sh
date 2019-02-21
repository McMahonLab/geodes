#!/bin/bash
#Classify contigs based on their gene's USEARCH hits provided by JGI
metaG=$(echo $1 | head -c 9)
cp /mnt/gluster/amlinz/metagenome_assemblies/phylogeny/$metaG.assembled.phylodist .
cat $1 | while read line
  do grep $line $metaG.assembled.phylodist;
  done > $1.phylodist

#add header
echo $'locus_tag\thomolog_gene_oid\thomolog_taxon_oid\tpercent_identity\tlineage' | cat - $1.phylodist > temp.phylodist && mv temp.phylodist $1.phylodist

tar xvzf python.tar.gz
#Update the path variable
mkdir home
export PATH=$(pwd)/python/bin:$PATH
export HOME=$(pwd)/home

chmod +x classifyWphylodist_contigs.py

python classifyWphylodist_contigs.py -pd $1.phylodist -pc .70 -hm 3 -conlen 18

rm *phylodist
rm *py
rm -rf python/
rm -r home/
rm python.tar.gz  
