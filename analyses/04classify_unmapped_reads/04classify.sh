#!/bin/bash
#Classify unmapped reads using kraken

#Copy files from gluster
cp $1 .
name=$(basename $1 |cut -d'.' -f1)
cp /mnt/gluster/amlinz/kraken.tar.gz .

#Unzip programs and files
tar zxvf  kraken.tar.gz

#Run kraken
cd kraken-0.10.5-beta/kraken_scripts
./kraken --threads 3 --preload --db ../../minikraken_20141208/ ../../$name.unaligned.fastq > $name.output

cp $name.output /mnt/gluster/amlinz/GEODES_kraken_results/

cd ../..
rm -r kraken-0.10.5-beta
rm -r minikraken_20141208
rm $name.unaligned.fastq
rm kraken.tar.gz
