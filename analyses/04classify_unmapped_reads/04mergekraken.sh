#!/bin/bash
#Concatenate and summarize kraken outputs

cp /mnt/gluster/amlinz/GEODES_kraken_results/$1* .
cp /mnt/gluster/amlinz/kraken.tar.gz .

#Unzip programs and files
tar zxvf  kraken.tar.gz

cat $1* > $1.output
cd kraken-0.10.5-beta/kraken_scripts

./kraken-report --db ../../minikraken_20141208/ ../../$1.output > ../../$1.report
cd ../..

#Remove spaces in the kraken report
awk '{gsub(" ","",$0)}1' $1.report > temp.txt && mv temp.txt $1.report

#Copy the output file to gluster
cp $1.report /mnt/gluster/amlinz/GEODES_kraken_concat/

#Cleanup
rm -r kraken-0.10.5-beta
rm -r minikraken_20141208
rm $1*.output
rm $1.report
rm kraken.tar.gz

