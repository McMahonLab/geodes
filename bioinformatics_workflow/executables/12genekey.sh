#!/bin/bash

gzip -d *readcounts.txt.gz
cp /mnt/gluster/amlinz/GEODES_ID90_genekey.txt .

tar xzf python.tar.gz
export PATH=$(pwd)/python/bin:$(pwd)/samtools/bin:$PATH
export HOME=$(pwd)/home
chmod +x genekey.py

python genekey.py Mendota_ID90_2018-01-07.readcounts.txt Mendota_ID90_genekey.csv
python genekey.py Trout_ID90_2018-01-07.readcounts.txt Trout_ID90_genekey.csv
python genekey.py Sparkling_ID90_2018-01-07.readcounts.txt Sparkling_ID90_genekey.csv

cp *csv /mnt/gluster/amlinz

rm *txt
rm *csv
rm *py
rm python.tar.gz
rm -rf python/
rm -rf home/
