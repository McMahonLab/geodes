#!/bin/bash
tar xvzf python.tar.gz
mkdir home
export PATH=$(pwd)/python/bin:$PATH
export HOME=$(pwd)/home

cp /mnt/gluster/amlinz/ref_MAGs_SAGs_key.txt .
chmod +x genekey.py

python genekey.py

rm GEODES_genes_2017-10-25.txt
rm ref_MAGs_SAGs_key.txt
rm genekey.py
rm -r home
rm -rf python
rm python.tar.gz
