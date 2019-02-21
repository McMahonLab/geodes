#!/bin/bash

for file in refseq_algae/fastas/*;do name=$(basename $file .fna.gz); echo $name; done > algae_list.txt

#For testing, uncomment the following lines:
#head -2 algae_list.txt > temp.txt
#mv temp.txt algae_list.txt
