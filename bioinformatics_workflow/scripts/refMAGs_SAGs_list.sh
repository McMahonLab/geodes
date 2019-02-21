#!/bin/bash

for file in ref_MAGs_SAGs/fastas/*;do name=$(basename $file | cut -d'.' -f1); echo $name; done > refMAGs_SAGs_list.txt

#For testing, uncomment the following lines:
#head -2 refMAGs_SAGs_list.txt > testrefMAGs_SAGs_list.txt
