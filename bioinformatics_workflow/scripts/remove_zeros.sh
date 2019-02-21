#!/bin/bash

awk 'NR > 1{s=0; for (i=3;i<=NF;i++) s+=$i; if (s!=0)print}' GEODES_ID75_2017-12-27.txt > GEODES_ID75_2017-12-27.readcounts.txt
awk '{print $1}' GEODES_ID75_2017-12-27.readcounts.txt > genes.txt
gzip GEODES_ID75_2017-12-27.readcounts.txt





