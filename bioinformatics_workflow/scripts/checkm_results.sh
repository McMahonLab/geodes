#!/bin/bash

cat *-contigs.txt > GEODES_binned_contigs.txt
echo -e "bin\ttaxonomy\tsize\tcompleteness\tcontamination" | cat - *-checkm.txt > GEODES_checkm_results.txt
