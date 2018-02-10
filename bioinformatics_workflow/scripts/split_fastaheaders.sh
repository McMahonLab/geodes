#!/bin/bash
mkdir split_fastaheaders
split -l 1500 -a 4 -d /mnt/gluster/amlinz/fasta_headers_nospace.txt split_fastaheaders/fastaheaders

#1500 should aim for just under 10000 jobs
# Make a list of files to run - only doing a couple to test
ls split_fastaheaders > splitfastaheaders.txt
