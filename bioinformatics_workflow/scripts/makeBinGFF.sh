#!/bin/bash

# This program makes a gff file for the individual bin from the assembly gff
# Usage: makeBinGFF contiglistFile path2assemblies

# Read in filename
filename=$1
# Get path to assemblies
assemPath=$2
# Put together name for assembly gff file
gffFilename=`echo "$filename" | cut -f1 -d- `
#gffFilename=${gffFilename##*/}
# Make bin GFF file name
gffOut="${filename%.contigs}".gff

# Grep each line of the filename (minus first character) from the gff file 
#    into new bin gff file
while read line
do
grep "${line:1:${#line}}" $assemPath/$gffFilename.assembled.gff >> $gffOut
done < $filename
