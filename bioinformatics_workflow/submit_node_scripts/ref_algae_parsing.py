###############################################################################
# ref_algae_parsing.py
# Copyright (c) 2017, Joshua J Hamilton and Alex Linz
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Make table of info from gff files
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from BCBio import GFF # Python package for working with GFF files
import pandas
import os
import sys

#%%#############################################################################
### Define input files
################################################################################

genome = sys.argv[1]
#genome = "GCF_000149405.2_ASM14940v2_genomic"
taxonFile = 'algae_phylogeny.txt'
inputGff = "sorted." + genome + '.gff'
outputGff = genome + '.parsed.gff'
outputTable = genome + '.table.txt'

#%%#############################################################################
### Update the inputGff file. Replace ID with ID + species name (otherwise genes are just numbered and you can't tell which genome they came from)
### Make a separate table with taxonomy and product name info
################################################################################

# Store the classification file as a dictionary

readme = pandas.read_table(taxonFile, names = ["Genome", "Taxonomy"])
taxonDict = readme.set_index('Genome').to_dict()['Taxonomy']
taxonomy  = taxonDict[genome]

# Read in the GFF file
# Extract product information
# Rename genes including species name because they are just numbered

inFile = open(inputGff, 'r')
outFile2 = open(outputTable, 'w')

limit_info = dict(gff_type = ["CDS"],)
for record in GFF.parse(inFile, limit_info = limit_info):
    for seq in record.features:
	seq.id = seq.qualifiers['ID'][0]
        product = seq.qualifiers['product'][0]
        outFile2.write(seq.id+'\t'+genome+'\t'+str(taxonomy)+'\t'+product+'\n')

inFile.close()
outFile2.close()
