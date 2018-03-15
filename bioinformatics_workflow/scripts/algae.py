###############################################################################
# CodeTitle.py
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
#genome = '2582580615'
taxonFile = 'algae_phylogeny.txt'
inputGff = genome + '.gff'
outputGff = genome + '.parsed.gff'
outputTable = genome + '.table.txt'

#%%#############################################################################
### Update the inputGff file. Replace ID with 'locus tag' field
### Make a separate table with taxonomy and product name info
################################################################################

# Store the classification file as a dictionary

readme = pandas.read_table(taxonFile, names = ["Genome", "Taxonomy"])
taxonDict = readme.set_index('Genome').to_dict()['Taxonomy']

# Read in the GFF file
# Each record contains all sequences belonging to the same contig
# For each sequence within the record, replace the ID with the locus_tag

inFile = open(inputGff, 'r')
outFile1 = open(outputGff, 'w')
outFile2 = open(outputTable, 'w')

for record in GFF.parse(inFile):
    for seq in record.features:
        seq.qualifiers['ID'][0] = seq.qualifiers['ID'][0] + genome
        if 'product' in seq.qualifiers.keys():
            product = seq.qualifiers['product'][0]
            del seq.qualifiers['product']
        else:
            product = 'None given'

        taxonomy = taxonDict[genome]
        outFile2.write(seq.id+'\t'+genome+'\t'+str(taxonomy)+'\t'+product+'\n')
    GFF.write([record], outFile1)

inFile.close()
outFile1.close()
outFile2.close()
