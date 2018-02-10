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
#genome = 'GEODES006'
taxonFile = genome + '.contig.classification.perc70.minhit3.txt'
productFile = genome + '.assembled.product_names'
inputGff = genome + '.assembled.gff'
outputGff = genome + '.parsed.gff'
outputTable = genome + '.table.txt'

#%%#############################################################################
### Update the inputGff file. Replace ID with 'locus tag' field
### Make a separate table with taxonomy and product name info
################################################################################

# Store the classification file as a dictionary

readme = pandas.read_csv(taxonFile, sep='\t', header = None)
readme.columns = ['Contig', 'TaxString']
taxonDict = readme.set_index('Contig').to_dict()['TaxString']

# Store the product info as a dictionary

product_table = pandas.read_csv(productFile, sep = '\t')
product_table.columns = ['Gene', 'Product', 'KO']
productDict = product_table.set_index('Gene').to_dict()['Product']

# Read in the GFF file
# Each record contains all sequences belonging to the same contig
# For each sequence within the record, replace the ID with the locus_tag

inFile = open(inputGff, 'r')
outFile1 = open(outputGff, 'w')
outFile2 = open(outputTable, 'w')

limit_info = dict(gff_type = ["CDS"],)
for record in GFF.parse(inFile, limit_info=limit_info):
    for seq in record.features:
        seq.id = seq.qualifiers['locus_tag'][0] # this is a list for some reason
        seq.qualifiers['ID'][0] = seq.id + "_" + genome
        del seq.qualifiers['locus_tag']

        taxonomy = taxonDict.get(seq.id[:18])
        product = productDict.get(seq.id)
        outFile2.write(seq.qualifiers['ID'][0]+'\t'+seq.id[:18]+'\t'+str(taxonomy)+'\t'+str(product)+'\n')
    GFF.write([record], outFile1)

inFile.close()
outFile1.close()
outFile2.close()
