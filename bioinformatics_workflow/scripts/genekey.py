#%%#############################################################################
### Import packages
################################################################################

import pandas
import os
import sys

#%%#############################################################################
### Define input files
################################################################################

readcounts = sys.argv[1]
genekey  = 'GEODES_ID90_genekey.txt'
outputkey = sys.argv[2]

#%%#############################################################################
### Keep only genes present in the read count table
################################################################################

# Open files

df = pandas.read_table(readcounts, sep = " ")
inKey = pandas.read_table(genekey, names = ["Gene", "Genome", "Taxonomy", "Product"])
outKey = open(outputkey, 'w')

# Remove genes in the key file that are no longer in the table
inKey = inKey[inKey['Gene'].isin(df['Geneid'])]
inKey.to_csv(outKey, index=False)
