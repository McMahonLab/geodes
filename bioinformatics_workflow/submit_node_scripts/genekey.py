#%%#############################################################################
### Import packages
################################################################################

import pandas
import os
import sys

#%%#############################################################################
### Define input files
################################################################################

readcounts = 'GEODES_genes_2017-10-25.txt'
genekey  = 'ref_MAGs_SAGs_key.txt'
outputcounts = 'GEODES2refMAGsSAGs_readcounts.txt'
outputkey = 'GEODES2refMAGsSAGS_geneinfo.txt'

#%%#############################################################################
### Remove rows of all 0s from the table of read counts and sort rows alphabetically
### Remove the same genes from the gene key and sort alphabetically
################################################################################

# Open files

df = pandas.read_table(readcounts)
inKey = pandas.read_table(genekey, header=None)
outTable = open(outputcounts, 'w')
outKey = open(outputkey, 'w')

# Remove rows that sum to 0
df = df.loc[(df.sum(axis=1) != 0),]

# Sort rows by gene name
df = df.sort_values(by = 'Geneid')

df.to_csv(outputcounts, index=False)

# Remove genes in the key file that are no longer in the table
# Or maybe keep only genes that are still in the table
inKey = inKey[inKey[0].isin(df['Geneid'])]
inKey = inKey.sort_values(by = 0)
inKey.to_csv(outKey, index=False)
