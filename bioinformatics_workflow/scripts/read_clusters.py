################################################################################
# Make table of info about each gene cluster
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import pandas
import os
import sys

#%%#############################################################################
### Define input files
################################################################################

inCluster = "test.clstr"
outTable1 = "cluster_info.txt"
inKey = "CDS_regions_genekey.txt"

f = open(inCluster, 'r')
s = open(outTable1, 'w')
cluster = None

for line in f:
  if line.startswith('>Cluster'):
      cluster = line.strip('\n')
      continue
  if ('nt' in line):
      split_line = line.split('\t')
      part2 = split_line[1].split('nt, >')
      length = part2[0]
      part3 = part2[1].split('... ')
      name = part3[0]
      threshold = part3[1]
  newline = cluster + '\t' + length +'\t' + name + '\t' + threshold
  print newline
  s.write(newline)

s.close()
  
genekey = pandas.read_csv(inKey, sep = '\t', names = ['Gene_Name', 'Genome', 'Taxonomy', 'Product'])
genekey['Truncated'] = genekey['Gene_Name'].str[0:19]
taxonDict = genekey.set_index('Truncated').to_dict()['Taxonomy']
productDict = genekey.set_index('Truncated').to_dict()['Product']


c = pandas.read_csv(outTable1, sep = '\t', names = ['Cluster', 'Length', 'Name', 'Threshold'])
c["Taxonomy"] = c["Name"].map(taxonDict)
c["Product"] = c["Name"].map(productDict)

pandas.c.to_csv(outTable1)
