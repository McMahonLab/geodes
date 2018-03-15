# coding: utf-8

"""makeBinPhylodist.py : make phyldist file for each bin"""

__author__ = "Sarah Stevens"
__email__ = "sstevens2@wisc.edu"

import pandas as pd, sys

def usage():
	print("Usage: makeBinPhylodist.py  gffFile path2assemblies")

if len(sys.argv) != 3:
	usage()
	sys.exit(2)

gffname=sys.argv[1]
assemName=gffname.split('-')[0]
path2assem=sys.argv[2]
phylodistname=path2assem+'/'+assemName+'.assembled.phylodist'

gff=pd.read_table(gffname, header=None, names=['seqname','source','features','start','end','score','strand','frame','attribute'])
#pulling out the locus tag to its own column
gff['locus_tag']=gff['attribute'].str.split('locus_tag=').str.get(1).str.slice(start=0,stop=-1)
phylodist=pd.read_table(phylodistname, header=None, names=['seq_id','homolog_gene_oid','homolog_taxon_oid','percent_identity','lineage'])
## Merging and writing out to file, did inner instead of left since some of the genes don't have a hit in the phylodist file
gff[['locus_tag']].merge(phylodist, how='inner', left_on='locus_tag', right_on='seq_id').to_csv(gffname.split('.gff')[0]+'.phylodist.txt',sep='\t', index=False)
