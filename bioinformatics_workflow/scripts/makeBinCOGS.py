# coding: utf-8

"""makeBinPhylodist.py : make cog file for each bin"""

__author__ = "Sarah Stevens"
__email__ = "sstevens2@wisc.edu"

import pandas as pd, sys

def usage():
	print("Usage: makeBinCOGs.py  gffFile path2assemblies")

if len(sys.argv) != 3:
	usage()
	sys.exit(2)

gffname=sys.argv[1]
assemName=gffname.split('-')[0]
path2assem=sys.argv[2]
cogname=path2assem+'/'+assemName+'.assembled.COG'

gff=pd.read_table(gffname, header=None, names=['seqname','source','features','start','end','score','strand','frame','attribute'])
#pulling out the locus tag to its own column
gff['locus_tag']=gff['attribute'].str.split('locus_tag=').str.get(1).str.slice(start=0,stop=-1)
cog=pd.read_table(cogname, header=None, names=['gene_id','cog_id','percent_identity','align_length','query_start','query_end','subj_start','subj_end','evalue','bit_score'])
## Merging and writing out to file, did inner instead of left since some of the genes don't have a hit in the cog file
gff[['locus_tag']].merge(cog, how='inner', left_on='locus_tag', right_on='gene_id').to_csv(gffname.split('.gff')[0]+'.cog.txt',sep='\t', index=False)
