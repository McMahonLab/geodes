"""makeBinPhylodist.py : make cog file for each bin, requires matching names for COG file and phylodist files, preferably made with makeBinXXX.py scripts"""

__author__ = "Sarah Stevens"
__email__ = "sstevens2@wisc.edu"

import pandas as pd, sys

def usage():
	print("Usage: filterPhyloCOGs.py  COGsFile")

if len(sys.argv) != 2:
	usage()
	sys.exit(2)

cogname=sys.argv[1]
prefix=cogname.split('.cog.txt')[0]
phyloname=prefix+'.phylodist.txt'

matchup=pd.read_table('Phylosift2COGs.tsv',header=0) # hard coded in path2file, add as argument?


cogs=pd.read_table(cogname,header=0)
phylodist=pd.read_table(phyloname,header=0)
phylocogs=cogs[cogs.cog_id.isin(matchup.COG_num)]
phylocogs.to_csv(prefix+'.phylocog.txt',sep='\t', index=False)

phylodist[phylodist.locus_tag.isin(phylocogs.locus_tag)].to_csv(prefix+'.phylodist.subphylocog.txt',sep='\t', index=False)
