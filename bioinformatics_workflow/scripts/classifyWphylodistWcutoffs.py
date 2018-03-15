import pandas as pd, sys, os

"""classifyWPhylodist.py  for each taxonomic level takes the classificaiton with the most hits,
	removes any not matching hits for the next leve, also returns the number and avg pid for those hits"""

__author__ = "Sarah Stevens"
__email__ = "sstevens2@wisc.edu"

def usage():
	print("Usage: classifyWPhylodistWcutoffs.py  phylodistFile")

if len(sys.argv) != 2:
	usage()
	sys.exit(2)

# Read in args
inname=sys.argv[1]
phylodist = pd.read_table(inname)
perc_co = .70
hit_min = 3
	

# Setting up output
outname=os.path.splitext(inname)[0]+'.perc70.minhit3.txt'
outname2=os.path.splitext(inname)[0]+'.perc70.minhit3.classonly.txt' # output file without the numbers in it
output=open(outname,'w')
output2=open(outname2,'w')
output.write(inname+'\t')
output2.write(inname+'\t')
if phylodist.empty: # if there are no hits in phylodist file
	print('Phylodist empty!  This maybe because it was subsetted into phylocogs and there were none.')
	output.write('NO CLASSIFICATION BASED ON GIVEN PHYLODIST\n')
	output2.write('NO CLASSIFICATION BASED ON GIVEN PHYLODIST\n')
	sys.exit(0)

phylodist['contig']=phylodist['locus_tag'].str[:17]
taxon_ranks=['Kingdom','Phylum','Class','Order','Family','Genus','Species','Taxon_name']
phylodist[taxon_ranks]=phylodist['lineage'].apply(lambda x: pd.Series(x.split(';')))
totGenes=len(phylodist)
genome_classification=list()


for rank in taxon_ranks:
	# Getting taxon, at that level, hit the most and the number of hits
	counts=phylodist.groupby(rank).size()
	max_ct_value=counts.max()
	perc=phylodist.groupby(rank).size()/float(totGenes) # percent of times each classification (added float for python2)
	max_pc_value=perc.max() # number of percentage the classification with the most hits had
	max_pc_taxon=perc.idxmax(axis=1) # classificaiton hit the most at this rank
	averages=phylodist.groupby(rank).mean().reset_index()[[rank,'percent_identity']] #getting average pid for best classifcation
	max_pc_avg=averages.loc[averages[rank]==max_pc_taxon, 'percent_identity'].iloc[0]
	result=(max_pc_taxon, round(max_pc_value, 2), round(max_pc_avg, 2), max_ct_value)
	genome_classification.append(result)
	rank_co_dict = {'Kingdom':.20,'Phylum':.45,'Class':.49,'Order':.53,'Family':.61,'Genus':.70,'Species':.90,'Taxon_name':.97}
	rank_co = rank_co_dict[rank]
	if (result[0]>.70 and result[3]>hit_min and result[1]>rank_co):
		output.write('{}({},{},{});'.format(result[0],result[1],result[2],result[3]))
		output2.write(result[0]+';')
	else:
		if rank == 'Kingdom':
			if result[3]<=hit_min:
				output.write('NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST')
				output2.write('NO CLASSIFICATION DUE TO FEW HITS IN PHYLODIST')
			else: 
				output.write('NO CLASSIFICATION DUE TO LOW PERCENT MATCHING')
				output2.write('NO CLASSIFICATION DUE TO LOW PERCENT MATCHING')
		output.write('\n')
		output2.write('\n')
		output.close()
		output2.close()
		sys.exit(0)
	phylodist=phylodist[phylodist[rank]==max_pc_taxon] # removing any hits which don't match the classification at this level

output.write('\n')
output.close()
output2.write('\n')
output2.close()
