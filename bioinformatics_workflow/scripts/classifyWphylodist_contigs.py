import pandas as pd, sys, os, argparse

"""classifyWPhylodist.py  classifys contigs, for each taxonomic level it takes the classificaiton with the most hits,
	removes any not matching hits for the next level, also returns the number and avg pid for those hits"""

__author__ = "Sarah Stevens"
__email__ = "sstevens2@wisc.edu"


# Read in args function
def parseArgs():
	parser = argparse.ArgumentParser(description='classifyWPhylodist_contig.py: classifies each contig from an IMG metagenome annotated phylodist file.')
	parser.add_argument('--phylodist','-pd' , action='store', dest='inname', type=str, required=True, metavar='PhylodistFile', help="Phylodist file from IMG annotation of metagenome")
	parser.add_argument('--percent_cutoff','-pc', action='store', dest='perc_co', default=.70, type=float, help='Minimum percentage of hits that much match at each taxon rank to use classification. DEFAULT: .70')
	parser.add_argument('--hit_minimum','-hm', action='store', dest='hit_min', default=3, type=int, help='Minium number of hits that must match at each taxon rank to use classification. DEFAULT: 3')
	parser.add_argument('--printVals', '-pv', action='store_true', dest='wvals', default=False, help='Include this flag if you would like it to print the percent matching, AAI, and num hits for each match.')
	parser.add_argument('--contignamelen','-conlen', action="store", dest='conlen', default=17, type=int, help='Num of characters from locus_tag to take pull for contig name. DEFAULT: 17')
	parser.add_argument('--alwaysPrintStop','-aps', action='store_true', dest='printStop', default=False, help='Include this flag if you would like it to always print the code, by default only prints stopcode if not classified to kingdom.')
	args=parser.parse_args()
	return args.inname, args.perc_co, args.hit_min, args.wvals, args.conlen, args.printStop

# Read in args
inname, perc_co, hit_min, wvals, conlen, printStop  = parseArgs()

phylodist = pd.read_table(inname)

# Setting contig column
phylodist['contig']=phylodist['locus_tag'].str[:conlen]

# Splitting taxonomy out into multipule cols and fixing the variation in number of cols
phylodist['numTaxons'] = phylodist['lineage'].str.count(';')
# check that 7 and 8 are the only #'s of occurrences
assert (phylodist[(phylodist.numTaxons != 8) & (phylodist.numTaxons != 7)].empty == True), \
		"lineage column has too many or two few ';' seps"
# new df for the 7 occurrences of ';', which is 8 values
eig_phylodist = phylodist[(phylodist.numTaxons == 7)]
if eig_phylodist.empty != True:
	eig_taxon_ranks=['Kingdom','Phylum','Class','Order','Family','Genus','Species','Taxon_name']
	eig_phylodist[eig_taxon_ranks]=eig_phylodist['lineage'].apply(lambda x: pd.Series(x.split(';')))
	eig_phylodist['substrain']=''
# new df for the 8 occurrences of ';', which is nine values
nin_phylodist = phylodist[(phylodist.numTaxons == 8)]
if nin_phylodist.empty != True:
	nin_taxon_ranks=['Kingdom','Phylum','Class','Order','Family','Genus','Species','Taxon_name','substrain']
	nin_phylodist[nin_taxon_ranks]=nin_phylodist['lineage'].apply(lambda x: pd.Series(x.split(';')))
# putting them back into one df or assigning the one to the phylodist df
if (eig_phylodist.empty != True) & (nin_phylodist.empty !=True):
	phylodist = pd.concat([eig_phylodist, nin_phylodist])
elif (eig_phylodist.empty == True) & (nin_phylodist.empty !=True):
	phylodist = nin_phylodist
elif (eig_phylodist.empty != True) & (nin_phylodist.empty == True):
	phylodist = eig_phylodist
else:
	print("You don't have any 8 or 9 taxon genomes")
	sys.exit(2)

# Setting up output
percname=format(float(perc_co), '.2f')[-2:]
outname=os.path.splitext(inname)[0]+'.contig.classification.perc{0}.minhit{1}.txt'.format(percname,hit_min)
if wvals:
	outname=os.path.splitext(inname)[0]+'.contig.classification.perc{0}.minhit{1}.wvals.txt'.format(percname,hit_min)
output=open(outname,'w')
#output.write(inname+'\n')
if phylodist.empty: # if there are no hits in phylodist file
	print('Phylodist empty!  This maybe because it was subsetted into phylocogs and there were none.')
	output.write('NO CLASSIFICATION BASED ON GIVEN PHYLODIST\n')
	sys.exit(0)

# Values for percent identity 
rank_co_dict = {'Kingdom':20,'Phylum':45,'Class':49,'Order':53,'Family':61,'Genus':70,'Species':90,'Taxon_name':97}


for contig in phylodist['contig'].unique():
	tempphylodist=phylodist[phylodist['contig'] == contig]
	contig_classification_str=""
	output.write(contig+'\t')
	totGenes=len(tempphylodist)
	for rank in eig_taxon_ranks:  # this will never get down to the extra last eukaryotic taxon that is sometimes in eukaryotic classifications
		# Getting taxon, at that level, hit the most and the number of hits
		perc=tempphylodist.groupby(rank).size()/float(totGenes) # percent of times each classification, had to add forced float for python2
		counts=tempphylodist.groupby(rank).size()
		max_ct_value=counts.max() # number of hits the classification with the most hits had
		max_ct_perc=perc.max() # percentage of hits the classification with the most hits had
		max_ct_taxon=perc.idxmax(axis=1) # classificaiton hit the most at this rank
		averages=tempphylodist.groupby(rank).mean().reset_index()[[rank,'percent_identity']] #getting average pid for best classifcation
		max_ct_avg=averages.loc[averages[rank]==max_ct_taxon, 'percent_identity'].iloc[0]
		result=(max_ct_taxon, round(max_ct_perc, 2), round(max_ct_avg, 2), max_ct_value)
		if wvals:
			contig_classification_str='{}({},{},{});'.format(result[0],result[1],result[2], result[3])
		else:
			contig_classification_str='{};'.format(result[0])
		tempphylodist=tempphylodist[tempphylodist[rank]==max_ct_taxon] # removing any hits which don't match the classification at this level
		rank_co = rank_co_dict[rank]
		if (max_ct_perc>perc_co and max_ct_value>hit_min and max_ct_avg>rank_co):
			output.write(contig_classification_str)
		else:
			if rank == 'Kingdom':
				stopcode = 'NO CLASSIFICATION '
				if result[3]<=hit_min:
					stopcode += 'MH' # *M*in *H*its too low aka 'NO CLASSIFICATION DUE TO TOO FEW HITS LEFT ON CONTIG AT THIS LEVEL'
				elif max_ct_avg>rank_co:
					stopcode += 'LP' # *L*ow *P*ID aka 'NO CLASSIFICATION DUE TO TOO LOW PID FOR CLASSICICATION AT THIS LEVEL'
				elif max_ct_perc>perc_co:
					stopcode += 'PM' # *P*ercent *M*atching too low aka NO CLASSIFICATION DUE TO TOO LOW PERCENT MATCHING FOR CLASSICICATION AT THIS LEVEL'
				output.write(stopcode)
			elif printStop: 
				stopcode = ''
				if result[3]<=hit_min:
					stopcode += 'MH' # *M*in *H*its too low aka 'NO CLASSIFICATION DUE TO TOO FEW HITS LEFT ON CONTIG AT THIS LEVEL'
				elif max_ct_avg>rank_co:
					stopcode += 'LP' # *L*ow *P*ID aka 'NO CLASSIFICATION DUE TO TOO LOW PID FOR CLASSICICATION AT THIS LEVEL'
				elif max_ct_perc>perc_co:
					stopcode += 'PM' # *P*ercent *M*atching too low aka NO CLASSIFICATION DUE TO TOO LOW PERCENT MATCHING FOR CLASSICICATION AT THIS LEVEL'
				output.write(stopcode)
			break
	output.write('\n')

output.close()
