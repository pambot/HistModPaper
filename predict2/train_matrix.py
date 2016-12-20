import sys
import pandas as pd
from itertools import groupby
from operator import itemgetter

cell = sys.argv[1]
nstates = 27
pm = 1000
print 'Starting', cell

loci = open('loci.txt', 'r').read().splitlines()
loci = pd.DataFrame(loci, columns=['loci'])

# populate train matrix
for n in range(1, nstates+1):
	loci['E'+str(n)+'_TSS'] = [0]*len(loci.index)
	loci['E'+str(n)+'_TTS'] = [0]*len(loci.index)

# convert loci to +/- 1kbp around locus bed intervals
def getInterval(loc, ts):
	l = loc.split(':')
	num, chrom, n1, n2, strand = l
	if (ts=='TSS' and strand=='+') or (ts=='TTS' and strand=='-'):
		return [chrom, int(n1)-pm, int(n1)+pm+1, num]
	elif (ts=='TSS' and strand=='-') or (ts=='TTS' and strand=='+'):
		return [chrom, int(n2)-pm, int(n2)+pm+1, num]


# perform intersection with posterior files
for ts in ['TSS', 'TTS']:
	locs = list(loci['loci'].map(lambda x: getInterval(x, ts)).values)
	bed = pd.DataFrame(locs, columns=['chrom', 'start', 'end', 'id'])
	bed['start'] = (bed['start']/200).round()*200
	bed['end'] = (bed['end']/200).round()*200
	loci['id'] = bed['id']
	for chrom in ['chr'+n for n in map(str, range(1,23))+['X','Y']]:
		intervals = bed[bed['chrom']==chrom]
		posterior = pd.read_csv('POSTERIOR/'+cell+'_'+str(nstates)+'_'+chrom+'_posterior.txt', skiprows=0, header=1, sep='\t')
		posterior.index = [i*200 for i in posterior.index]
		posterior.columns = [c+'_'+ts for c in posterior.columns]
		for n in intervals.index.values:
			v = intervals.loc[n]
			scores = pd.DataFrame(posterior.loc[v.start:v.end].sum(axis=0)).T
			loci.ix[n, scores.columns] = scores
			if intervals['id'].loc[n] != loci['id'].loc[n]:
				print 'Error: unmatched IDs'
		print 'Done', chrom, ts


loci = loci.drop('id', axis=1)
loci.to_csv('train/'+cell+'.matrix', index=False, sep='\t')


