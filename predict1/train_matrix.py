import sys
import pandas as pd
from itertools import groupby
from operator import itemgetter

try:
    cell = sys.argv[1]
except IndexError:
    cell = 'Gm12878'

pm = 1000
print 'Starting', cell

loci = open('loci.txt', 'r').read().splitlines()

histones = ['H2az', 'H3k27ac', 'H3k27me3', 'H3k36me3', 'H3k4me1', 'H3k4me2', 'H3k4me3', 'H3k79me2', 'H3k9ac', 'H3k9me3', 'H4k20me1']
columns = []
for hist in histones:
    for ts in ['TSS', 'TTS']:
        for n in range(10):
            columns.append('{0}_{1}_{2}'.format(hist, ts, n))

matrix = pd.DataFrame(index=range(len(loci)), columns=['loci']+columns)
matrix['loci'] = loci

# convert loci to +/- 1kbp around locus bed intervals
def getInterval(loc, ts):
	l = loc.split(':')
	num, chrom, n1, n2, strand = l
	if (ts=='TSS' and strand=='+') or (ts=='TTS' and strand=='-'):
		return [chrom, int(n1)-pm, int(n1)+pm+1, num]
	elif (ts=='TSS' and strand=='-') or (ts=='TTS' and strand=='+'):
		return [chrom, int(n2)-pm, int(n2)+pm+1, num]

# fill with binary data
bed = {'TSS':None, 'TTS':None}
for ts in ['TSS', 'TTS']:
    bed[ts] = pd.DataFrame(list(matrix['loci'].map(lambda x: getInterval(x, ts)).values), columns=['chrom', 'start', 'end', 'id'])
    bed[ts]['start'] = (bed[ts]['start']/200).round()*200
    bed[ts]['end'] = (bed[ts]['end']/200).round()*200
    bed[ts].index = bed[ts]['id']

matrix.index = bed['TSS'].index

for chrom in ['chr'+n for n in map(str, range(1,23))+['X','Y']]:
    binary = pd.read_csv('binary/'+cell+'_'+chrom+'_binary.txt', skiprows=0, header=1, sep='\t')
    binary.index = [i*200 for i in binary.index]
    for ts in ['TSS', 'TTS']:
        chrom_bed = bed[ts][bed[ts]['chrom']==chrom]
        for tid in chrom_bed.index:
            r = chrom_bed.ix[tid]
            region = pd.DataFrame(binary.loc[r.start:(r.start+pm*2-200)])
            region.index = range(pm*2/200)
            for hist in histones:
                for n in range(pm*2/200):
                    col = '{0}_{1}_{2}'.format(hist, ts, n)
                    matrix.ix[tid, col] = region.ix[n, hist]
    print 'Done', cell, chrom

matrix.to_csv('train/'+cell+'.matrix', index=False, sep='\t')


