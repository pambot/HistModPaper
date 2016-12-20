import os
import sys
import glob
import re
import pandas as pd
from itertools import groupby
from operator import itemgetter

try:
    cell = sys.argv[1]
except IndexError:
    cell = 'Gm12878'

pm = 1000
print 'Starting', cell

matrix = open('loci.txt', 'r').read().splitlines()
fnames = glob.glob('summary/*'+cell+'*.summary')

collect = []
for hist in ['H2az', 'H3k0?9me3', 'H3k27ac', 'H3k27me3', 'H3k36me3', 'H3k0?4me1', 'H3k0?4me2', 'H3k0?4me3', 'H3k79me2', 'H3k0?9ac', 'H4k20me1']:
    matrix = pd.DataFrame(index=range(len(matrix)))
    hname = hist.replace('0?', '')
    for ts in ['TSS', 'TTS']:
        for n in range(pm*2/200):
            matrix[hname+'_'+ts+'_'+str(n)] = [0]*len(matrix.index)
        
        summaries = [fn for fn in fnames if re.search(hist, fn) and ts in fn]
        for sm in summaries:
            summary = pd.read_csv(sm, header=None, sep=' ')
            summary.replace(inplace=True, to_replace='n/a', value='0')
            matrix.loc[:, [hname+'_'+ts+'_'+str(n) for n in range(pm*2/200)]] = summary.loc[:, 3:13].values
    
    collect.append(matrix)
    print 'Done', cell, hname

matrix = pd.read_csv('loci.txt', header=None, names=['loci'])
collect.insert(0, matrix)
merged = pd.concat(collect, axis=1)

fmask = '../predict1/train/'+cell+'.matrix'
if os.path.isfile(fmask):
    mask = pd.read_csv(fmask, sep='\t', header=0)
    merged[mask==0] = 0
    print 'Done masking'
else:
    print 'No masking applied. File not found.'

merged.to_csv('train/'+cell+'.matrix', index=False, sep='\t')



