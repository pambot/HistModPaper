from __future__ import division
import sys
import glob
import numpy as np
import pandas as pd

# get name stubs
try:
    cell = sys.argv[1]
except IndexError:
    cell = 'Gm12878'

fnames = glob.glob('multiraw/*'+cell+'*.multicov')
fstubs = [fname.replace('multiraw/', '').replace('.multicov', '') for fname in fnames]
fdict = dict(zip(fnames, fstubs))

scounts = pd.read_csv('samcounts.txt', header=None, sep='\t', index_col=0)
loci = pd.read_csv('loci.txt', header=None)
multicov = pd.DataFrame(index=range(len(loci.index)), columns=fstubs)

for fn in fnames:
    fs = fdict[fn]
    mc = pd.read_csv(fn, sep='\t', header=None)
    multicov[fs] = mc[6]*(10**6)/scounts.ix[fs+'.bam', 1]

multicov.index = loci[0].values
multicov.to_csv('multicov/'+cell+'.multicov', sep='\t', header=True, index=True)


