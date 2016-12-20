import sys
import glob
import numpy as np
import pandas as pd

folder = 'multicov'
files = glob.glob(folder+'/*.multicov')
loci = pd.read_csv('../predict2/loci.txt', header=None)

cells = [f.replace('multicov', '')[1:-1] for f in files if '.log.' not in f]
expr = pd.DataFrame(columns=cells, index=loci[0])

# extract info from multicov, log transform, resave
for cell in cells:
    mc = pd.read_csv('multicov/'+cell+'.multicov', sep='\t', header=0, index_col=0)
    expr[cell] = [1 if val else 0 for val in mc.sum(axis=1).values]

expr.to_csv('expr.matrix', sep='\t', header=True, index=True)


