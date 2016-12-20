import sys
import glob
import numpy as np
import pandas as pd

try:
    cell = sys.argv[1]
except IndexError:
    cell = 'Hepg2'

print 'Starting', cell

categories = ['Expr', 'Coding', 'Specificity', 'Pap']
loci = open('loci.txt', 'r').read().splitlines()
loci = pd.DataFrame(loci, columns=['loci'])
loci = pd.concat([loci, pd.DataFrame(columns=categories)])
loci.fillna(0, inplace=True)

expr = pd.read_csv('expr.matrix', sep='\t', header=0)

mc = pd.read_csv('multicov/'+cell+'.multicov', sep='\t', header=0, index_col=0)
mc_sum = np.log(mc.sum(axis=1).values + 1) 
pap = [1 if v else 0 for v in mc[[c for c in mc.columns if 'Pap' in c]].sum(axis=1).values]

loci['Multicov'] = mc_sum
loci['Pap'] = pap
loci['Coding'] = [1 if item[0]=='C' else 0 for item in loci['loci']]
loci['Expr'] = expr[cell]
loci['Specificity'] = expr.sum(axis=1).values*1.0/len(expr.columns)

loci = loci[['loci'] + categories + ['Multicov']]
loci.to_csv('train/'+cell+'.labels', sep='\t', header=True, index=False)



