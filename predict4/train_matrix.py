import sys
import pandas as pd

try:
    cell = sys.argv[1]
except IndexError:
    cell = 'Gm12878'

froot = '../predict{0}/train/'

collect = []
for n in range(1, 4):
    df = pd.read_csv(froot.format(n)+cell+'.matrix', sep='\t', header=0, index_col=0)
    df.columns = [c+'_{0}'.format(n) for c in df.columns]
    collect.append(df)

matrix = pd.concat(collect, axis=1)
matrix.to_csv('train/'+cell+'.matrix', header=True, index=True, sep='\t')


