import sys
import numpy as np
import pandas as pd

cell = sys.argv[1]
rtype = sys.argv[2]

multicov = pd.read_csv('../cage/'+cell+'/'+rtype+'.multicov', sep='\t', header=0, index_col=0)
multicov.columns = [cell+c.split(cell)[-1] for c in multicov.columns]

for ttype in ['tss', 'tts']:
	bed = pd.read_csv(rtype.lower().title()+ttype.upper()+'.bed', sep='\t', header=None)
	for col in multicov.columns:
		nzind = np.nonzero(multicov[col])[0]
		zind = np.where(multicov[col]==0)[0]
		nzlocs = bed.loc[nzind]
		zlocs = bed.loc[zind]
		nzlocs.to_csv('COORDS/hg19/'+col+rtype.title()+ttype.upper()+'Expr.bed', sep='\t', header=False, index=False)
		zlocs.to_csv('COORDS/hg19/'+col+rtype.title()+ttype.upper()+'Zero.bed', sep='\t', header=False, index=False)
