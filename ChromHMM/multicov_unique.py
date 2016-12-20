import sys
import numpy as np
import pandas as pd

rtype = sys.argv[1]
cells = ['Gm12878', 'Helas3', 'Hepg2', 'H1hesc', 'Huvec', 'K562', 'Nhek']

specdf = pd.DataFrame()
for cell in cells:
	mc = pd.read_csv('../cage/'+cell+'/'+rtype+'.multicov', sep='\t', header=0, index_col=0)
	specdf[cell] = mc.sum(axis=1)


for cell in cells:
	spec = []
	for n in range(len(specdf.index)):
		s = specdf.loc[n]
		# test for cell-specific expression
		if s[cell] and s[cell] == s.sum():
			spec.append(n)
	spec = np.array(spec)
	nonedf = specdf.sum(axis=1)
	none = np.where(nonedf==0)[0]
	for ttype in ['TSS', 'TTS']:
		bed = pd.read_csv(rtype+ttype+'.bed', sep='\t', header=None)
		bed.loc[spec].to_csv(cell+rtype+ttype+'Spec.bed', sep='\t', header=None, index=None)


for ttype in ['TSS', 'TTS']:
	bed = pd.read_csv(rtype+ttype+'.bed', sep='\t', header=None)
	bed.loc[none].to_csv(rtype+ttype+'None.bed', sep='\t', header=None, index=None)




