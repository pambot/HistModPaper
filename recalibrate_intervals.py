import sys
import re
import numpy as np
import pandas as pd

filename = sys.argv[1]
rtype = filename.split('_')[0]
ttype = filename.split('_')[1].split('.')[-2]
nr = int(sys.argv[2])
bfile = pd.read_csv(filename, sep='\t', header=None)
tdict = {'g1':'TSS', 'g2':'TTS'}
ttype = tdict[ttype]
drop = []
for n in range(len(bfile.index)):
	b = bfile.loc[n]
	gr = b[3].split(':')
	g1 = int(gr[2])
	g2 = int(gr[3])
	if g2 - g1 < 2000:
		drop.append(n)
	if b[5] == '+':
		gb = g1
		ge = g2
	elif b[5] =='-':
		gb = g2
		ge = g1
	if ttype.lower() == 'tss':
		v = gb
	elif ttype.lower() == 'tts':
		v = ge
	bfile.loc[n, 1] = v - nr
	bfile.loc[n, 2] = v + nr + 1

	
bfile = bfile.drop(drop, axis=0)
bfile.to_csv(rtype.lower().title()+ttype.upper()+sys.argv[3]+'.bed', sep='\t', header=False, index=False)
