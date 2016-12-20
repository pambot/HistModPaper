import sys
import pandas as pd

rtype = sys.argv[1]
if rtype not in ['lncRNA', 'coding']:
	print 'lncRNA or coding'
	sys.exit()

df1 = pd.read_csv(rtype+'_tss_filt.bed', sep='\t', header=None, names=['chrom', 'cStart', 'cEnd', 'name', 'score', 'strand'], dtype='object')
df2 = pd.read_csv(rtype+'_tts_filt.bed', sep='\t', header=None, names=['chrom', 'cStart', 'cEnd', 'name', 'score', 'strand'], dtype='object')

ind1 = [int(it.split(':')[0]) for it in df1['name']]
ind2 = [int(it.split(':')[0]) for it in df2['name']]
dict1 = dict(zip(ind1, df1.index))
dict2 = dict(zip(ind2, df2.index))
ind = set([int(it) for it in set(ind1) & set(ind2)])

drop_df1 = []
for n in ind1:
	if n not in ind:
		drop_df1.append(dict1[n])

drop_df2 = []
for n in ind2:
	if n not in ind:
		drop_df2.append(dict2[n])	

df1 = df1.drop(df1.index[drop_df1])
df2 = df2.drop(df2.index[drop_df2])

assert len(df1) == len(df2)

df1.to_csv(rtype+'_tss.bed', sep='\t', header=False, index=False)
df2.to_csv(rtype+'_tts.bed', sep='\t', header=False, index=False)
