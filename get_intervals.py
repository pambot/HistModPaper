# get the TSS from all the GENCODE transcripts
# eliminate duplicates
import pandas as pd
import sys

fname = sys.argv[1]
pm = int(sys.argv[2])

with open(fname, 'r') as f:
	tscripts = [l.split('\t') for l in f.read().splitlines()]

bed_pos = [[t[0], t[1], t[2], t[5]] for t in tscripts]
bed = pd.DataFrame(bed_pos, columns=['chrom', 'chromStart', 'chromEnd', 'strand'])

tss_bins = []
tts_bins = []
for n in range(len(bed.index)):
	# tss of transcript
	t = bed.iloc[n]
	if t[0]=='chrM':
		continue
	tss = int(t[1])
	# name is index:chrom:chromStart:chromEnd:strand
	name = str(n)+':'+t[0]+':'+t['chromStart']+':'+t['chromEnd']+':'+t[3]
	tss_bins.append([t[0], tss-pm, tss+1+pm, name, '0', t[3]])
	# tts of transcript
	tts = int(t[2])
	name = str(n)+':'+t[0]+':'+t['chromStart']+':'+t['chromEnd']+':'+t[3]
	tts_bins.append([t[0], tts-pm, tts+1+pm, name, '0', t[3]])

def write_bins(bin_var, fwname):
	bins = pd.DataFrame(bin_var, columns=None)#.drop_duplicates()
	bins.to_csv(fwname+'.bed', sep='\t', index=None, header=False)

write_bins(tss_bins, fname.split('.')[0]+'_tss')
write_bins(tts_bins, fname.split('.')[0]+'_tts')

