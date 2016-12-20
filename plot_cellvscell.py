# change environment
import sys
sys.path.insert(1, '/ifs/home/pw801/bin/python')

# load modules
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from itertools import izip

try:
    folder, fname = sys.argv[1:3]
except (IndexError, ValueError):
    folder, fname = 'predict1/', 'rgrSpecValidR.heatmap'

matrix = pd.read_csv(folder + 'results/'+fname, sep='\t', header=0, index_col=0)

def show_values(pc, fmt='%.2f', **kw):
	pc.update_scalarmappable()
	ax = pc.get_axes()
	for p, color, value in izip(pc.get_paths(), pc.get_facecolors(), pc.get_array()):
		x, y = p.vertices.mean(0)+0.1
		if np.all(color[:3] > 0.2):
			color = '#000000'
		else:
			color = '#FFFFFF'
		ax.text(x, y, fmt % value, ha='center', va='center', color=color, **kw)
	return

fig, ax = plt.subplots(figsize=(6,6))
m = ax.pcolormesh(matrix.values, cmap=plt.get_cmap('YlOrRd'), vmin=0.1, vmax=0.9)
show_values(m, fontsize=14, fontweight='bold')

ax.invert_yaxis()
ax.set_yticks(np.arange(0, len(matrix.index))+0.5)
ax.set_xticks(np.arange(0, len(matrix.columns))+0.5)
ax.set_yticklabels(matrix.index, fontsize=16, color='k')
ax.set_xticklabels(matrix.columns, fontsize=16, color='k', rotation='vertical')
ax.tick_params(direction='out', pad=10)
ax.set_frame_on(True)

for tk1, tk2 in zip(ax.xaxis.get_major_ticks(), ax.yaxis.get_major_ticks()):
	tk1.tick1On, tk2.tick1On, tk1.tick2On, tk2.tick2On = [False]*4

plt.savefig(folder + 'figures/'+fname.split('.')[0]+'.png', bbox_inches='tight')


