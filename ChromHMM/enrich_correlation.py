from __future__ import division
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scipy
import scipy.cluster.hierarchy
import scipy.spatial.distance

# for cell in Gm12878 Helas3 Hepg2 H1hesc Huvec K562 Nhek ; do for ttype in TSS TTS ; do python enrich_correlation.py $cell $ttype ; done ; done

cell = sys.argv[1]
ttype = sys.argv[2]
fname = 'Hist_load/Model_20/'+cell+'_'+ttype+'_enrich.txt'
enrich = pd.read_csv(fname, sep='\t', header=0, index_col=0)
dM = scipy.spatial.distance.pdist(enrich.T, metric='correlation')
dSM = scipy.spatial.distance.squareform(dM)
labels = [c.replace(cell, '').split('.')[0] for c in enrich.columns]


fig = plt.figure(figsize=(15.0, 15.0), dpi=100)
cmap = matplotlib.cm.YlOrRd

ax_dendro = plt.subplot2grid((10,10), (0, 2), rowspan=2, colspan=8)
lM = scipy.cluster.hierarchy.linkage(dSM, method='centroid')
dendro = scipy.cluster.hierarchy.dendrogram(lM, orientation='top', distance_sort='ascending', link_color_func=lambda k: 'k')
leaves = dendro['leaves']
X = dSM[leaves, :]
X = X[:, leaves]
X = 1 - X
relabel = dict(zip(range(len(labels)), labels))
newlabels = [relabel[x] for x in leaves]
ax_dendro.axis('off')

ax_plot = plt.subplot2grid((10,10), (2,2), rowspan=8, colspan=8)
ax_plot.pcolormesh(X, cmap=cmap, vmin=0, vmax=1)
ax_plot.set_xlim([0, len(newlabels)])
ax_plot.set_ylim([0, len(newlabels)])
ax_plot.set_xticks(np.arange(0.5, len(newlabels)+0.5, 1))
ax_plot.set_yticks(np.arange(0.5, len(newlabels)+0.5, 1))
ax_plot.set_xticklabels(['']*len(newlabels))
ax_plot.set_yticklabels(newlabels, fontsize=12, rotation='horizontal')
ax_plot.tick_params(which='both', bottom='off', top='off', left='off', right='off')

ax_color = fig.add_axes([0.15, 0.8, 0.1, 0.06])
norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
cb = matplotlib.colorbar.ColorbarBase(ax_color, cmap=cmap, norm=norm, orientation='horizontal')

ax_color.set_xticks(np.arange(0,10,5))
ax_color.set_xticklabels([0.0, '', '', '', '', 0.5, '', '', '', '', 1.0], fontsize=12)
ax_color.tick_params(which='both', bottom='on', top='off', left='off', right='off')
plt.savefig('figures/EnrichHM'+cell+ttype+'.png', figsize=(15.0, 15.0), dpi=100)





