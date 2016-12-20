# load modules
import sys
import glob
import numpy as np
import pandas as pd
import scipy.stats
from statsmodels.sandbox.stats.multicomp import multipletests
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cPickle as pickle

plt.style.use('ggplot')

try:
    condition, test = sys.argv[1:3]
except (IndexError, ValueError):
    condition, test = 'Coding', 'S2'

cells = ['Gm12878', 'H1hesc', 'Helas3', 'Hepg2', 'Huvec', 'K562', 'Nhek']

folder = {
1: 'predict1/results/',
2: 'predict2/results/',
3: 'predict3/results/',
}

def confidence_interval(data, conf=0.95):
	a = 1.0*np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * scipy.stats.t._ppf((1+conf)/2., n-1)
	return m, h

expr_int = ['0.1', '2.0', '4.0', '6.0', '20.0']
sc_types = ['auc', 'fscore']
sc_labels = ['AUC', 'F-score']

for i, sc in enumerate(sc_types):
    collect_m = {1:[], 2:[], 3:[]}
    collect_c = {1:[], 2:[], 3:[]}
    thresh_bins = []
    for thresh1, thresh2 in zip(expr_int, expr_int[1:]):
        thresh_bins.append(thresh1+'-'+thresh2)
        tukey_collect1 = []
        tukey_collect2 = []
        temp_m = {1:[], 2:[], 3:[]}
        temp_c = {1:[], 2:[], 3:[]}
        for n in range(1, 4):
            for cell in cells:
                adj_scores = pickle.load(
                    open(folder[n]+'scores'+cell+condition+'_'+thresh1+'-'+thresh2+'.pkl', 'rb'))
                mean_sc, pm_sc = confidence_interval(adj_scores[sc], conf=0.95)
                temp_m[n].append(mean_sc)
            collect_m[n].append(np.mean(temp_m[n]))
            collect_c[n].append(np.std(temp_m[n]))
            
            for v, k in enumerate(range(1, 4)):
                tukey_collect1.extend(temp_m[k])
                tukey_collect2.extend([v+1] * len(temp_m[k]))
        
        print sc, 'threshold:', thresh1+'-'+thresh2
        print pairwise_tukeyhsd(np.array(tukey_collect1), np.array(tukey_collect2))
        print '\n'
    
    plt.figure()
    
    plt.errorbar(range(1, len(expr_int)), collect_m[1], yerr=collect_c[1], 
        label='B', color='#000000', fmt='-o', linewidth=2)
    plt.errorbar(range(1, len(expr_int)), collect_m[2], yerr=collect_c[2], 
        label='C', color='#990033', fmt='-o', linewidth=2)
    plt.errorbar(range(1, len(expr_int)), collect_m[3], yerr=collect_c[3], 
        label='H', color='#FF9900', fmt='-o', linewidth=2)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size':14})
    plt.xticks(np.arange(len(thresh_bins))+1, thresh_bins, fontsize=16)
    plt.yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize=16)
    plt.xlim([0.5, len(thresh_bins)+0.5])
    plt.ylim([0.5, 1.0])
    plt.xlabel('Expression Bins', fontsize=20)
    plt.ylabel('F-Score', fontsize=20)
    plt.savefig('figures/byexpr'+condition+sc.title()+'.png', bbox_inches='tight')





