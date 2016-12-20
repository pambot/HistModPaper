# load modules
import sys
import glob
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cPickle as pickle


plt.style.use('ggplot')

cells = ['Gm12878', 'H1hesc', 'Helas3', 'Hepg2', 'Huvec', 'K562', 'Nhek']

folder = {
1:'../predict1/results/',
2:'../predict2/results/',
3:'../predict3/results/',
}

plt.figure()
cdict = {'full': '#003366', 'res':'#900000'}
for ttype in ['full', 'res']:
    v_means = []
    for fol in range(1,4):
        vals = []
        for cell in cells:
            data = pickle.load(open(folder[fol]+'rgrExpr_'+ttype+'_'+cell+'_1.0.pkl', 'rb'))
            vals.append(data['r_value'])
        
        v_means.append(np.array(vals))
    
    if ttype=='full':
        pos = np.array([0.5, 2.0, 3.5])
    else:
        pos = np.array([1.0, 2.5, 4.0])
    violin = plt.violinplot(v_means, pos, widths=0.3, showmeans=True, showextrema=False)
    plt.setp(violin['bodies'], facecolor=cdict[ttype], edgecolor=cdict[ttype])
    for key in ['cmeans']:
	    plt.setp(violin[key], color=cdict[ttype], linewidth='3', alpha=0.5)

blue_patch = mpatches.Patch(color='#003366', label='Full data') 
red_patch = mpatches.Patch(color='#900000', label='Residuals')
plt.legend(handles=[blue_patch, red_patch])
plt.xticks([0.75, 2.25, 3.75], ['Binary', 'ChromStates', 'Signals'], fontsize=16, color='k')
plt.yticks(np.arange(3, 11)/10.0, fontsize=16, color='k')
plt.ylabel('Pearson R', fontsize=18, color='k')
plt.title('Expression', fontsize=18)
plt.savefig('figures/r2Expr.png', bbox_inches='tight')


plt.figure()
cdict = {'full': '#003366', 'res':'#900000'}
for ttype in ['full', 'res']:
    v_means = []
    for fol in range(1,4):
        vals = []
        for cell in cells:
            data = pickle.load(open(folder[fol]+'rgrSpec_'+ttype+'_'+cell+'_1.0.pkl', 'rb'))
            vals.append(data['r_value'])
        
        v_means.append(np.array(vals))
    
    if ttype=='full':
        pos = np.array([0.5, 2.0, 3.5])
    else:
        pos = np.array([1.0, 2.5, 4.0])
    violin = plt.violinplot(v_means, pos, widths=0.3, showmeans=True, showextrema=False)
    plt.setp(violin['bodies'], facecolor=cdict[ttype], edgecolor=cdict[ttype])
    for key in ['cmeans']:
	    plt.setp(violin[key], color=cdict[ttype], linewidth='3', alpha=0.5)

blue_patch = mpatches.Patch(color='#003366', label='Full data') 
red_patch = mpatches.Patch(color='#900000', label='Residuals')
plt.legend(handles=[blue_patch, red_patch])
plt.xticks([0.75, 2.25, 3.75], ['Binary', 'ChromStates', 'Signals'], fontsize=16, color='k')
plt.yticks(np.arange(3, 11)/10.0, fontsize=16, color='k')
plt.ylabel('Pearson R', fontsize=18, color='k')
plt.title('Specificity', fontsize=18)
plt.savefig('figures/r2Spec.png', bbox_inches='tight')


from scipy.stats import ttest_ind
ttest_ind(v_means[1], v_means[2])


