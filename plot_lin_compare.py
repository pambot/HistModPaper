# load modules
import sys
import glob
import numpy as np
import pandas as pd
import scipy.stats
from statsmodels.sandbox.stats.multicomp import multipletests
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import seaborn as sns
import matplotlib.pyplot as plt
import cPickle as pickle

cells = ['Gm12878', 'H1hesc', 'Helas3', 'Hepg2', 'Huvec', 'K562', 'Nhek']

def confidence_interval(data, conf=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t._ppf((1+conf)/2., n-1)
    return m, h

def get_vals(folder, fname_format, ttype, layer=None):
    vals = []
    for cell in cells:
        if layer:
            fname = fname_format.format(cell, layer)
        else:
            fname = fname_format.format(cell)
        data = pickle.load(open(folder + 'results/' + fname, 'rb'))
        vals.append(np.mean(data[ttype]))
    return np.array(vals)

# plot all feature types in order with horizontal bars of the mean and CI of the combined linear
def plot_linear(exp, fname_format, ttype, yticks):
    keys = ['B', 'C', 'H']
    
    v_data = dict.fromkeys(keys, None)
    for n, k in zip(range(1, 4), keys):
        folder = 'predict{0}/'.format(n)
        v_data[k] = get_vals(folder, fname_format, ttype)
    
    ld = len(v_data.keys())
    
    c_data = get_vals('predict4/', fname_format, ttype)
    c_mean, c_ci = confidence_interval(c_data)
    
    fdr_collect = []
    for k in keys:
        fdr_collect.append(scipy.stats.ttest_ind(v_data[k], c_data)[1])
    
    fdr_adj = multipletests(fdr_collect, method='fdr_bh')
    
    if any(x > 0.05 for x in fdr_adj[1]):
        pass_or_fail = 'FAIL'
    else:
        pass_or_fail = 'PASS'
    
    print exp, ttype, ':', ' '.join(['{0:.10f}'.format(a) for a in fdr_adj[1]]), pass_or_fail
    
    tukey_collect1 = []
    tukey_collect2 = []
    for v, k in enumerate(keys):
        tukey_collect1.extend(v_data[k])
        tukey_collect2.extend([v+1] * len(v_data[k]))
    
    print pairwise_tukeyhsd(np.array(tukey_collect1), np.array(tukey_collect2))
    print '\n'
    
    rename = {
    'fscore': 'F-Score',
    'auc': 'ROC AUC',
    'r_value': 'Pearson R',
    }
    
    v_df = pd.DataFrame(v_data)
    v_df = pd.melt(v_df)
    v_df.columns = ['Feature Type', rename[ttype]]
    
    plt.figure(figsize=(3,6))
    plt.axhline(y=c_mean, color='k', linewidth=2, linestyle='--')
    plt.fill_between([-1, ld+1.0], c_mean-c_ci, c_mean+c_ci, facecolor='k', alpha=0.2)
    
    sns.set(style='darkgrid')
    pal = sns.cubehelix_palette(len(v_data.keys()), rot=0.4, dark=0.3)
    ax = sns.boxplot(x=v_df[v_df.columns[0]], y=v_df[v_df.columns[1]], data=v_df, palette=pal)
    plt.xticks(fontsize=16)
    plt.yticks(yticks, fontsize=16)
    plt.xlabel(v_df.columns[0], fontsize=20)
    plt.ylabel(v_df.columns[1], fontsize=20)
    plt.savefig('figures/lin_{0}_{1}.png'.format(exp, ttype), bbox_inches='tight')
    return

# run code
plot_linear('Coding', 'scores{0}Coding_0.0-20.0.pkl', 'fscore', (np.arange(5, 11)/10.0))
plot_linear('Coding', 'scores{0}Coding_0.0-20.0.pkl', 'auc', (np.arange(5, 11)/10.0))

plot_linear('Expr_full', 'rgrExpr_full_{0}_0.0.pkl', 'r_value', (np.arange(5, 11)/10.0))
plot_linear('Spec_full', 'rgrSpec_full_{0}_0.0.pkl', 'r_value', (np.arange(5, 11)/10.0))

plot_linear('Expr_res', 'rgrExpr_res_{0}_0.0.pkl', 'r_value', (np.arange(2, 11)/10.0))
plot_linear('Spec_res', 'rgrSpec_res_{0}_0.0.pkl', 'r_value', (np.arange(2, 11)/10.0))




