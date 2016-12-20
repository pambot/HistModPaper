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

# plot all number of layers in order with horizontal bars of the mean and CI of the combined MLP
def plot_layers(exp, folder, fname_mlp, ttype, yticks):
    fdict = {'predict1/':'B', 'predict2/':'C', 'predict3/':'H'}
    f_str = fdict[folder]
    
    v_data = dict.fromkeys(range(1, 6), None)
    for n in range(1, 6):
        v_data[n] = get_vals(folder, fname_mlp, ttype, layer=n)
    
    ld = len(v_data.keys())
    
    c_data = get_vals('predict4/', fname_mlp, ttype, layer=n)
    c_mean, c_ci = confidence_interval(c_data)
    
    collect = []
    for n in range(1, 6):
        collect.append(scipy.stats.ttest_ind(v_data[n], c_data)[1])
    
    adj_coll = multipletests(collect, method='fdr_bh')
    print exp, f_str, ttype, ':', ' '.join(['{0:.10f}'.format(a) for a in adj_coll[1]])
    
    rename = {
    'fscore': 'F-Score',
    'auc': 'ROC AUC',
    'r_value': 'Pearson R',
    }
    
    v_df = pd.DataFrame(v_data)
    v_df = pd.melt(v_df)
    v_df.columns = ['Layer', rename[ttype]]
    
    plt.figure(figsize=(4,6))
    plt.axhline(y=c_mean, color='k', linewidth=2, linestyle='--')
    plt.fill_between([-1, ld+1.0], c_mean-c_ci, c_mean+c_ci, facecolor='k', alpha=0.2)
    
    sns.set(style='darkgrid')
    pal = sns.cubehelix_palette(len(v_data.keys()), rot=-0.4, dark=0.3)
    ax = sns.boxplot(x=v_df[v_df.columns[0]], y=v_df[v_df.columns[1]], data=v_df, palette=pal)
    plt.xticks(fontsize=16)
    plt.yticks(yticks, fontsize=16)
    plt.xlabel(v_df.columns[0], fontsize=20)
    plt.ylabel(v_df.columns[1], fontsize=20)
    plt.savefig('figures/mlp_{0}_{1}_{2}.png'.format(exp, f_str, ttype), bbox_inches='tight')
    #plt.show()
    return

#plot_layers('Coding', 'predict{0}/'.format(1), 'mlp{0}Coding_{1}.pkl', 'fscore', (np.arange(90, 102)/100.0)[::2])


def plot_lin_vs_mlp(condition, fname_lin, fname_mlp, ttype, yticks):
    keys1 = ['B', 'C', 'H']
    keys2 = ['Lin', 'MLP']
    
    v_data = {}
    for n, k in zip(range(1, 4), keys1):
        folder = 'predict{0}/'.format(n)
        for k2 in keys2:
            if k2=='Lin':
                v_data[k+'-'+k2] = get_vals(folder, fname_lin, ttype)
            elif k2=='MLP':
                v_data[k+'-'+k2] = get_vals(folder, fname_mlp, ttype, layer=1)
    
    collect = []
    for k in keys1:
        collect.append(scipy.stats.ttest_ind(v_data[k+'-Lin'], v_data[k+'-MLP'])[1])
    
    adj_coll = multipletests(collect, method='fdr_bh')
    print 'Tests for B C H, Lin vs. MLP:', ' '.join(['{0:.10f}'.format(a) for a in adj_coll[1]])
    
    rename = {
    'fscore': 'F-Score',
    'auc': 'ROC AUC',
    'r_value': 'Pearson R',
    }
    
    v_df = pd.DataFrame(v_data)
    v_df = pd.melt(v_df)
    v_df.columns = ['Feature Set', rename[ttype]]
    
    l_or_m = []
    for n in v_df.index:
        if 'Lin' in v_df.ix[n, 'Feature Set']:
            l_or_m.append('Linear')
        else:
            l_or_m.append('MLP')
    
    v_df['Linear/MLP'] = l_or_m
    v_df['Feature Set'] = [l.split('-')[0] for l in v_df['Feature Set']]
    
    plt.figure(figsize=(4,6))
    sns.set(style='darkgrid')
    ax = sns.boxplot(x=v_df[v_df.columns[0]], y=v_df[v_df.columns[1]], hue=v_df[v_df.columns[2]], data=v_df, palette='Paired')
    plt.xticks(fontsize=16)
    plt.yticks(yticks, fontsize=16)
    plt.xlabel(v_df.columns[0], fontsize=20)
    plt.ylabel(v_df.columns[1], fontsize=20)
    plt.savefig('figures/lin_mlp_compare'+condition+'.png', bbox_inches='tight')
    #plt.show()
    return

plot_lin_vs_mlp('Coding', 'scores{0}Coding_0.0-20.0.pkl', 'mlp{0}Coding_{1}.pkl', 'fscore', (np.arange(5, 11)/10.0))
plot_lin_vs_mlp('Expr', 'rgrExpr_full_{0}_0.0.pkl', 'mlp{0}MulticovFull_{1}.pkl', 'r_value', (np.arange(5, 11)/10.0))
plot_lin_vs_mlp('Spec', 'rgrSpec_full_{0}_0.0.pkl', 'mlp{0}SpecificityFull_{1}.pkl', 'r_value', (np.arange(5, 11)/10.0))


for n in range(1, 4):
    plot_layers('Coding', 'predict{0}/'.format(n), 'mlp{0}Coding_{1}.pkl', 'fscore', (np.arange(90, 102)/100.0)[::2])

for n in range(1, 4):
    plot_layers('Coding', 'predict{0}/'.format(n), 'mlp{0}Coding_{1}.pkl', 'auc', (np.arange(50, 84)/100.0)[::4])

for n in range(1, 4):
    plot_layers('Expr_full', 'predict{0}/'.format(n), 'mlp{0}MulticovFull_{1}.pkl', 'r_value', (np.arange(5, 11)/10.0))

for n in range(1, 4):
    plot_layers('Spec_full', 'predict{0}/'.format(n), 'mlp{0}SpecificityFull_{1}.pkl', 'r_value', (np.arange(5, 11)/10.0))

for n in range(1, 4):
    plot_layers('Expr_res', 'predict{0}/'.format(n), 'mlp{0}MulticovRes_{1}.pkl', 'r_value', (np.arange(5, 11)/10.0))

for n in range(1, 4):
    plot_layers('Spec_res', 'predict{0}/'.format(n), 'mlp{0}SpecificityRes_{1}.pkl', 'r_value', (np.arange(5, 11)/10.0))




