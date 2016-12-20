import sys
import math
import pandas as pd
import numpy as np
from scipy.stats import hypergeom, spearmanr
from statsmodels.sandbox.stats.multicomp import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
import cPickle as pickle
import pdb

emissions = pd.read_csv('emissions_27.txt', sep='\t', header=0, index_col=0)

# LIN
# process feature selection by histone
def process_histones(hist_scores):
    histones = emissions.columns
    features = hist_scores.keys()
    fs_scores = np.array(hist_scores.values())
    hist_tss_indexes = []
    hist_tts_indexes = []
    for hist in histones:
        tss_indexes = []
        tts_indexes = []
        for f_ind, feat in enumerate(features):
            if hist in feat and 'TSS' in feat:
                tss_indexes.append(f_ind)
            elif hist in feat and 'TTS' in feat:
                tts_indexes.append(f_ind)
        hist_tss_indexes.append(tss_indexes)
        hist_tts_indexes.append(tts_indexes)
    hist_tss_indexes = np.array(hist_tss_indexes)
    hist_tts_indexes = np.array(hist_tts_indexes)
    hist_tss_scores = {}
    hist_tts_scores = {}
    for hist, hist_tss_ind, hist_tts_ind in zip(histones, hist_tss_indexes, hist_tts_indexes):
        hist_tss_scores[hist] = fs_scores[hist_tss_ind]
        hist_tts_scores[hist] = fs_scores[hist_tts_ind]
    return hist_tss_scores, hist_tts_scores

def decompose_states(hist_scores):
    histones = emissions.columns
    features = hist_scores.keys()
    fs_scores = np.array(hist_scores.values())
    hist_tss_scores = dict.fromkeys(emissions.columns, np.array([]))
    hist_tts_scores = dict.fromkeys(emissions.columns, np.array([]))
    for state in emissions.index:
        state_em = emissions.ix[state]
        state_tss_ind = [i for i,v in enumerate(features) if v.startswith('E{0}_TSS'.format(state))][0]
        state_tts_ind = [i for i,v in enumerate(features) if v.startswith('E{0}_TTS'.format(state))][0]
        for hist in histones:
            hist_tss_scores[hist] = np.append(hist_tss_scores[hist], state_em[hist] * fs_scores[state_tss_ind])
            hist_tts_scores[hist] = np.append(hist_tts_scores[hist], state_em[hist] * fs_scores[state_tts_ind])
    return hist_tss_scores, hist_tts_scores

def process_or_combine(folder, hist_scores):
    if folder in ['predict1/', 'predict3/']:
        hist_tss_scores, hist_tts_scores = process_histones(hist_scores)
    elif folder == 'predict2/':
        hist_tss_scores, hist_tts_scores = decompose_states(hist_scores)
    elif folder == 'predict4/':
        hist_tss_scores = {}
        hist_tts_scores = {}
        hist_tss_scores1, hist_tts_scores1 = process_histones(hist_scores)
        hist_tss_scores2, hist_tts_scores2 = decompose_states(hist_scores)
        for hist in histones:
            hist_tss_scores[hist] = np.append(hist_tss_scores1[hist], hist_tss_scores2[hist])
            hist_tts_scores[hist] = np.append(hist_tts_scores1[hist], hist_tts_scores2[hist])
    return hist_tss_scores, hist_tts_scores

def combine_hist_cells(folder, condition, ttype):
    tss_scores = dict.fromkeys(emissions.columns, np.array([]))
    tts_scores = dict.fromkeys(emissions.columns, np.array([]))
    cells = ['Gm12878', 'H1hesc', 'Helas3', 'Hepg2', 'Huvec', 'K562', 'Nhek']
    for cell in cells:
        hist_scores = pickle.load(open(folder + 
            'results/histScores'+condition+cell+ttype+'.pkl', 'rb'))
        hist_tss, hist_tts = process_or_combine(folder, hist_scores)
        for hist in emissions.columns:
            tss_scores[hist] = np.append(tss_scores[hist], np.sum(hist_tss[hist]))
            tts_scores[hist] = np.append(tts_scores[hist], np.sum(hist_tts[hist]))
    tss_scores, tts_scores = pd.melt(pd.DataFrame(tss_scores)), pd.melt(pd.DataFrame(tts_scores))
    tss_scores.columns, tts_scores.columns = [['Histone Modification', 'Sum of Stability Scores']]*2
    tss_scores['Cell'] = cells * len(emissions.columns)
    tts_scores['Cell'] = cells * len(emissions.columns)
    return tss_scores, tts_scores

# plot average stabilities across all histones
def set_cmap(hist_list):
    histones = emissions.columns
    cmaps = sns.color_palette('Set1', 8) + sns.color_palette('Set2', 8)
    cdict = dict(zip(histones, cmaps))
    cmap_list = [cdict[h] for h in hist_list]
    return cmap_list

def ranked_hist_boxplot(scores, (ymin, ymax), ranked_hist, name):
    scores['HistCat'] = pd.Categorical(scores[scores.columns[0]], categories=ranked_hist, ordered=True)
    scores[scores.columns[1]].fillna(0, inplace=True)
    scores = scores.sort_values('HistCat')
    if 'lin' in name:
        add = 1
    else:
        add = 0
    plt.figure()
    sns.set(style='darkgrid')
    ax = sns.boxplot(x=scores[scores.columns[0]], y=scores[scores.columns[1]], data=scores, palette=set_cmap(ranked_hist))
    plt.xticks(rotation='vertical', fontsize=16)
    plt.yticks(np.linspace(ymin, ymax+add, num=5), fontsize=16)
    plt.ylim([ymin, ymax+add])
    plt.xlabel(scores.columns[0], fontsize=20)
    plt.ylabel(scores.columns[1], fontsize=20)
    plt.savefig('figures/{0}.png'.format(name), bbox_inches='tight')

def plot_lin_ranks(condition, ttype):
    sdict = dict.fromkeys(['P2TSS', 'P3TSS', 'P2TTS', 'P3TTS'])
    sdict['P2TSS'], sdict['P2TTS'] = combine_hist_cells('predict2/', condition, ttype)
    sdict['P3TSS'], sdict['P3TTS'] = combine_hist_cells('predict3/', condition, ttype)
    var_col = sdict['P2TSS'].columns[0]
    val_col = sdict['P2TSS'].columns[1]
    ymin = np.concatenate((sdict['P2TSS'][val_col].values, sdict['P3TSS'][val_col].values)).min()
    ymax = np.concatenate((sdict['P2TSS'][val_col].values, sdict['P3TSS'][val_col].values)).max()
    ymin, ymax = math.floor(ymin), math.ceil(ymax)
    ranked_hist = {}
    for ts in ['TSS', 'TTS']:
        ranked_hist[ts] = sdict['P3'+ts].groupby(var_col).median().sort_values(val_col, ascending=False).index
    for ptype in ['P2', 'P3']:
        for ts in ['TSS', 'TTS']:
            scores = sdict[ptype+ts]
            name = 'linHistRank' + condition + ttype + ptype + ts
            ranked_hist_boxplot(scores, (ymin, ymax), ranked_hist[ts], name)

# calculate whether there's evidence of preferential feature selection in combined models
def calc_feat_overrep_lin(condition, ttype):
    cells = ['Gm12878', 'H1hesc', 'Helas3', 'Hepg2', 'Huvec', 'K562', 'Nhek']
    keys = ['H', 'B', 'C']
    nonzeros = dict.fromkeys(keys, np.zeros(len(cells)))
    f_len = dict.fromkeys(keys, None)
    nz_len = dict.fromkeys(keys, None)
    for c, cell in enumerate(cells):
        hist_scores = pickle.load(open('predict4/' + 
            'results/histScores'+condition+cell+ttype+'.pkl', 'rb'))
        features = hist_scores.keys()
        fs_scores = np.array(hist_scores.values())
        ind = dict.fromkeys(keys, None)
        ind['H'] = np.array([i for i,v in enumerate(features) if v.startswith('H') and v.endswith('3')])
        ind['B'] = np.array([i for i,v in enumerate(features) if v.startswith('H') and v.endswith('1')])
        ind['C'] = np.array([i for i,v in enumerate(features) if v.startswith('E')])
        for k in keys:
            nonzeros[k][c] = len(np.nonzero(fs_scores[ind[k]])[0])
    for k in keys:
        f_len[k] = len(fs_scores[ind[k]])
        nz_len[k] = math.floor(nonzeros[k].mean())
    t_nz = sum([nz_len[k] for k in keys])
    print condition, ttype
    pvals = []
    for k in keys:
        h_score = hypergeom.sf(nz_len[k], len(features), t_nz, f_len[k])
        pvals.append(h_score)
        print k, h_score
    return pvals

# MLP
# modify above for mlp
def garson(coefs):
    nl = len(coefs)
    w, v = coefs
    P = abs(w.T * v)
    Q = P / np.sum(P, axis=1).reshape((-1, 1))
    S = np.nansum(Q, axis=0)
    RI = S / np.sum(S)
    return RI

def combine_mlp_coefs(folder, condition, ttype):
    tss_scores = dict.fromkeys(emissions.columns, np.array([]))
    tts_scores = dict.fromkeys(emissions.columns, np.array([]))
    cells = ['Gm12878', 'H1hesc', 'Helas3', 'Hepg2', 'Huvec', 'K562', 'Nhek']
    for cell in cells:
        mlp_scores = pickle.load(open(folder + 
            'results/mlpMaskCoefs'+cell+condition+ttype+'.pkl', 'rb'))
        fs_scores = garson(mlp_scores)
        with open(folder + 'train/Gm12878.matrix', 'r') as f:
            features = f.readline().rstrip('\n').split('\t')[1:]
        hist_scores = dict(zip(features, fs_scores))
        hist_tss, hist_tts = process_or_combine(folder, hist_scores)
        for hist in emissions.columns:
            tss_scores[hist] = np.append(tss_scores[hist], np.sum(hist_tss[hist]))
            tts_scores[hist] = np.append(tts_scores[hist], np.sum(hist_tts[hist]))
    tss_scores, tts_scores = pd.melt(pd.DataFrame(tss_scores)), pd.melt(pd.DataFrame(tts_scores))
    tss_scores.columns, tts_scores.columns = [['Histone Modification', 'Sum of Stability Scores']]*2
    tss_scores['Cell'] = cells * len(emissions.columns)
    tts_scores['Cell'] = cells * len(emissions.columns)
    return tss_scores, tts_scores

def plot_mlp_ranks(condition, ttype):
    sdict = dict.fromkeys(['P2TSS', 'P3TSS', 'P2TTS', 'P3TTS'])
    sdict['P2TSS'], sdict['P2TTS'] = combine_mlp_coefs('predict2/', condition, ttype)
    sdict['P3TSS'], sdict['P3TTS'] = combine_mlp_coefs('predict3/', condition, ttype)
    var_col = sdict['P2TSS'].columns[0]
    val_col = sdict['P2TSS'].columns[1]
    ymin = np.concatenate([sdict[k][val_col].values for k in ['P2TSS', 'P3TSS', 'P2TTS', 'P3TTS']]).min()
    ymax = np.concatenate([sdict[k][val_col].values for k in ['P2TSS', 'P3TSS', 'P2TTS', 'P3TTS']]).max()
    ranked_hist = {}
    for ts in ['TSS', 'TTS']:
        ranked_hist[ts] = sdict['P3'+ts].groupby(var_col).median().sort_values(val_col, ascending=False).index
    for ptype in ['P2', 'P3']:
        for ts in ['TSS', 'TTS']:
            scores = sdict[ptype+ts]
            name = 'mlpHistRank' + condition + ttype + ptype + ts
            ranked_hist_boxplot(scores, (ymin, ymax), ranked_hist[ts], name)

def plot_lin_mlp_ranks(ptype, condition, ttype):
    keys = [ptype+'TSS', ptype+'TSS']
    if ptype == 'P2':
        folder = 'predict2/'
    elif ptype == 'P3':
        folder = 'predict3/'
    lin_sdict = dict.fromkeys(keys)
    mlp_sdict = dict.fromkeys(keys)
    lin_sdict[ptype+'TSS'], lin_sdict[ptype+'TTS'] = combine_hist_cells(folder, condition, ttype)
    mlp_sdict[ptype+'TSS'], mlp_sdict[ptype+'TTS'] = combine_mlp_coefs(folder, condition, ttype)
    var_col = lin_sdict[ptype+'TSS'].columns[0]
    val_col = lin_sdict[ptype+'TSS'].columns[1]
    for ts in ['TSS', 'TTS']:
        sr = spearmanr(lin_sdict[ptype+ts][val_col], mlp_sdict[ptype+ts][val_col])
        scatter_data = lin_sdict[ptype+ts]
        scatter_data['MLP'] = mlp_sdict[ptype+ts][val_col]
        fname = 'LvsMHistRank' + condition + ttype + ptype + ts
        sns.set(style='darkgrid')
        ax = sns.lmplot(x=val_col, y='MLP', data=scatter_data, hue=var_col, fit_reg=False, legend=False, palette=sns.color_palette('Set1')[:-1]+sns.color_palette('dark'))
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel(val_col+' (Linear)', fontsize=20)
        plt.ylabel(val_col+' (MLP)', fontsize=20)
        if sr[1] < 0.0005:
            sig = '***'
        elif sr[1] < 0.005:
            sig = '**'
        elif sr[1] < 0.05:
            sig = '*'
        else:
            sig = 'NS'
        plt.title('SR = {0:.5f} ({1})'.format(sr[0], sig), fontsize=14, fontweight='bold')
        #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
        plt.savefig('figures/{0}.png'.format(fname), bbox_inches='tight')
    return

#plot_lin_mlp_ranks('P2', 'Multicov', 'full')

# calculate whether there's evidence of preferential feature selection in combined models
def calc_feat_overrep_mlp(condition, ttype):
    cells = ['Gm12878', 'H1hesc', 'Helas3', 'Hepg2', 'Huvec', 'K562', 'Nhek']
    keys = ['H', 'B', 'C']
    nonzeros = dict.fromkeys(keys, np.zeros(len(cells)))
    f_len = dict.fromkeys(keys, None)
    nz_len = dict.fromkeys(keys, None)
    for c, cell in enumerate(cells):
        folder = 'predict4/'
        mlp_scores = pickle.load(open(folder + 
            'results/mlpMaskCoefs'+cell+condition+ttype+'.pkl', 'rb'))
        fs_scores = garson(mlp_scores)
        with open(folder + 'train/Gm12878.matrix', 'r') as f:
            features = f.readline().rstrip('\n').split('\t')[1:]
        hist_scores = dict(zip(features, fs_scores))
        features = hist_scores.keys()
        fs_scores = np.array(hist_scores.values())
        ind = dict.fromkeys(keys, None)
        ind['H'] = np.array([i for i,v in enumerate(features) if v.startswith('H') and v.endswith('3')])
        ind['B'] = np.array([i for i,v in enumerate(features) if v.startswith('H') and v.endswith('1')])
        ind['C'] = np.array([i for i,v in enumerate(features) if v.startswith('E')])
        for k in keys:
            nonzeros[k][c] = len(np.nonzero(fs_scores[ind[k]])[0])
    for k in keys:
        f_len[k] = len(fs_scores[ind[k]])
        nz_len[k] = math.floor(nonzeros[k].mean())
    t_nz = sum([nz_len[k] for k in keys])
    print condition, ttype
    pvals = []
    for k in keys:
        h_score = hypergeom.sf(nz_len[k], len(features), t_nz, f_len[k])
        pvals.append(h_score)
        print k, h_score
    return pvals

# spearman's on the histone rankings for both lin and mlp
def spear_rank_histones1(condition, ttype, mlp=False):
    keys = ['P2TSS', 'P3TSS', 'P2TTS', 'P3TTS']
    sdict = dict.fromkeys(keys)
    if not mlp:
        sdict['P2TSS'], sdict['P2TTS'] = combine_hist_cells('predict2/', condition, ttype)
        sdict['P3TSS'], sdict['P3TTS'] = combine_hist_cells('predict3/', condition, ttype)
    else:
        sdict['P2TSS'], sdict['P2TTS'] = combine_mlp_coefs('predict2/', condition, ttype)
        sdict['P3TSS'], sdict['P3TTS'] = combine_mlp_coefs('predict3/', condition, ttype)
    var_col = sdict['P2TSS'].columns[0]
    val_col = sdict['P2TSS'].columns[1]
    assert all(x1==x2 and x2==x3 and x3==x4 for x1, x2, x3, x4 in zip(sdict['P2TSS'][var_col], sdict['P3TSS'][var_col], sdict['P2TTS'][var_col], sdict['P3TTS'][var_col]))
    assert all(x1==x2 and x2==x3 and x3==x4 for x1, x2, x3, x4 in zip(sdict['P2TSS']['Cell'], sdict['P3TSS']['Cell'], sdict['P2TTS']['Cell'], sdict['P3TTS']['Cell']))
    if not mlp:
        mlp_status = 'Linear'
    else:
        mlp_status = 'MLP'
    print 'Spearman R for {0} {1} {2}'.format(condition, ttype, mlp_status)
    pvals = []
    for ts in ['TSS', 'TTS']:
        sr = spearmanr(sdict['P3'+ts][val_col], sdict['P2'+ts][val_col])
        pvals.append(sr[1])
        print '['+ts+'] R: {0:.5f}\tP-value: {1:.10f}'.format(sr[0], sr[1])
    return pvals

def spear_rank_histones2(condition, ttype):
    keys = ['P2TSS', 'P3TSS', 'P2TTS', 'P3TTS']
    lin_sdict = dict.fromkeys(keys)
    mlp_sdict = dict.fromkeys(keys)
    lin_sdict['P2TSS'], lin_sdict['P2TTS'] = combine_hist_cells('predict2/', condition, ttype)
    mlp_sdict['P2TSS'], mlp_sdict['P2TTS'] = combine_mlp_coefs('predict2/', condition, ttype)
    lin_sdict['P3TSS'], lin_sdict['P3TTS'] = combine_hist_cells('predict3/', condition, ttype)
    mlp_sdict['P3TSS'], mlp_sdict['P3TTS'] = combine_mlp_coefs('predict3/', condition, ttype)
    var_col = lin_sdict['P2TSS'].columns[0]
    val_col = lin_sdict['P2TSS'].columns[1]
    assert all(x1==x2 and x2==x3 and x3==x4 for x1, x2, x3, x4 in zip(lin_sdict['P2TSS'][var_col], lin_sdict['P3TSS'][var_col], mlp_sdict['P2TSS'][var_col], mlp_sdict['P3TSS'][var_col]))
    assert all(x1==x2 and x2==x3 and x3==x4 for x1, x2, x3, x4 in zip(lin_sdict['P2TTS'][var_col], lin_sdict['P3TTS'][var_col], mlp_sdict['P2TTS'][var_col], mlp_sdict['P3TTS'][var_col]))
    print 'Spearman R for {0} {1}'.format(condition, ttype)
    pvals = []
    for k in keys:
        sr = spearmanr(lin_sdict[k][val_col], mlp_sdict[k][val_col])
        pvals.append(sr[1])
        print '['+k+'] R: {0:.5f}\tP-value: {1:.10f}'.format(sr[0], sr[1])
    return pvals

"""
plot_lin_ranks('Coding', '')
plot_lin_ranks('Multicov', 'full')
plot_lin_ranks('Multicov', 'res')
plot_lin_ranks('Specificity', 'full')
plot_lin_ranks('Specificity', 'res')

plot_mlp_ranks('Coding', '')
plot_mlp_ranks('Multicov', 'full')
#plot_mlp_ranks('Multicov', 'res')
plot_mlp_ranks('Specificity', 'full')
#plot_mlp_ranks('Specificity', 'res')
"""

for ptype in ['P2', 'P3']:
    plot_lin_mlp_ranks(ptype, 'Coding', '')
    plot_lin_mlp_ranks(ptype, 'Multicov', 'full')
    #plot_lin_mlp_ranks(ptype, 'Multicov', 'res')
    plot_lin_mlp_ranks(ptype, 'Specificity', 'full')
    #plot_lin_mlp_ranks(ptype, 'Specificity', 'res')

pval_collect1 = []
print 'Hypergeom SF scores for Linear:'
pval_collect1.extend(calc_feat_overrep_lin('Coding', ''))
pval_collect1.extend(calc_feat_overrep_lin('Multicov', 'full'))
pval_collect1.extend(calc_feat_overrep_lin('Multicov', 'res'))
pval_collect1.extend(calc_feat_overrep_lin('Specificity', 'full'))
pval_collect1.extend(calc_feat_overrep_lin('Specificity', 'res'))

pval_collect2 = []
print 'Hypergeom SF scores for MLP:'
pval_collect2.extend(calc_feat_overrep_mlp('Coding', ''))
pval_collect2.extend(calc_feat_overrep_mlp('Multicov', 'full'))
pval_collect2.extend(calc_feat_overrep_mlp('Multicov', 'res'))
pval_collect2.extend(calc_feat_overrep_mlp('Specificity', 'full'))
pval_collect2.extend(calc_feat_overrep_mlp('Specificity', 'res'))

print 'Adjusted p-values for hypergeom Lin'
fdr_adj1 = multipletests(pval_collect1, method='fdr_bh')
print ' '.join(['{0:.10f}'.format(q) for q in fdr_adj1[1]])

print 'Asjusted p-values for hypergeom MLP'
fdr_adj2 = multipletests(pval_collect2, method='fdr_bh')
print ' '.join(['{0:.10f}'.format(q) for q in fdr_adj2[1]])


pval_collect1 = []
for b in [False, True]:
    pval_collect1.extend(spear_rank_histones1('Coding', '', mlp=b))
    pval_collect1.extend(spear_rank_histones1('Multicov', 'full', mlp=b))
    #pval_collect1.extend(spear_rank_histones1('Multicov', 'res', mlp=b))
    pval_collect1.extend(spear_rank_histones1('Specificity', 'full', mlp=b))
    #pval_collect1.extend(spear_rank_histones1('Specificity', 'res', mlp=b))

pval_collect2 = []
pval_collect2.extend(spear_rank_histones2('Coding', ''))
pval_collect2.extend(spear_rank_histones2('Multicov', 'full'))
#pval_collect2.extend(spear_rank_histones2('Multicov', 'res'))
pval_collect2.extend(spear_rank_histones2('Specificity', 'full'))
#pval_collect2.extend(spear_rank_histones2('Specificity', 'res'))

print 'Adjusted p-values for P2 vs. P3'
fdr_adj1 = multipletests(pval_collect1, method='fdr_bh')
print ' '.join(['{0:.10f}'.format(q) for q in fdr_adj1[1]])

print 'Asjusted p-values for Lin vs. MLP'
fdr_adj2 = multipletests(pval_collect2, method='fdr_bh')
print ' '.join(['{0:.10f}'.format(q) for q in fdr_adj2[1]])



