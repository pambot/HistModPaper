import sys
import numpy as np
import pandas as pd
from scipy.stats import levene, ttest_ind
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.pyplot as plt

test_pairs = [
('predict3/results/rgrExprValid_res.heatmap', 'predict2/results/rgrExprValid_res.heatmap', 'ChromStates vs. HistSig Expr res'),
('predict3/results/rgrSpecValid_res.heatmap', 'predict2/results/rgrSpecValid_res.heatmap', 'ChromStates vs. HistSig Spec res'),

('predict3/results/rgrExprValid_full.heatmap', 'predict2/results/rgrExprValid_full.heatmap', 'ChromStates vs. HistSig Expr full'),
('predict3/results/rgrSpecValid_full.heatmap', 'predict2/results/rgrSpecValid_full.heatmap', 'ChromStates vs. HistSig Spec full'),

('predict3/results/cellvsMeanCodingFscore.heatmap', 'predict2/results/cellvsMeanCodingFscore.heatmap', 'ChromStates vs. HistSig Coding fscore'),
]

def nondiagonal(df):
    m = df.values
    m_nd = np.triu(m, k=1) + np.tril(m, k=-1)
    return [n for n in m_nd.flatten() if n != 0]

pvals = []
for pair in test_pairs:
    df1 = pd.read_csv(pair[0], sep='\t', header=0, index_col=0)
    df2 = pd.read_csv(pair[1], sep='\t', header=0, index_col=0)
    df1_vals, df2_vals = map(nondiagonal, (df1, df2))
    levene_res = levene(df1_vals, df2_vals, center='median')
    pvals.append(levene_res[1])

fdr_adj = multipletests(pvals, method='fdr_bh')

print 'Cell vs. cell for HistSig & ChromStates'
for pair, pv in zip(test_pairs, fdr_adj[1]):
    print '{0}: {1:.10f}'.format(pair[2], pv)


test_sets = [
('predict3/results/rgrExprValid_res.heatmap', 'H Expr Res'), 
('predict2/results/rgrExprValid_res.heatmap', 'C Expr Res'), 
('predict3/results/rgrSpecValid_res.heatmap', 'H Spec Res'), 
('predict2/results/rgrSpecValid_res.heatmap', 'C Spec Res'), 
('predict3/results/rgrExprValid_full.heatmap', 'H Expr Full'), 
('predict2/results/rgrExprValid_full.heatmap', 'C Expr Full'), 
('predict3/results/rgrSpecValid_full.heatmap', 'H Spec Full'), 
('predict2/results/rgrSpecValid_full.heatmap', 'C Spec Full'), 
('predict3/results/cellvsMeanCodingFscore.heatmap', 'H Coding'), 
('predict2/results/cellvsMeanCodingFscore.heatmap', 'C Coding'), 
]

pvals = []
for pair in test_sets:
    df = pd.read_csv(pair[0], sep='\t', header=0, index_col=0)
    diag = df.values.diagonal()
    nondiag = nondiagonal(df)
    ttest_res = ttest_ind(diag, nondiag, equal_var=False)
    pvals.append(ttest_res[1])

fdr_adj = multipletests(pvals, method='fdr_bh')

print '\nCell vs. cell for Diag vs. Non-diag'
for sets, pv in zip(test_sets, fdr_adj[1]):
    print '{0}: {1:.10f}'.format(sets[1], pv)

