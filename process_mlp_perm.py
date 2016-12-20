import sys
import numpy as np
import scipy.stats
import cPickle as pickle

try:
    folder, cell, condition, ttype = sys.argv[1:5]
except (IndexError, ValueError):
    folder, cell, condition, ttype = 'predict1/', 'Gm12878', 'Coding', ''

if condition == 'Coding':
    ttype = ''

n_perm = 10
n_seed = 20
mlp_perms = []
for seed in range(n_seed):
    coefs = pickle.load(open(folder + 'results/permCoefs{0}{1}{2}_{3}x{4}.pkl'.format(cell, condition, ttype, seed, n_perm)))
    mlp_perms.extend(coefs['perm'])

def compare_perm_coef(d1, d2, mlp_perms, mlp_coefs, c_ind):
    distr = []
    for perm in mlp_perms:
        distr.append(perm[c_ind][d1, d2])
    coef_val = mlp_coefs[c_ind][d1, d2]
    distr.append(coef_val)
    z_scores = scipy.stats.zscore(distr)
    norm_coef_val = z_scores[-1]
    if abs(norm_coef_val) > 1.96:
        return coef_val
    else:
        return 0

mlp_coefs = coefs['real']
f_size, l_size = mlp_perms[0][0].shape
masked_coefs = [np.zeros((f_size, l_size)), np.zeros((l_size, 1))]

for f in range(f_size):
    for l in range(l_size):
        masked_coefs[0][f, l] = compare_perm_coef(f, l, mlp_perms, mlp_coefs, 0)

for l in range(l_size):
    masked_coefs[1][l, 0] = compare_perm_coef(l, 0, mlp_perms, mlp_coefs, 1)

pickle.dump(masked_coefs, open(folder + 'results/mlpMaskCoefs{0}{1}{2}.pkl'.format(cell, condition, ttype), 'wb'))





