import sys
sys.path.insert(0, '/ifs/home/pw801/bin/venv/lib/python2.7/site-packages')

import sklearn
if sklearn.__version__ != '0.18.dev0':
    print 'Sklearn version:', sklearn.__version__
    sys.exit()


# load modules
import pandas as pd
import numpy as np
import scipy
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import cPickle as pickle

try:
    folder, cell, n_layers = sys.argv[1:4]
except (IndexError, ValueError):
    folder, cell, n_layers = 'predict1/', 'Gm12878', '1'

condition = 'Coding'
n_layers = int(n_layers)

# load the data
X = pd.read_csv(folder + 'train/'+cell+'.matrix', sep='\t', header=0).drop('loci', axis=1)
features = np.array(X.columns)
X = X.values

# load the targets
y = pd.read_csv(folder + 'train/'+cell+'.labels', sep='\t', header=0)

# customize the targets
def yCustom(y, condition):
    pass_expr = y['Expr'] == 1
    yi = np.array([i for i in y.index if pass_expr[i]])
    y = y.ix[yi, condition].values
    return y, yi

y, yi = yCustom(y, condition)
X = X[yi]

scaler = StandardScaler()
X = scaler.fit_transform(X) 
clf =  MLPClassifier(
    hidden_layer_sizes=tuple([50] * n_layers),
    alpha=0.001,
    learning_rate_init=0.01,
    activation='logistic',
    random_state=0,
    shuffle=True,
)

cv = StratifiedKFold(y, n_folds=2, shuffle=True, random_state=0)
for train, test in cv:
    pass

def confidence_interval(data, conf=0.95):
	a = 1.0*np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * scipy.stats.t._ppf((1+conf)/2., n-1)
	return m, h

mlp = clf.fit(X[train], y[train])
y_pred = mlp.predict(X[test])
y_resub = mlp.predict(X[train])

resub = {
'auc': roc_auc_score(y[train], y_resub),
'precision': precision_score(y[train], y_resub),
'recall': recall_score(y[train], y_resub),
'fscore': f1_score(y[train], y_resub),
}

n_boot = 100
rng_seed = 0
boot_scores = {'auc': [], 'precision': [], 'recall': [], 'fscore': []}
rng = np.random.RandomState(rng_seed)
while len(boot_scores['auc']) < n_boot:
	b_inds = rng.random_integers(0, high=(len(y_pred)-1), size=len(y_pred))
	if len(np.unique(y[test][b_inds])) < 2:
		continue
	
	boot_scores['auc'].append(roc_auc_score(y[test][b_inds], y_pred[b_inds]))
	boot_scores['precision'].append(precision_score(y[test][b_inds], y_pred[b_inds]))
	boot_scores['recall'].append(recall_score(y[test][b_inds], y_pred[b_inds]))
	boot_scores['fscore'].append(f1_score(y[test][b_inds], y_pred[b_inds]))

# .632 adjustment for bootstrapped AUC
adj_scores = {'auc': None, 'precision': None, 'recall': None, 'fscore': None}
ci_scores = {'auc': None, 'precision': None, 'recall': None, 'fscore': None}
for sc in adj_scores.keys():
    adj_scores[sc] = .632 * np.array(boot_scores[sc]) + .368 * resub[sc]
    mean_sc, pm_sc = confidence_interval(adj_scores[sc], conf=0.95)
    ci_scores[sc] = {'mean': mean_sc, 'pm': pm_sc}

# display
print 'Scores for', cell, condition, n_layers
for sc in ci_scores.keys():
    print sc+':', ci_scores[sc]['mean']

print '\n'

# save coefs
#adj_scores['coefs'] = mlp.coefs_

# what to keep
pickle.dump(adj_scores, open(folder + 'results/mlp'+cell+condition+'_'+str(n_layers)+'.pkl', 'wb'))





