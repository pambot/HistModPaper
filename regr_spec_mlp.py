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
from scipy.stats import linregress
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score
from sklearn.model_selection import LabelKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import cPickle as pickle

try:
    folder, cell, ttype, threshold, n_layers = sys.argv[1:6]
except (IndexError, ValueError):
    folder, cell, ttype, threshold, n_layers = 'predict1/', 'Gm12878', 'res', '0.0', '1'

condition = 'Specificity'
threshold = float(threshold)
n_layers = int(n_layers)

# load the data
X = pd.read_csv(folder + 'train/'+cell+'.matrix', sep='\t', header=0).drop('loci', axis=1)
features = np.array(X.columns)
X = X.values

scaler = StandardScaler()
X = scaler.fit_transform(X)

# load the targets
def center_data(data):
    return (data-data.mean())/data.std()

def residual(x, y, slope, intercept):
    y_hat = np.dot(x, slope) + intercept
    return y - y_hat

# y is the residual of the correlation between specificity and expression
labels = pd.read_csv(folder + 'train/'+cell+'.labels', sep='\t', header=0)

def transform_Xy(X, y, threshold=threshold):
    pass_vals = y['Multicov'].values > threshold
    new_y = y[pass_vals]
    new_X = X[pass_vals]
    return new_X, new_y

X, labels = transform_Xy(X, labels)

expr = center_data(labels['Multicov'].values)
spec = center_data(labels['Specificity'].values)
slope, intercept, r_value, p_value, stderr = linregress(expr, spec)
residues = residual(expr, spec, slope, intercept)

if ttype=='full':
    y = spec
elif ttype=='res':
    y = residues

clf =  MLPRegressor(
    hidden_layer_sizes=tuple([50] * n_layers),
    alpha=0.001,
    learning_rate_init=0.01,
    activation='logistic',
    random_state=0,
    shuffle=True,
)

cv = LabelKFold(y, n_folds=2)
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

slope, intercept, r_value, p_value, stderr = linregress(y[train], y_resub)
resub = {
'r_value': r_value,
}

n_boot = 100
rng_seed = 0
boot_scores = {'r_value': []}
rng = np.random.RandomState(rng_seed)
while len(boot_scores['r_value']) < n_boot:
	b_inds = rng.random_integers(0, high=(len(y_pred)-1), size=len(y_pred))
	if len(np.unique(y[test][b_inds])) < 2:
		continue
	
	slope, intercept, r_value, p_value, stderr = linregress(y[test][b_inds], y_pred[b_inds])
	boot_scores['r_value'].append(r_value)

# .632 adjustment for bootstrapped AUC
adj_scores = {'r_value': None}
ci_scores = {'r_value': None}
for sc in adj_scores.keys():
    adj_scores[sc] = .632 * np.array(boot_scores[sc]) + .368 * resub[sc]
    mean_sc, pm_sc = confidence_interval(adj_scores[sc], conf=0.95)
    ci_scores[sc] = {'mean': mean_sc, 'pm': pm_sc}

# display
print 'Scores for', cell, condition, ttype, n_layers
for sc in ci_scores.keys():
    print sc+':', ci_scores[sc]['mean']

print '\n'

# save coefs
#adj_scores['coefs'] = mlp.coefs_

# what to keep
pickle.dump(adj_scores, open(folder + 'results/mlp'+cell+condition+ttype.title()+'_'+str(n_layers)+'.pkl', 'wb'))





