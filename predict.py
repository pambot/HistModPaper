# change environment
import sys
sys.path.insert(1, '/ifs/home/pw801/bin/python')

# load modules
import pandas as pd
import numpy as np
import scipy
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import cPickle as pickle

try:
    folder, cell, condition, test, thresh1, thresh2 = sys.argv[1:7]
except (IndexError, ValueError):
    folder, cell, condition, test, thresh1, thresh2 = 'predict1/', 'Gm12878', 'Coding', 'S2', '0.0', '20.0'

# load the data
X = pd.read_csv(folder + 'train/'+cell+'.matrix', sep='\t', header=0).drop('loci', axis=1)
features = np.array(X.columns)
X = X.values

# load the targets
y = pd.read_csv(folder + 'train/'+cell+'.labels', sep='\t', header=0)

# customize the targets
def yCustom(y, condition, test, thresh1, thresh2):
    thresh1, thresh2 = float(thresh1), float(thresh2)
    pass_expr = y['Expr'] == 1
    pass_thresh1 = y['Multicov'].values > thresh1
    pass_thresh2 = y['Multicov'].values <= thresh2
    if test == 'S1':
        yi = y.index
    elif test == 'S2':
        yi = np.array([i for i in y.index 
            if pass_expr[i] and pass_thresh1[i] and pass_thresh2[i]])
    y = y.ix[yi, condition].values
    return y, yi

y, yi = yCustom(y, condition, test, thresh1, thresh2)

# tests
if not yi.any():
    print 'Selected indexes empty.'
    sys.exit()

if len(y[y==1]) < 2:
	print 'Scores for', cell, condition, 'between', thresh1, '-', thresh2
	print 'Not enough positive cases found.\n'
	sys.exit()

# get relevant X
X = X[yi]

# set classifier
scaler = StandardScaler()
X = scaler.fit_transform(X)

clf = LogisticRegression(
	C=1, 
	penalty='l2',
	class_weight='balanced', 
	random_state=0, 
	solver='liblinear'
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

lgr = clf.fit(X[train], y[train])
y_pred = lgr.predict(X[test])
y_resub = lgr.predict(X[train])

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
print 'Scores for', cell, condition, 'between', thresh1, '-', thresh2
for sc in ci_scores.keys():
    print sc+':', ci_scores[sc]['mean']

print '\n'

# save coefs
adj_scores['coefs'] = lgr.coef_

# what to keep
pickle.dump(adj_scores, open(folder + 'results/scores'+cell+condition+'_'+thresh1+'-'+thresh2+'.pkl', 'wb'))





