# change environment
import sys
sys.path.insert(1, '/ifs/home/pw801/bin/python')

# load modules
import pandas as pd
import numpy as np
import scipy
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import cPickle as pickle

try:
    folder, condition, test, show_sc = sys.argv[1:5]
except (IndexError, ValueError):
    folder, condition, test, show_sc = 'predict1/', 'Coding', 'S2', 'fscore'

cells = ['Gm12878', 'H1hesc', 'Helas3', 'Hepg2', 'Huvec', 'K562', 'Nhek']
m_matrix = pd.DataFrame(index=cells, columns=cells)
c_matrix = pd.DataFrame(index=cells, columns=cells)

def yCustom(y, condition, test):
    pass_expr = y['Expr'] == 1
    if test == 'S1':
        yi = y.index
    elif test == 'S2':
        yi = np.array([i for i in y.index if pass_expr[i]])
    y = y.ix[yi, condition].values
    return y, yi

def load(cell):
	X = pd.read_csv(folder + 'train/'+cell+'.matrix', sep='\t', header=0).drop('loci', axis=1)
	features = np.array(X.columns)
	X = X.values
	y = pd.read_csv(folder + 'train/'+cell+'.labels', sep='\t', header=0)
	y, yi = yCustom(y, condition, test)
	X = X[yi]
	return X, y

scaler = 'scaler'
mstep = 'logit'
pipe = [
	(scaler, StandardScaler()), 
	(mstep, LogisticRegression(
		C=1, 
		penalty='l2',
		class_weight='balanced', 
		random_state=0, 
		solver='liblinear'
		))
	]

clf = Pipeline(pipe)

def confidence_interval(data, conf=0.95):
	a = 1.0*np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * scipy.stats.t._ppf((1+conf)/2., n-1)
	return m, h

for cell1 in cells:
    for cell2 in cells:
        if cell1 == cell2:
            adj_scores = pickle.load(open(folder + 'results/scores'+cell1+condition+'_0.0-20.0.pkl', 'rb'))
            mean_sc, pm_sc = confidence_interval(adj_scores[show_sc], conf=0.95)
            m_matrix.ix[cell1, cell2] = mean_sc
            c_matrix.ix[cell1, cell2] = pm_sc
            continue
        
        X_train, y_train = load(cell1)
        X_test, y_test = load(cell2)
        
        y_resub = clf.fit(X_train, y_train).predict(X_train)
        y_pred = clf.fit(X_train, y_train).predict(X_test)
        
        resub = {
        'auc': roc_auc_score(y_train, y_resub),
        'precision': precision_score(y_train, y_resub),
        'recall': recall_score(y_train, y_resub),
        'fscore': f1_score(y_train, y_resub),
        }
        
        n_boot = 100
        rng_seed = 0
        boot_scores = {'auc': [], 'precision': [], 'recall': [], 'fscore': []}
        rng = np.random.RandomState(rng_seed)
        while len(boot_scores['auc']) < n_boot:
	        b_inds = rng.random_integers(0, high=(len(y_pred)-1), size=len(y_pred))
	        if len(np.unique(y_test[b_inds])) < 2:
		        continue
	        
	        boot_scores['auc'].append(roc_auc_score(y_test[b_inds], y_pred[b_inds]))
	        boot_scores['precision'].append(precision_score(y_test[b_inds], y_pred[b_inds]))
	        boot_scores['recall'].append(recall_score(y_test[b_inds], y_pred[b_inds]))
	        boot_scores['fscore'].append(f1_score(y_test[b_inds], y_pred[b_inds]))
        
        # .632 adjustment for bootstrapped AUC
        adj_scores = {'auc': None, 'precision': None, 'recall': None, 'fscore': None}
        ci_scores = {'auc': None, 'precision': None, 'recall': None, 'fscore': None}
        for sc in adj_scores.keys():
            adj_scores[sc] = .632 * np.array(boot_scores[sc]) + .368 * resub[sc]
            mean_sc, pm_sc = confidence_interval(adj_scores[sc], conf=0.95)
            ci_scores[sc] = {'mean': mean_sc, 'pm': pm_sc}
        
        m_matrix.ix[cell1, cell2] = ci_scores[show_sc]['mean']
        print 'Done', condition+':', cell1, 'vs.', cell2


m_matrix.to_csv(folder + 'results/cellvsMean' + condition + show_sc[0].upper() + show_sc[1:] + '.heatmap', sep='\t', header=True, index=True)





