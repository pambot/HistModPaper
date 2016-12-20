# change environment
import sys
sys.path.insert(1, '/ifs/home/pw801/bin/python')

# load modules
import pandas as pd
import numpy as np
import scipy
from scipy.stats import linregress
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
import cPickle as pickle

try:
    folder, residuals = sys.argv[1:3]
except (IndexError, ValueError):
    folder, residuals = 'predict1/', '0'

residuals = int(residuals)
cells = ['Gm12878', 'H1hesc', 'Helas3', 'Hepg2', 'Huvec', 'K562', 'Nhek']
m_matrix = pd.DataFrame(index=cells, columns=cells)
c_matrix = pd.DataFrame(index=cells, columns=cells)

def center_data(data):
    return (data-data.mean())/data.std()

def residual(x, y, slope, intercept):
    y_hat = np.dot(x, slope) + intercept
    return y - y_hat

def transform_Xy(X, y, threshold=0):
    pass_vals = y['Multicov'].values > threshold
    new_y = y[pass_vals]
    new_X = X[pass_vals]
    return new_X, new_y

def load(cell, residuals=False):
    X = pd.read_csv(folder + 'train/'+cell+'.matrix', sep='\t', header=0).drop('loci', axis=1)
    features = np.array(X.columns)
    X = X.values
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
    labels = pd.read_csv(folder + 'train/'+cell+'.labels', sep='\t', header=0)
    X, labels = transform_Xy(X, labels)
    expr = center_data(labels['Multicov'].values)
    spec = center_data(labels['Specificity'].values)
    slope, intercept, r_value, p_value, stderr = linregress(expr, spec)
    residues = residual(expr, spec, slope, intercept)
    if residuals:
        return X, residues
    else:
        return X, spec

def confidence_interval(data, conf=0.95):
	a = 1.0*np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * scipy.stats.t._ppf((1+conf)/2., n-1)
	return m, h

for cell1 in cells:
    for cell2 in cells:
        if cell1 == cell2:
            if residuals:
                rgr_scores = pickle.load(open(folder + 'results/rgrSpec_res_'+cell1+'_0.0.pkl', 'rb'))
            else:
                rgr_scores = pickle.load(open(folder + 'results/rgrSpec_full_'+cell1+'_0.0.pkl', 'rb'))
            m_matrix.ix[cell1, cell2] = rgr_scores['r_value']
            print 'Done', cell1, 'vs.', cell2
            continue
        
        X_train, y_train = load(cell1, residuals)
        X_test, y_test = load(cell2, residuals)
        
        rgr = LinearRegression()
        y_pred = rgr.fit(X_train, y_train).predict(X_test)
        y_resub = rgr.fit(X_train, y_train).predict(X_train)
        
        slope, intercept, r_value, p_value, stderr = linregress(y_train, y_resub)
        
        resub = {
        'r_value': r_value,
        }
        
        n_boot = 100
        rng_seed = 0
        boot_scores = {'r_value': []}
        rng = np.random.RandomState(rng_seed)
        while len(boot_scores['r_value']) < n_boot:
	        b_inds = rng.random_integers(0, high=(len(y_pred)-1), size=len(y_pred))
	        if len(np.unique(y_test[b_inds])) < 2:
		        continue
	        
	        slope, intercept, r_value, p_value, stderr = linregress(y_test[b_inds], y_pred[b_inds])
	        boot_scores['r_value'].append(r_value)
        
        # .632 adjustment for bootstrapped AUC
        adj_scores = {'r_value': None}
        ci_scores = {'r_value': None}
        for sc in adj_scores.keys():
            adj_scores[sc] = .632 * np.array(boot_scores[sc]) + .368 * resub[sc]
            mean_sc, pm_sc = confidence_interval(adj_scores[sc], conf=0.95)
            m_matrix.ix[cell1, cell2] = r_value
        
        print 'Done', cell1, 'vs.', cell2

if residuals:
    m_matrix.to_csv(folder + 'results/rgrSpecValid_res.heatmap', sep='\t', header=True, index=True)
else:
    m_matrix.to_csv(folder + 'results/rgrSpecValid_full.heatmap', sep='\t', header=True, index=True)





