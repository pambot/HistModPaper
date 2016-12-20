# change environment
import sys
sys.path.insert(1, '/ifs/home/pw801/bin/python')

# load modules
import pandas as pd
import numpy as np
import scipy
from scipy.stats import linregress
from sklearn.cross_validation import train_test_split
from sklearn.linear_model import RandomizedLasso
from sklearn.preprocessing import StandardScaler
import cPickle as pickle

try:
    folder, cell, exp = sys.argv[1:4]
except (IndexError, ValueError):
    folder, cell, exp = 'predict4/', 'Gm12878', 'full'

threshold = 0.0
condition = 'Specificity'

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
    pass_vals = y['Multicov'].values >= threshold
    new_y = y[pass_vals]
    new_X = X[pass_vals]
    return new_X, new_y

X, labels = transform_Xy(X, labels)

expr = center_data(labels['Multicov'].values)
spec = center_data(labels['Specificity'].values)
slope, intercept, r_value, p_value, stderr = linregress(expr, spec)
residues = residual(expr, spec, slope, intercept)

if exp == 'full':
    y = spec
elif exp == 'res':
    y = residues

X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.5, random_state=0)

rl = RandomizedLasso()
fs = rl.fit(X_train, y_train)
hist_scores = dict(zip(features, fs.scores_))

pickle.dump(hist_scores, open(folder + 'results/histScores'+condition+cell+exp+'.pkl', 'wb'))


