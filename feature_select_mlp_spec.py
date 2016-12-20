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
from sklearn.neural_network import MLPRegressor, MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import cPickle as pickle

try:
    folder, cell, ttype, n_layers, seed = sys.argv[1:6]
except (IndexError, ValueError):
    folder, cell, ttype, n_layers, seed = 'predict1/', 'Gm12878', 'full', '1', '0'

condition = 'Specificity'
threshold = 0.0
n_layers = int(n_layers)
seed = int(seed)

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

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=0)

coefs = {'real':None, 'perm':None}

# mlp real fit
mlp = clf.fit(X_train, y_train)
coefs['real'] = mlp.coefs_

# mlp permutations
np.random.seed(seed=seed)
n_perm = 10
perm_coefs = []
for n in range(n_perm):
    y_perm = np.random.permutation(y_train)
    mlp_perm = clf.fit(X_train, y_perm)
    perm_coefs.append(mlp_perm.coefs_)

coefs['perm'] = perm_coefs

pickle.dump(coefs, open(folder + 'results/permCoefs{0}{1}{2}_{3}x{4}.pkl'.format(cell, condition, ttype, seed, n_perm), 'wb'))
    

