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
    folder, cell, n_layers, seed = sys.argv[1:5]
except (IndexError, ValueError):
    folder, cell, n_layers, seed = 'predict1/', 'Gm12878', '1', '0'

condition = 'Coding'
n_layers = int(n_layers)
seed = int(seed)

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

pickle.dump(coefs, open(folder + 'results/permCoefs{0}{1}_{2}x{3}.pkl'.format(cell, condition, seed, n_perm), 'wb'))
    

