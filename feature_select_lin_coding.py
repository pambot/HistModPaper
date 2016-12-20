# change environment
import sys
sys.path.insert(1, '/ifs/home/pw801/bin/python')

# load modules
import pandas as pd
import numpy as np
import scipy
from sklearn.linear_model import RandomizedLogisticRegression
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import StandardScaler
import cPickle as pickle

try:
    folder, cell = sys.argv[1:3]
except (IndexError, ValueError):
    folder, cell = 'predict2/', 'Gm12878'

condition = 'Coding'

# load the data
X = pd.read_csv(folder + 'train/'+cell+'.matrix', sep='\t', header=0).drop('loci', axis=1)
features = np.array(X.columns)
X = X.values

# load the targets
y = pd.read_csv(folder + 'train/'+cell+'.labels', sep='\t', header=0)

# customize the targets
def yCustom(y, condition, thresh1=0.0, thresh2=20.0):
    thresh1, thresh2 = float(thresh1), float(thresh2)
    pass_expr = y['Expr'] == 1
    pass_thresh1 = y['Multicov'].values > thresh1
    pass_thresh2 = y['Multicov'].values <= thresh2
    yi = np.array([i for i in y.index 
        if pass_expr[i] and pass_thresh1[i] and pass_thresh2[i]])
    y = y.ix[yi, condition].values
    return y, yi

y, yi = yCustom(y, condition)

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

X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.5, random_state=0)

rl = RandomizedLogisticRegression()
fs = rl.fit(X_train, y_train)
hist_scores = dict(zip(features, fs.scores_))

pickle.dump(hist_scores, open(folder + 'results/histScores'+condition+cell+'.pkl', 'wb'))




