# change environment
import sys
sys.path.insert(1, '/ifs/home/pw801/bin/python')

# load modules
import pandas as pd
import numpy as np
import scipy
from scipy.stats import linregress
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
from sklearn.metrics import *
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import cPickle as pickle

plt.style.use('ggplot')

try:
    folder, cell, threshold = sys.argv[1:4]
except (IndexError, ValueError):
    folder, cell, threshold = 'predict1/', 'Gm12878', '0.0'

threshold = float(threshold)

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

y = {}
y['full'] = spec
y['res'] = residues

# regression
def regression(X, y):
    rgr = LinearRegression()
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.5, random_state=0)
    y_pred = rgr.fit(X_train, y_train).predict(X_test)
    return y_test, y_pred

def get_coefs(X, y):
    rgr = LinearRegression()
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.5, random_state=0)
    coefs = rgr.fit(X_train, y_train).coef_
    return coefs

for ttype in ['full', 'res']:
    y_test, y_pred = regression(X, y[ttype])
    coefs = get_coefs(X, y[ttype])
    
    slope, intercept, r_value, p_value, stderr = linregress(y_test, y_pred)
    
    """
    f = plt.figure()
    ax = f.add_subplot(111)
    plt.scatter(y_test, y_pred, alpha=0.1, color='k')
    plt.xticks(fontsize=16, color='k')
    plt.yticks(fontsize=16, color='k')
    if ttype == 'res':
        label_note = 'Residuals'
        plt.xlim([-1.1, 1.6])
        plt.ylim([-1.1, 1.6])
    else:
        label_note = 'Specificity'
        plt.xlim([-0.1, 1.1])
        plt.ylim([-0.1, 1.1])
    plt.xlabel(label_note, fontsize=20, color='k')
    plt.ylabel('Predicted '+label_note, fontsize=20, color='k')
    plt.text(0.05, 0.95, 'r = {0:.2f}'.format(r_value), fontsize=12, fontweight='semibold', transform=ax.transAxes)
    plt.savefig(folder + 'figures/rgrSpec_'+ttype+'_'+cell+'_'+str(threshold)+'.png', bbox_inches='tight')
    plt.clf()
    """
    
    saveres = {'r_value':r_value, 'slope':slope, 'intercept':intercept, 
        'p_value':p_value, 'stderr':stderr, 'coefs':coefs}
    pickle.dump(saveres, open(folder + 'results/rgrSpec_'+ttype+'_'+cell+'_'+str(threshold)+'.pkl', 'wb'))

print 'Done Spec', folder+cell


