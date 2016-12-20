# load modules
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.optimize import curve_fit
import cPickle as pickle

plt.style.use('ggplot')

# specificity vs expression
def center_data(data):
    return (data-data.mean())/data.std()

def exponential(x, a, b, c):
    return a * np.exp(-b * x) + c

cell = 'Gm12878'
labels = pd.read_csv('predict3/train/'+cell+'.labels', sep='\t', header=0)
spec = center_data(labels['Specificity'].values)
expr = center_data(labels['Multicov'].values)

popt, pcov = curve_fit(exponential, spec, expr)

plt.figure()
plt.scatter(center_data(spec), center_data(expr), alpha=0.1, color='k')
plt.plot(np.sort(center_data(spec)), exponential(np.sort(center_data(spec)), *popt), linewidth=2)
plt.xticks(fontsize=16, color='k')
plt.yticks(fontsize=16, color='k')
plt.xlabel('Centered Specificity', fontsize=20, color='k')
plt.ylabel('Centered CAGE Tags', fontsize=20, color='k')
plt.savefig('figures/explore_specexpr.png', bbox_inches='tight')

print 'R2 =', pearsonr(spec, expr)[0]**2

















































































































































































































































































































# histogram of expression
cell = 'Gm12878'
y = pd.read_csv('predict3/train/'+cell+'.labels', sep='\t', header=0)

plt.figure()
hist = plt.hist(y['Multicov'].values, bins=100, color='k')
plt.xticks(fontsize=16, color='k')
plt.yscale('log')
plt.yticks(fontsize=16, color='k')
plt.xlabel('CAGE Tags', fontsize=20, color='k')
plt.ylabel('Frequency', fontsize=20, color='k')
plt.savefig('figures/explore_cagedistr.png', bbox_inches='tight')

