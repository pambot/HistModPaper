import glob
import math as m

nsegs = 15478507 - 2*24 # wc -l <bin files> minus headers
prefix = 'Hist-2'
folder = prefix + '_load/'
files = glob.glob(folder + '*/model_*.txt')

models = []
collect = []
for model in files:
	with open(model, 'r') as f:
		l = f.readline().rstrip('\n').split('\t')
	
	k = int(model.split('_')[-1].split('.')[0])
	logl = float(l[3])
	bic = -2 * logl + k * m.log(nsegs)
	models.append(k)
	collect.append(bic)

best = collect.index(min(collect))
bestm = models[best]
print "The best number of states is", bestm

import matplotlib.pyplot as plt
plt.style.use('ggplot')

f = plt.figure(1, figsize=(5, 5))
ax = f.add_subplot(111)
plt.scatter(models, collect, c='k')
ax.tick_params(axis='both', which='major', labelsize=16, colors='k')
plt.xlabel('Model Size', fontsize=20, color='k')
plt.ylabel('BIC Score', fontsize=20, color='k')
plt.savefig(prefix + '_bic.png', bbox_inches='tight')


