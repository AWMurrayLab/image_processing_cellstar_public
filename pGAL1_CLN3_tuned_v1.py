import multisizer_import as g
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cPickle as pickle
import scipy
from scipy import stats

bp = '/scratch/lab/multisizer/'
columns = ['Mean (fL)', 'Standard Deviation (fL)', 'CV', 'Mode Volume (fL)', 'Median Volume (fL)', \
           'Date', 'Condition', 'Carbon Source', 'Galactose Concentration', 'Strain']
colors=sns.color_palette()
# print colors
full_params =[]
# note that we have to import experiments done on separate days separately because the labeling is often
# inconsistent

# 0
params = [['190820', 'yFB86', '2% Raff', '200uMGal', 1]]
for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_1_0'+str(params[i0][4]))
full_params+=params

# 1
params = [['190828', 'yFB86', '2% Raff', '200uMGal', 1]]
for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_1_0'+str(params[i0][4]))
full_params+=params

# 2
params = [['191119', 'yFB86', '2% Raff', '200uMGal', 1]]
for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_1_0'+str(params[i0][4]))
full_params+=params

# 3
params = [['180531', 'yFB43','2% Raff', '125uMGal', 1]]
for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_'+str(params[i0][4])+'_01')
full_params+=params
# 4
params = [['190807', 'yFB43','2% Raff', '125uMGal', 1]]
for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_1_0'+str(params[i0][4]))
full_params+=params
# 5
params = [['191003', 'yFB43','2% Raff', '125uMGal', 1]]
for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_1_0'+str(params[i0][4]))
full_params+=params
sns.set(font_scale=1.5)

print len(full_params)
sns.set(font_scale=2)
strain_labels = {'yFB43':'WT', 'yFB86':'$P_{GAL1}-CLN3$, $200\mu M$ Gal'}
# color_guide = {'$P_{GAL1}-WHI5$':0,'$P_{GAL1}-CLN3$':1,'$P_{GAL1}-WHI5$, $\Delta bck2$':0}
# PLOTTING THE AVERAGED STATISTICS OF RELEVANT STRAINS FOR REPEATS OF EACH EXPT
df = pd.DataFrame(columns=columns)
temp_strains = range(len(full_params))
temp_strains=[2,4]
df1,fig = g.import_and_plot_v5(df, full_params, temp_strains, bp,strain_labels=strain_labels) # just getting the data in pandas format
sns.set(font_scale=3)
# plt.show(fig)
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/spread_scaling/CLN3_tuned_CDF.png',dpi=300,bbox_inches='tight')