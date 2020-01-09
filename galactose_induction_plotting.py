import multisizer_import as g
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import cPickle as pickle

bp = '/scratch/lab/multisizer/'
columns = ['Mean (fL)', 'Standard Deviation (fL)', 'CV', 'Mode Volume (fL)', 'Median Volume (fL)', \
           'Date', 'Condition', 'Carbon Source', 'Galactose Concentration', 'Strain']
colors=sns.color_palette()
# print colors
full_params =[]
# note that we have to import experiments done on separate days separately because the labeling is often
# inconsistent

# 0-16
params = [['190820', 'yFB29','2% Raff', '0uMGal', 1],    ['190820', 'yFB29','2% Raff', '0uMGal', 2],
          ['190820', 'yFB29', '2% Raff', '100uMGal', 1], ['190820', 'yFB29', '2% Raff', '100uMGal', 2],
          ['190820', 'yFB29', '2% Raff', '200uMGal', 1], ['190820', 'yFB29', '2% Raff', '200uMGal', 2],
          ['190820', 'yFB29', '2% Raff', '400uMGal', 1], ['190820', 'yFB29', '2% Raff', '400uMGal', 2],
          ['190820', 'yFB29', '2% Raff', '800uMGal', 1], ['190820', 'yFB86', '2% Raff', '0uMGal', 1],
          ['190820', 'yFB86', '2% Raff', '0uMGal', 2],   ['190820', 'yFB86', '2% Raff', '100uMGal', 1],
          ['190820', 'yFB86', '2% Raff', '100uMGal', 2], ['190820', 'yFB86', '2% Raff', '200uMGal', 1],
          ['190820', 'yFB86', '2% Raff', '400uMGal', 1], ['190820', 'yFB86', '2% Raff', '400uMGal', 2],
          ['190820', 'yFB86', '2% Raff', '800uMGal', 1]]
for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_1_0'+str(params[i0][4]))
full_params+=params

# 17-33
params = [['190828', 'yFB29','2% Raff', '0uMGal', 1],    ['190828', 'yFB29','2% Raff', '100uMGal', 1],
          ['190828', 'yFB29', '2% Raff', '200uMGal', 1], ['190828', 'yFB29', '2% Raff', '400uMGal', 1],
          ['190828', 'yFB29', '2% Raff', '800uMGal', 1], ['190829', 'yFB30','2% Raff', '0uMGal', 1],
          ['190828', 'yFB30', '2% Raff', '100uMGal', 1], ['190828', 'yFB30', '2% Raff', '200uMGal', 3],
          ['190828', 'yFB30', '2% Raff', '400uMGal', 2],
          ['190828', 'yFB30', '2% Raff', '800uMGal', 2], ['190828', 'yFB30', '2% Raff', '800uMGal', 1],
          ['190828', 'yFB30', '2% Raff', '800uMGal', 3], ['190828', 'yFB86', '2% Raff', '0uMGal', 1],
          ['190828', 'yFB86', '2% Raff', '100uMGal', 3], ['190828', 'yFB86', '2% Raff', '200uMGal', 1],
          ['190828', 'yFB86', '2% Raff', '400uMGal', 1], ['190828', 'yFB86', '2% Raff', '800uMGal', 1]]
for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_1_0'+str(params[i0][4]))
full_params+=params

# 34-50
params = [['191119', 'yFB29','2% Raff', '0uMGal', 1],    ['191119', 'yFB29','2% Raff', '100uMGal', 1],
          ['191119', 'yFB29', '2% Raff', '200uMGal', 1], ['191119', 'yFB29', '2% Raff', '400uMGal', 1],
          ['191119', 'yFB29', '2% Raff', '800uMGal', 1], ['191119', 'yFB30','2% Raff', '0uMGal', 1],
          ['191119', 'yFB30', '2% Raff', '100uMGal', 1], ['191119', 'yFB30', '2% Raff', '200uMGal', 1],
          ['191119', 'yFB30', '2% Raff', '400uMGal', 1],
          ['191119', 'yFB30', '2% Raff', '800uMGal', 1], ['191119', 'yFB30', '2% Raff', '800uMGal', 1],
          ['191119', 'yFB30', '2% Raff', '800uMGal', 1], ['191119', 'yFB86', '2% Raff', '0uMGal', 1],
          ['191119', 'yFB86', '2% Raff', '100uMGal', 1], ['191119', 'yFB86', '2% Raff', '200uMGal', 1],
          ['191119', 'yFB86', '2% Raff', '400uMGal', 1], ['191119', 'yFB86', '2% Raff', '800uMGal', 1]]
for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_1_0'+str(params[i0][4]))
full_params+=params

print len(full_params)
sns.set(font_scale=1.5)
color_guide = {'$P_{GAL1}-WHI5$':0,'$P_{GAL1}-CLN3$':1,'$P_{GAL1}-WHI5$, $\Delta bck2$':0}
# Bar graphs for mutant data
df = pd.DataFrame(columns=columns)
excl_strains = range(22,29,1)+range(39,46,1)
temp_strains = list(set(range(51))-set(excl_strains))
df1,figs = g.import_and_plot_v4(df, full_params, temp_strains, bp, good_copies=True, fs=1.25,lcutoff=10) # just getting the data in pandas format
sns.set(font_scale=1.5)
temp_vals = [0,100,200,400,800,0,100,200,400,800]
xlabels = [str(ind) +' $\mu$M' for ind in temp_vals]
plt.clf()
df1['Genotype'] = [g.strain_db[df1.Strain[i0]] for i0 in temp_strains]
df1['Setup'] = [g.strain_db[df1.Strain[i0]]+' '+ df1['Galactose Concentration'][i0] for i0 in temp_strains]
sns.stripplot(x='Setup', y=g.names['std'],  data=df1, size=5, jitter=False, linewidth=2.0,palette=colors,hue='Genotype')
ax=sns.barplot(x='Setup', y=g.names['std'], data=df1, ci='sd',alpha=0.4,capsize=0.25,palette=colors,dodge=False,hue='Genotype')
# xlabels = [df1['Galactose Concentration'][i0][:-5] +' $\mu$M' for i0 in temp_strains]
# print xlabels
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[2:], labels[2:],loc=2)
ax.set_xticklabels(xlabels, rotation=90)
plt.xlabel('')
plt.ylim(ymax=120,ymin=0)
plt.ylabel(r'$\sigma(V)$ (fL)')
# plt.title(r'$\sigma(V)$, variable Galactose conc.')
fig = ax.get_figure()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/spread_scaling/sd_limited.png',bbox_inches='tight',dpi=300)

plt.clf()
sns.stripplot(x='Setup', y=g.names['med'],  data=df1, size=5, jitter=False, linewidth=2.0,palette=colors,hue='Genotype')
ax=sns.barplot(x='Setup', y=g.names['med'], data=df1, ci='sd',alpha=0.4,capsize=0.25,palette=colors,dodge=False,hue='Genotype')
# xlabels = [df1['Galactose Concentration'][i0][:-5] +' $\mu$M' for i0 in temp_strains]
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[2:], labels[2:],loc=2)
ax.set_xticklabels(xlabels, rotation=90)
plt.xlabel('')
plt.ylim(ymax=140,ymin=0)
plt.ylabel(r'Median volume (fL)')
# plt.title(r'Median Volume, variable Galactose conc.')
fig = ax.get_figure()
# plt.show()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/spread_scaling/med_limited.png',bbox_inches='tight',dpi=300)
plt.clf()
sns.stripplot(x='Setup', y=g.names['mean'],  data=df1, size=5, jitter=False, linewidth=2.0,palette=colors,hue='Genotype')
ax=sns.barplot(x='Setup', y=g.names['mean'], data=df1, ci='sd',alpha=0.4,capsize=0.25,palette=colors,dodge=False,hue='Genotype')
# xlabels = [df1['Galactose Concentration'][i0][:-5] +' $\mu$M' for i0 in temp_strains]
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[2:], labels[2:],loc=2)
ax.set_xticklabels(xlabels, rotation=90)
plt.xlabel('')
plt.ylim(ymax=140,ymin=0)
plt.ylabel(r'$\langle V\rangle$ (fL)')
# plt.title(r'$\langle V\rangle$, variable Galactose conc.')
fig = ax.get_figure()
# plt.show()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/spread_scaling/mean_limited.png',bbox_inches='tight',dpi=300)
plt.clf()
sns.stripplot(x='Setup', y=g.names['cv'],  data=df1, size=5, jitter=False, linewidth=2.0,palette=colors,hue='Genotype')
ax=sns.barplot(x='Setup', y=g.names['cv'], data=df1, ci='sd',alpha=0.4,capsize=0.25,palette=colors,dodge=False,hue='Genotype')
# xlabels = [df1['Galactose Concentration'][i0][:-5] +' $\mu$M' for i0 in temp_strains]
handles, labels = ax.get_legend_handles_labels()
l = plt.legend(handles[2:], labels[2:],loc=2)
ax.set_xticklabels(xlabels, rotation=90)
plt.xlabel('')
plt.ylim(ymax=1.0,ymin=0)
plt.ylabel(r'CV')
# plt.title(r'$CV(V)$, variable Galactose conc.')
fig = ax.get_figure()
# plt.show()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/spread_scaling/cv_limited.png',bbox_inches='tight',dpi=300)
plt.clf()
