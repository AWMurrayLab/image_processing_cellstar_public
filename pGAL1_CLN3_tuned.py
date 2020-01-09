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

print len(full_params)

# color_guide = {'$P_{GAL1}-WHI5$':0,'$P_{GAL1}-CLN3$':1,'$P_{GAL1}-WHI5$, $\Delta bck2$':0}
# PLOTTING THE AVERAGED STATISTICS OF RELEVANT STRAINS FOR REPEATS OF EACH EXPT
df = pd.DataFrame(columns=columns)
temp_strains = range(len(full_params))
excl_strains = []
# excl_strains=[25]
temp_strains=list(set(temp_strains)-set(excl_strains))
df1,figs = g.import_and_plot_v4(df, full_params, temp_strains, bp, good_copies=True, fs=1.25,lcutoff=10) # just getting the data in pandas format
sns.set(font_scale=3)
plt.clf()
df1['Genotype'] = [g.strain_db[df1.Strain[i0]] for i0 in temp_strains]
df1['Setup'] = [g.strain_db[df1.Strain[i0]]+' '+ df1['Galactose Concentration'][i0] for i0 in temp_strains]
pickle_in = open("./multisizer_data/pGAL1_WHI5_mCherry_tuned.pkl","rb")
df_temp = pickle.load(pickle_in)
df_temp1 = df_temp[df_temp.Genotype=='$P_{WHI5}-WHI5$']
df_temp1.Setup=['WT', 'WT', 'WT']
df_temp1=df_temp1.replace('$P_{WHI5}-WHI5$', 'WT')

df2=pd.concat([df1,df_temp1],axis=0,sort=True)

temp_df= df2.groupby(['Genotype', 'Date'])[g.names['cv']].mean()
print temp_df.head(20)
labels = temp_df.index.unique().levels[0]
# print labels
# print temp_df.head(20)
print 'P values for '+g.names['cv']
for i in range(len(labels)):
    for j in range(len(labels)):
        if i<j:
            print labels[i], labels[j]
            vals = scipy.stats.ttest_ind(temp_df.loc[labels[i]],temp_df.loc[labels[j]],equal_var=False)
            if vals[1]<0.1:
                print 'Statistically significant difference'
                print np.around(vals,3)
            else:
                print 'No difference'
temp_df= df2.groupby(['Genotype', 'Date'])[g.names['std']].mean()
print temp_df.head(20)
labels = temp_df.index.unique().levels[0]
# print labels
# print temp_df.head(20)
print 'P values for '+g.names['std']
for i in range(len(labels)):
    for j in range(len(labels)):
        if i<j:
            print labels[i], labels[j]
            vals = scipy.stats.ttest_ind(temp_df.loc[labels[i]],temp_df.loc[labels[j]],equal_var=False)
            if vals[1]<0.1:
                print 'Statistically significant difference'
                print np.around(vals,3)
            else:
                print 'No difference'
temp_df= df2.groupby(['Genotype', 'Date'])[g.names['mean']].mean()
print temp_df.head(20)
labels = temp_df.index.unique().levels[0]
# print labels
# print temp_df.head(20)
print 'P values for '+g.names['mean']
for i in range(len(labels)):
    for j in range(len(labels)):
        if i<j:
            print labels[i], labels[j]
            vals = scipy.stats.ttest_ind(temp_df.loc[labels[i]],temp_df.loc[labels[j]],equal_var=False)
            if vals[1]<0.1:
                print 'Statistically significant difference'
                print np.around(vals,3)
            else:
                print 'No difference'
print df2.tail(10)
plot = sns.barplot(x='Genotype', y=g.names['std'], data=df2, ci='sd',alpha=0.4,capsize=0.25,palette=colors,dodge=False)
ax=sns.stripplot(x='Genotype', y=g.names['std'],  data=df2, size=5, jitter=False, linewidth=2.0,palette=colors)
# labels = [df2['Galactose Concentration'][i0][:-5] +' $\mu$M' for i0 in temp_strains]
# plt.legend(loc=[1.0,0.2])
# plot.set_xticklabels(rotation=90)
plt.xlabel('')
plt.ylabel(r'$\sigma(V)$ (fL)')
# plt.title(r'$\sigma(V)$')
fig = plot.get_figure()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/spread_scaling/CLN3_tuned_sd.png',bbox_inches='tight',dpi=300)
print df1.Genotype
plt.clf()
plot = sns.barplot(x='Genotype', y=g.names['med'], data=df2, ci='sd',alpha=0.4,capsize=0.25,dodge=False)
ax=sns.stripplot(x='Genotype', y=g.names['med'],  data=df2, size=5, jitter=False, linewidth=2.0)
# labels = [df1['Galactose Concentration'][i0][:-5] +' $\mu$M' for i0 in temp_strains]
# plt.legend(loc=[1.0,0.2])
# plot.set_xticklabels(labels, rotation=90)
plt.xlabel('')
plt.ylabel(r'Median volume (fL)')
# plt.title(r'Median Volume, variable Galactose conc.')
fig = plot.get_figure()
# plt.show()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/spread_scaling/CLN3_tuned_med.png',bbox_inches='tight',dpi=300)
plt.clf()
plot = sns.barplot(x='Genotype', y=g.names['mean'], data=df2, ci='sd',alpha=0.4,capsize=0.25,palette=colors,dodge=False)
ax=sns.stripplot(x='Genotype', y=g.names['mean'],  data=df2, size=5, jitter=False, linewidth=2.0,palette=colors)
# labels = [df1['Galactose Concentration'][i0][:-5] +' $\mu$M' for i0 in temp_strains]
# plt.legend(loc=[1.0,0.2])
# plot.set_xticklabels(labels, rotation=90)
plt.xlabel('')
plt.ylabel(r'$\langle V\rangle$ (fL)')
# plt.title(r'$\langle V\rangle$, variable Galactose conc.')
fig = plot.get_figure()
# plt.show()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/spread_scaling/CLN3_tuned_mean.png',bbox_inches='tight',dpi=300)
plt.clf()
plot = sns.barplot(x='Genotype', y=g.names['cv'], data=df2, ci='sd',alpha=0.4,capsize=0.25,palette=colors,dodge=False)
ax=sns.stripplot(x='Genotype', y=g.names['cv'],  data=df2, size=5, jitter=False, linewidth=2.0,palette=colors)
# labels = [df1['Galactose Concentration'][i0][:-5] +' $\mu$M' for i0 in temp_strains]
# plt.legend(loc=[1.0,0.2])
# plot.set_xticklabels(labels, rotation=90)
plt.xlabel('')
plt.ylabel(r'CV')
# plt.title(r'$CV(V)$, variable Galactose conc.')
fig = plot.get_figure()
# plt.show()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/spread_scaling/CLN3_tuned_cv.png',bbox_inches='tight',dpi=300)
plt.clf()
