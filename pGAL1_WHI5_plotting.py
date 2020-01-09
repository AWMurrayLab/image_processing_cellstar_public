import multisizer_import as g
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
# import cPickle as pkl

bp = '/scratch/lab/multisizer/'
columns = ['Mean (fL)', 'Standard Deviation (fL)', 'CV', 'Mode Volume (fL)', 'Median Volume (fL)', \
           'Date', 'Condition', 'Carbon Source', 'Galactose Concentration', 'Strain']


full_params =[]
# note that we have to import experiments done on separate days separately because the labeling is often
# inconsistent

# 0-6
params = [['180522', 'yFB29','2% Raff', '0uMGal', 1], ['180531', 'yFB29','2% Raff', '100uMGal', 1],
          ['180531', 'yFB29','2% Raff', '125uMGal', 1], ['180531', 'yFB29','2% Raff', '150uMGal', 1],
          ['180531', 'yFB30','2% Raff', '125uMGal', 1], ['180531', 'yFB43','2% Raff', '125uMGal', 1],
          ['180531', 'yFB46','2% Raff', '125uMGal', 1]]
for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_'+str(params[i0][4])+'_01')
full_params+=params
# 7-10
params = [['180725', 'yFB29','2% Raff', '800uMGal', 1], ['180725', 'yFB30','2% Raff', '125uMGal', 1], \
          ['180725', 'yFB41','2% Raff', '800uMGal', 1], ['180725', 'yFB45','2% Raff', '125uMGal', 1]]
for i0 in range(len(params)):
    params[i0].append(params[i0][0]+'_'+params[i0][1]+'_'+params[i0][3]+'_'+str(params[i0][4])+'_01')
full_params+=params
# 11-15
params = [['180524', 'yFB29','2% Raff', '100uMGal', 1], ['180524', 'yFB29','2% Raff', '125uMGal', 1], \
          ['180524', 'yFB29', '2% Raff', '125uMGal', 2], ['180524', 'yFB29','2% Raff', '150uMGal', 1], \
          ['180524', 'yFB29', '2% Raff', '150uMGal', 2]]
for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3][:3]+'g_V'+str(params[i0][4])+'_01')
full_params+=params
# 16-23
params = [['190807', 'yFB29','2% Raff', '125uMGal', 1], ['190807', 'yFB30','2% Raff', '125uMGal', 1],
          ['190807', 'yFB43','2% Raff', '125uMGal', 1], ['190807', 'yFB78', '2% Raff', '0uMGal', 1],
          ['190807', 'yFB46', '2% Raff', '125uMGal', 2], ['190807', 'yFB43', '2% Raff', '125uMGal', 2],
          ['190807', 'yFB29', '2% Raff', '125uMGal', 2], ['190807', 'yFB29', '2% Raff', '125uMGal', 3]]

for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_1_0'+str(params[i0][4]))
full_params+=params
# 24-27
params = [['191003', 'yFB29','2% Raff', '125uMGal', 1], ['191003', 'yFB30','2% Raff', '125uMGal', 1],
          ['191003', 'yFB43','2% Raff', '125uMGal', 1], ['191003', 'yFB46', '2% Raff', '125uMGal', 1]]

for i0 in range(len(params)):
    params[i0].append(params[i0][1]+'_'+params[i0][3]+'_1_0'+str(params[i0][4]))
full_params+=params
sns.set(font_scale=1.5)

# PLOTTING THE AVERAGED STATISTICS OF RELEVANT STRAINS FOR REPEATS OF EACH EXPT
df = pd.DataFrame(columns=columns)
temp_strains = [5, 18, 26, 2, 22, 24, 6, 20, 27, 4, 17, 25]
df1,figs = g.import_and_plot_v4(df, full_params, temp_strains, bp, good_copies=True, fs=1.25, label = r'2% Raff, $125\mu$M Galactose',lcutoff=10.0) # just getting the data in pandas format
plt.clf()
sns.set(font_scale=1.0)
df1['Genotype'] = [g.strain_db[df1.Strain[i0]] for i0 in temp_strains]
temp_df= df1.groupby(['Genotype', 'Date'])[g.names['cv']].mean()
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
temp_df= df1.groupby(['Genotype', 'Date'])[g.names['mean']].mean()
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
temp_df= df1.groupby(['Genotype', 'Date'])[g.names['std']].mean()
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
sns.set(font_scale=2)
plot = sns.barplot(x='Genotype', y=g.names['std'], data=df1, ci='sd',alpha=0.5,capsize=0.25)
ax = sns.stripplot(x='Genotype', y=g.names['std'], data=df1, size=7, jitter=False, linewidth=2.0)
labels = ['$P_{WHI5}-WHI5$', '$P_{GAL1}-WHI5$','$P_{WHI5}-WHI5$,\n$\Delta bck2$', '$P_{GAL1}-WHI5$,\n$\Delta bck2$']
# ax.set_xticklabels(labels, rotation=90)
plt.xticks(rotation=90)
plt.xlabel('')
plt.ylabel(r'Standard Deviation $\sigma_V$ (fL)')
# plt.title(r'Standard Deviation, 125$\mu$M Galactose')
fig = plot.get_figure()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/pGAL1_WHI5_std_tuned.png',bbox_inches='tight',dpi=500)
plt.clf()
plot = sns.barplot(x='Genotype', y=g.names['med'], data=df1, ci='sd',alpha=0.4)
plt.xlabel('Genotype')
plt.xticks(rotation=90)
plt.ylabel(r'Median cell volume $\sigma_V$ (fL)')
# plt.title(r'Median Volume, 2% Raffinose, 125$\mu$M Galactose')
fig = plot.get_figure()
# plt.show()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/pGAL1_WHI5_med_tuned.png',bbox_inches='tight',dpi=500)
plt.clf()
plot = sns.barplot(x='Genotype', y=g.names['mean'], data=df1, ci='sd',alpha=0.5,capsize=0.25)
ax = sns.stripplot(x='Genotype', y=g.names['mean'], data=df1, size=7, jitter=False, linewidth=2.0)
labels = ['$P_{WHI5}-WHI5$', '$P_{GAL1}-WHI5$','$P_{WHI5}-WHI5$,\n$\Delta bck2$', '$P_{GAL1}-WHI5$,\n$\Delta bck2$']
# ax.set_xticklabels(labels, rotation=90)
plt.xticks(rotation=90)
plt.xlabel('')
plt.ylabel(r'Average cell volume (fL)')
# plt.title(r'Average Volume, 125$\mu$M Galactose')
fig = plot.get_figure()
# plt.show()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/pGAL1_WHI5_av_tuned.png',bbox_inches='tight',dpi=500)
plt.clf()
plot = sns.barplot(x='Genotype', y=g.names['cv'], data=df1, ci='sd',alpha=0.5,capsize=0.25)
ax = sns.stripplot(x='Genotype', y=g.names['cv'], data=df1, size=7, jitter=False, linewidth=2.0)
labels = ['$P_{WHI5}-WHI5$', '$P_{GAL1}-WHI5$','$P_{WHI5}-WHI5$,\n$\Delta bck2$', '$P_{GAL1}-WHI5$,\n$\Delta bck2$']
# ax.set_xticklabels(labels, rotation=90)
plt.xticks(rotation=90)
plt.xlabel('')
plt.ylabel(r'CV in cell volume')
# plt.title(r'CV in cell volume, 125$\mu$M Galactose')
fig = plot.get_figure()
# plt.show()
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/pGAL1_WHI5_cv_tuned.png',bbox_inches='tight',dpi=500)
plt.clf()

df1.to_pickle('./multisizer_data/pGAL1_WHI5_mCherry_tuned.pkl')
sns.set(font_scale=1.5)
# plotting CDF for each paired data set
df = pd.DataFrame(columns=columns)
temp_strains = [5,2,6,4]
fig1, df1 = g.import_and_plot_v1(df, full_params, temp_strains, bp, good_copies=True, fs=1.7, title = r'2% Raff, $125\mu$M Galactose',cutoff=2.0) # just getting the data in pandas format
fig1.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/pGAL1_WHI5_CDF_tuned_1.png',bbox_inches='tight',dpi=500)
df = pd.DataFrame(columns=columns)
temp_strains = [18,22,20,17]
fig1, df1 = g.import_and_plot_v1(df, full_params, temp_strains, bp, good_copies=True, fs=1.7, title = r'2% Raff, $125\mu$M Galactose',cutoff=2.0) # just getting the data in pandas format
fig1.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/pGAL1_WHI5_CDF_tuned_2.png',bbox_inches='tight',dpi=500)
df = pd.DataFrame(columns=columns)
temp_strains = [26,24,27,25]
fig1, df1 = g.import_and_plot_v1(df, full_params, temp_strains, bp, good_copies=True, fs=1.7, title = r'2% Raff, $125\mu$M Galactose',cutoff=2.0) # just getting the data in pandas format
fig1.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/pGAL1_WHI5_CDF_tuned_3.png',bbox_inches='tight',dpi=500)

# temp_strains = [5, 9, 18, 2, 13, 22, 6, 10, 20, 4, 8, 17]

#
# # trying again
# sns.set(font_scale=1.5)
# df1['Genotype'] = [g.strain_db[df1.Strain[i0]] for i0 in temp_strains]
# # print df.Genotype
# plot = sns.barplot(x='Genotype', y=g.names['std'], data=df1, ci='sd',alpha=0.4)
# # plt.ylim(ymin=0.0, ymax=180)
# plt.xlabel('Genotype')
# plt.xticks(rotation=90)
# plt.ylabel(r'$\sigma_V$ (fL)')
# plt.title(r'Standard Deviation $\sigma_V$')
# fig = plot.get_figure()
# # plt.show()
# fig.savefig('./multisizer_plots/pGAL1_WHI5/pGAL1_WHI5_std_limited_1_v1.png',bbox_inches='tight')
# plt.clf()
# plot = sns.barplot(x='Genotype', y=g.names['med'], data=df1, ci='sd',alpha=0.4)
# # plot = df1.plot(x='Genotype',y=g.names['med'], kind='bar',figsize=[15,5],legend=None,alpha=0.6)
# # plt.ylim(ymin=0.0, ymax=180)
# plt.xlabel('Genotype')
# plt.xticks(rotation=90)
# plt.ylabel(r'Median cell volume $(fL)$')
# plt.title(r'Median Volume, 2% Raffinose, 125$\mu$M Galactose')
# fig = plot.get_figure()
# # plt.show()
# fig.savefig('./multisizer_plots/pGAL1_WHI5/pGAL1_WHI5_med_limited_1_v1.png',bbox_inches='tight')
# plt.clf()
# plot = sns.barplot(x='Genotype', y=g.names['cv'], data=df1, ci='sd',alpha=0.4)
# # plot = df1.plot(x='Genotype',y=g.names['med'], kind='bar',figsize=[15,5],legend=None,alpha=0.6)
# # plt.ylim(ymin=0.0, ymax=180)
# plt.xlabel('Genotype')
# plt.xticks(rotation=90)
# plt.ylabel(r'CV in cell volume ')
# plt.title(r'CV in cell volume, 2% Raffinose, 125$\mu$M Galactose')
# fig = plot.get_figure()
# # plt.show()
# fig.savefig('./multisizer_plots/pGAL1_WHI5/pGAL1_WHI5_cv_limited_1_v1.png',bbox_inches='tight')
# plt.clf()
#
#
#
# # PLOTTING THE CORRECTLY TUNED WHI5 EXPRESSION DATA
# sel_strains = [5, 2, 6, 4]
# df = pd.DataFrame(columns=columns)
# df1 = g.import_and_plot_v3(df, full_params, temp_strains, bp, good_copies=True, fs=1.25, title = r'2% Raff, $125\mu$M Galactose') # just getting the data in pandas format
# df1['Genotype'] = [g.strain_db[df1.Strain[i0]] for i0 in temp_strains]
# for ind in range(1,len(sel_strains)+1):
#     temp_strains = sel_strains[:ind]
#     # CDF
#     fig, df = g.import_and_plot_v1(df, full_params, temp_strains, bp, good_copies=True, fs=1.5, title = r'2% Raff, $125\mu$M Galactose')
#     fig.savefig('./multisizer_plots/pGAL1_WHI5/pGAL1_WHI5_CDF_tuned_v{0}.png'.format(ind),bbox_inches='tight')
#     # PDF
#     fig, df = g.import_and_plot_v2(df, full_params, temp_strains, bp, good_copies=True, fs=1.25, title = r'2% Raff, $125\mu$M Galactose')
#     fig.savefig('./multisizer_plots/pGAL1_WHI5/pGAL1_WHI5_PDF_tuned_v{0}.png'.format(ind),bbox_inches='tight')
#
# # PLOTTING THE COMPARISON OF 125UM GAL WHI5 EXPRESSION DATA
# temp_strains = [2, 12, 13]
# df = pd.DataFrame(columns=columns)
# fig, df = g.import_and_plot_v1(df, full_params, temp_strains, bp, good_copies=True, fs=1.0, title = r'2% Raff, $125\mu$M Galactose')
# fig.savefig('./multisizer_plots/pGAL1_WHI5/pGAL1_WHI5_CDF_tuned_comparison_.png',bbox_inches='tight')
#
# ###############################################################################
# # Plotting tagged WHI5 overexpression data
# saving_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/'
# group = 'induction_validation'
# df = pd.DataFrame(columns=columns)
# temp_strains = [2,7,5,9]
# fig, df = g.import_and_plot_v1(df, full_params, temp_strains, bp, good_copies=True, fs=1.5, title = r'2% Raff')
# fig.savefig(saving_path+group+'/pGAL1_WHI5-mVenNB_CDF_tuned_v1.png', bbox_inches='tight')
# plt.clf()
#
# df = pd.DataFrame(columns=columns)
# temp_strains = [2,13,7,5,9]
# fig, df = g.import_and_plot_v1(df, full_params, temp_strains, bp, good_copies=True, fs=1.5, title = r'2% Raff')
# df1 = g.import_and_plot_v3(df, full_params, temp_strains, bp, good_copies=True, fs=1.25, title = r'2% Raff') # just getting the data in pandas format
# sns.set(font_scale=1.5)
# df1['Genotype'] = [g.strain_db[df1.Strain[i0]] for i0 in temp_strains]
# df1['Cond'] = [df1['Genotype'][i0] +', '+ df1['Galactose Concentration'][i0] for i0 in temp_strains]
# plt.clf()
# plot = sns.barplot(x='Cond', y=g.names['med'], data=df1, ci='sd',alpha=0.4)
# plt.xlabel('Condition')
# plt.xticks(rotation=90)
# plt.ylabel(r'Median Volume (fL)')
# plt.title(r'Median cell volume')
# fig = plot.get_figure()
# fig.savefig(saving_path+group+'/pGAL1_Whi5_mVenNB_overexpression_med.png',bbox_inches='tight')
# plt.clf()
#
# plot = sns.barplot(x='Cond', y=g.names['cv'], data=df1, ci='sd',alpha=0.4)
# plt.xlabel('Condition')
# plt.xticks(rotation=90)
# plt.ylabel(r'CV (fL)')
# plt.title(r'CV in volume')
# fig = plot.get_figure()
# fig.savefig(saving_path+group+'/pGAL1_Whi5_mVenNB_overexpression_cv.png',bbox_inches='tight')
# plt.clf()
#
# plot = sns.barplot(x='Cond', y=g.names['std'], data=df1, ci='sd',alpha=0.4)
# plt.xlabel('Condition')
# plt.xticks(rotation=90)
# plt.ylabel(r'Standard Deviation in volume (fL)')
# plt.title(r'Standard Deviation in volume')
# fig = plot.get_figure()
# fig.savefig(saving_path+group+'/pGAL1_Whi5_mVenNB_overexpression_std.png',bbox_inches='tight')
# plt.clf()