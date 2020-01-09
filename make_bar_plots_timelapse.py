import pandas as pd
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pandas.plotting import scatter_matrix
import scipy
from scikits import bootstrap as boot
sns.set(font_scale=1.0)
from pandas.plotting import table
import custom_image_toolkit as c
from scipy.stats import ttest_ind

cell_type_index={'M':'Mothers','D':'Daughters'}
# yFB78 raffinose on 190606
expt_id0 = '/190606_yFB78_60X_Raff_125uMGal'
df0 = pd.read_pickle("./expt_ids"+expt_id0+".pkl")
df0['expt'] = expt_id0
df0['label'] = 'pGAL1-WHI5 06/06'
df0['genotype'] = 'pGAL1-WHI5'
# # yFB79 Raffinose experiment on 181207
expt_id1 = '/181207_yFB79_60X_Raff_125uMGal'
df1 = pd.read_pickle("./expt_ids"+expt_id1+".pkl")
df1['expt'] = expt_id1
df1['label'] = 'pWHI5-WHI5 12/07'
df1['genotype'] = 'pWHI5-WHI5'
# yFB79 Raffinose experiment on 190417
expt_id2 = '/190417_yFB79_60X_Raff_125uMGal'
df2= pd.read_pickle("./expt_ids"+expt_id2+".pkl")
df2['expt'] = '/190417_yFB79_60X_Raff_125uMGal'
df2['label'] = 'pWHI5-WHI5 4/17'
df2['genotype'] = 'pWHI5-WHI5'
# yFB78 expt 190607, 12 min timestep
expt_id3 = '/190607_yFB78_60X_Raff_125uMGal'
df3= pd.read_pickle("./expt_ids"+expt_id3+".pkl")
df3['expt'] = '/190607_yFB78_60X_Raff_125uMGal'
df3['label'] = 'pGAL1-WHI5 6/07'
df3['genotype'] = 'pGAL1-WHI5'
# yFB78 expt 190725
expt_id4 = '/190725_yFB78_60X_2Raff_125uMGal'
df4= pd.read_pickle("./expt_ids"+expt_id4+".pkl")
df4['expt'] = expt_id4
df4['label'] = 'pGAL1-WHI5 7/25'
df4['genotype'] = 'pGAL1-WHI5'
# yFB79 Raffinose experiment on 190612
expt_id6 = '/190612_yFB79_timelapse'
df6= pd.read_pickle("./expt_ids"+expt_id6+".pkl")
df6['expt'] = expt_id6
df6['label'] = 'pWHI5-WHI5 6/12'
df6['genotype'] = 'pWHI5-WHI5'
df0=c.processing(df0)
df1=c.processing(df1)
df=df0.append(df1);
df2=c.processing(df2)
df=df.append(df2);
df3=c.processing(df3)
df=df.append(df3);
df4=c.processing(df4)
df=df.append(df4);
df6=c.processing(df6)
df=df.append(df6);
df = c.normalize(df);
df['Cell Type'] = [cell_type_index[df.iloc[i0].type] for i0 in range(len(df))]
df['Condition'] = [cell_type_index[df.iloc[i0].type]+', '+df.iloc[i0].genotype for i0 in range(len(df))]


def make_bar_plot_cv(input_df,var,ylabel,out_path,ylims=None):
    fig = plt.figure(figsize=[5, 5])
    temp_df = input_df.groupby(['Condition', 'expt'])[var].std()/input_df.groupby(['Condition', 'expt'])[var].mean()
    temp_df1 = temp_df.groupby(['Condition']).mean()
    temp_df2 = temp_df.groupby(['Condition']).std(ddof=1) / np.sqrt(temp_df.groupby(['Condition']).describe()['count'])
    sns.set(font_scale=2)
    temp_df1.plot(kind='bar', yerr=temp_df2, capsize=20.0, figsize=[10, 8], grid=True, legend=False, alpha=0.5)
    # temp_df.reset_index().groupby(['Condition']).plot(kind='bar', x='Condition',y=var)
    # fig=plt.figure(figsize=[8,5])
    g = sns.stripplot(x='Condition', y=var, data=temp_df.reset_index(), size=7, jitter=False, linewidth=2.0)
    labels = ['Daughters,\n $P_{GAL1}-WHI5$', 'Daughters,\n $P_{WHI5}-WHI5$', 'Mothers,\n $P_{GAL1}-WHI5$',
              'Mothers,\n $P_{WHI5}-WHI5$']
    # g.set_xticklabels(labels, rotation=90)
    # g.errorbar(yerr=temp_df1, linewidth=0, capsize=15)
    plt.ylabel(ylabel)
    plt.xlabel('')
    if not(ylims is None):
        plt.ylim(ymin=ylims[0], ymax=ylims[1])
    fig.savefig(out_path, dpi=500, bbox_inches='tight')
    # plt.show()
    plt.clf()



def make_bar_plot_average(input_df,var,ylabel,out_path,ylims=None):
    fig = plt.figure(figsize=[5, 5])
    temp_df = input_df.groupby(['Condition', 'expt'])[var].mean()
    temp_df1 = temp_df.groupby(['Condition']).mean()
    temp_df2 = temp_df.groupby(['Condition']).std(ddof=1) / np.sqrt(temp_df.groupby(['Condition']).describe()['count'])
    sns.set(font_scale=2)
    temp_df1.plot(kind='bar', yerr=temp_df2, capsize=20.0, figsize=[10, 8], grid=True, legend=False, alpha=0.5)
    # temp_df.reset_index().groupby(['Condition']).plot(kind='bar', x='Condition',y=var)
    # fig=plt.figure(figsize=[8,5])
    g = sns.stripplot(x='Condition', y=var, data=temp_df.reset_index(), size=7, jitter=False, linewidth=2.0)
    labels = ['Daughters,\n $P_{GAL1}-WHI5$', 'Daughters,\n $P_{WHI5}-WHI5$', 'Mothers,\n $P_{GAL1}-WHI5$',
              'Mothers,\n $P_{WHI5}-WHI5$']
    # g.set_xticklabels(labels, rotation=90)
    # g.errorbar(yerr=temp_df1, linewidth=0, capsize=15)
    plt.ylabel(ylabel)
    plt.xlabel('')
    if not(ylims is None):
        plt.ylim(ymin=ylims[0], ymax=ylims[1])
    fig.savefig(out_path, dpi=500, bbox_inches='tight')
    # plt.show()
    plt.clf()

# Figure 2 generation

# average volumes

temp_var = "$V_{b,ell}$"
temp_ylabel='Average volume at birth (fL)'
temp_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/Vb_av_strip_timelapse.png'
temp_ylims=[0,110]
make_bar_plot_average(df, temp_var, temp_ylabel, temp_path,temp_ylims)
temp_var = "$V_{s,ell}$"
temp_ylabel='Average volume at Start (fL)'
temp_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/Vs_av_strip_timelapse.png'
temp_ylims=[0,120]
make_bar_plot_average(df, temp_var, temp_ylabel, temp_path,temp_ylims)
temp_var = "$V_{div,ell}$"
temp_ylabel='Average volume at division (fL)'
temp_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/Vd_av_strip_timelapse.png'
temp_ylims=[0,160]
make_bar_plot_average(df, temp_var, temp_ylabel, temp_path,temp_ylims)

# CV in volumes

temp_var = "$V_{b,ell}$"
temp_ylabel='CV in volume at birth'
temp_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/Vb_cv_strip_timelapse.png'
make_bar_plot_cv(df, temp_var, temp_ylabel, temp_path)
temp_var = "$V_{s,ell}$"
temp_ylabel='CV in volume at Start'
temp_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/Vs_cv_strip_timelapse.png'
make_bar_plot_cv(df, temp_var, temp_ylabel, temp_path)
temp_var = "$V_{div,ell}$"
temp_ylabel='CV in volume at division'
temp_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_distribution_perturbation/Vd_cv_strip_timelapse.png'
make_bar_plot_cv(df, temp_var, temp_ylabel, temp_path)

# Figure S5 generation

temp_var = "$c1_{b,seg}$"
temp_ylabel='Average [Whi5] at birth (a.u.)'
temp_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/induction_validation/c1b_av_strip_timelapse.png'
make_bar_plot_average(df, temp_var, temp_ylabel, temp_path)
temp_var = "$c2_{b,seg}$"
temp_ylabel='Average [mCherry] at birth (a.u.)'
temp_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/induction_validation/c2b_av_strip_timelapse.png'
make_bar_plot_average(df, temp_var, temp_ylabel, temp_path)

# figure S7 generation

# average timings

temp_var = "$t_{G1}$"
temp_ylabel='Time in G1 (mins)'
temp_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/microscopy_validation/tG1_av_strip_timelapse.png'
make_bar_plot_average(df, temp_var, temp_ylabel, temp_path)
temp_var = "$t_{budded}$"
temp_ylabel='Time between Start and division (mins)'
temp_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/microscopy_validation/tbudded_av_strip_timelapse.png'
make_bar_plot_average(df, temp_var, temp_ylabel, temp_path)
temp_var = "$t_{div}$"
temp_ylabel='Time from birth to division (mins)'
temp_path = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/microscopy_validation/tdiv_av_strip_timelapse.png'
make_bar_plot_average(df, temp_var, temp_ylabel, temp_path)
