import pandas as pd
import numpy as np
import cPickle as pickle
import seaborn as sns
import matplotlib.pyplot as plt
import os

# which experiment we want to plot
# expt_id = '/190403_timepoint'
# expt_id = '/190322_timepoint'
expt_id = '/190417_timepoint'

# loading data
pickle_in = open("./expt_ids"+expt_id+'.pkl',"rb")
df1 = pickle.load(pickle_in)

# additional manipulations
df1['pixel_thresh_fluor_vals_av'] = df1['pixel_thresh_fluor_vals']/df1['pixel_thresh_vol']
df1['pixel_thresh_fluor_vals_av_c2'] = df1['pixel_thresh_fluor_vals_c2']/df1['pixel_thresh_vol']
df1['nucl_cyt_fluor_ratio'] = df1['nuclear_fluor_int']/df1['cytoplasmic_fluor_int']
df1['nucl_cyt_fluor_ratio_c2'] = df1['nuclear_fluor_int_c2']/df1['cytoplasmic_fluor_int_c2']
df1['nucl_cyt_vol_ratio'] =  df1['nuclear_vol']/(df1['pixel_thresh_vol']-df1['nuclear_vol'])

if not os.path.exists('./expt_ids'+expt_id):
    os.makedirs('./expt_ids'+expt_id)

# plotting
sns_plot = sns.lmplot(x='pixel_thresh_vol', y='pixel_thresh_fluor_vals_av', hue='Strain', data=df1, col='nuclear_whi5', scatter_kws={'alpha': 0.3})
plt.ylim(ymax=5000)
sns_plot.savefig('./expt_ids'+expt_id+'/fig_1.png')

sns_plot = sns.lmplot(x='pixel_thresh_vol', y='nuclear_fluor_av', col='Strain', data=df1[df1['nuclear_whi5']==1], scatter_kws={'alpha': 0.3})
sns_plot.savefig('./expt_ids'+expt_id+'/fig_2.png')

sns_plot = sns.catplot(x='Strain', y='pixel_thresh_vol',data=df1[df1['nuclear_whi5']==1],kind='violin')
sns_plot.savefig('./expt_ids'+expt_id+'/fig_3.png')

sns_plot = sns.catplot(x='Strain', y='nucl_cyt_fluor_ratio',data=df1[df1['nuclear_whi5']==1],kind='violin')
plt.ylim(ymax=1.0)
sns_plot.savefig('./expt_ids'+expt_id+'/fig_4.png')

sns_plot = sns.catplot(x='Strain', y='nucl_cyt_vol_ratio',data=df1[df1['nuclear_whi5']==1], kind='violin')
plt.ylim(ymax=0.2)
sns_plot.savefig('./expt_ids'+expt_id+'/fig_5.png')

sns_plot = sns.catplot(x='Strain', y='nucl_cyt_fluor_ratio_c2',data=df1[df1['nuclear_whi5']==1],kind='violin')
sns_plot.savefig('./expt_ids'+expt_id+'/fig_6.png')

sns_plot = sns.lmplot(x='pixel_thresh_vol', y='nuclear_fluor_av_c2', col='Strain', data=df1[df1.nuclear_whi5==1], scatter_kws={'alpha': 0.3})
sns_plot.savefig('./expt_ids'+expt_id+'/fig_7.png')

sns_plot = sns.lmplot(x='pixel_thresh_vol', y='pixel_thresh_fluor_vals_av_c2', hue='Strain', data=df1, col='nuclear_whi5', scatter_kws={'alpha': 0.3})
sns_plot.savefig('./expt_ids'+expt_id+'/fig_8.png')