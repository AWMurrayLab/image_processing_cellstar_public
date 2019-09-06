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


sns.set_style("white")

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
# yFB78 raffinose on 190606
expt_id0 = '/190606_yFB78_60X_Raff_125uMGal'
df0 = pd.read_pickle("./expt_ids"+expt_id0+".pkl")
df0['expt'] = expt_id0
df0['label'] = 'pGAL1-WHI5 06/06'
df0['genotype'] = 'pGAL1-WHI5'
# yFB78 expt 190607
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

df0=c.processing(df0)
df1=c.processing(df1)
df=df0.append(df1);
df2=c.processing(df2)
df=df.append(df2);
df3=c.processing(df3)
df=df.append(df3);
df4=c.processing(df4)
df=df.append(df4);
df = c.normalize(df);

bp = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/size_regulation_perturbation/'
titles={'pGAL1-WHI5':r'$P_{GAL1}$ perturbed cells','pWHI5-WHI5':r'$P_{WHI5}$ unperturbed cells'}
#-------------------------------------------
y1,z1='D','pWHI5-WHI5'
temp_xv,temp_yv,save_ext='$V_{b,ell}$','$t_{G1}$','tg1_vb'
temp_xlab, temp_ylab = 'Volume at birth (fL)', 'Time in G1 (mins)'
bin_range=(15,45)
y,z,save_path=df.type==y1,df.genotype==z1,bp+y1+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df[y&z],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],temp_cmap='Oranges')
fig.savefig(save_path,bbox_inches='tight',dpi=500)
xlim = plt.xlim()  # we want each comparison plot to look the same so we get the axes from the first and give it to the second
ylim = plt.ylim()
plt.clf()
#-------------------------------------------
y1,z1='D','pGAL1-WHI5'
temp_xv,temp_yv,save_ext='$V_{b,ell}$','$t_{G1}$','tg1_vb'
temp_xlab, temp_ylab = 'Volume at birth (fL)', 'Time in G1 (mins)'
bin_range=(15,45)
y,z,save_path=df.type==y1,df.genotype==z1,bp+y1+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df[y&z],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],temp_cmap='Blues', temp_xlim=xlim, temp_ylim=ylim,gridsize=[20,12])
fig.savefig(save_path,bbox_inches='tight',dpi=500)
plt.clf()
#-------------------------------------------
y1,z1='D','pWHI5-WHI5'
temp_xv,temp_yv,save_ext='$V_{s,ell}$','$t_{budded}$','tbud_vs'
temp_xlab, temp_ylab = 'Volume at Start (fL)', 'Budded duration (mins)'
bin_range=(30,90)
y,z,save_path=df.type==y1,df.genotype==z1,bp+y1+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df[y&z],bin_range,xlab=temp_xlab,ylab=temp_ylab,gridsize=12,temp_title=titles[z1],temp_cmap='Oranges')
fig.savefig(save_path,bbox_inches='tight',dpi=500)
xlim = plt.xlim()  # we want each comparison plot to look the same so we get the axes from the first and give it to the second
ylim = plt.ylim()
plt.clf()
#-------------------------------------------
y1,z1='D','pGAL1-WHI5'
temp_xv,temp_yv,save_ext='$V_{s,ell}$','$t_{budded}$','tbud_vs'
temp_xlab, temp_ylab = 'Volume at Start (fL)', 'Budded duration (mins)'
bin_range=(30,90)
y,z,save_path=df.type==y1,df.genotype==z1,bp+y1+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df[y&z],bin_range,xlab=temp_xlab,ylab=temp_ylab,gridsize=[15,27],temp_title=titles[z1],temp_cmap='Blues', temp_xlim=xlim, temp_ylim=ylim)
fig.savefig(save_path,bbox_inches='tight',dpi=500)
plt.clf()
#-------------------------------------------
y1,z1='D','pWHI5-WHI5'
temp_xv,temp_yv,save_ext='$V_{b,ell}$','$V_{d,ell}$','vb_vd'
temp_xlab, temp_ylab = 'Volume at birth (fL)', 'Volume at division (fL)'
bin_range=(15,45)
y,z,save_path=df.type==y1,df.genotype==z1,bp+y1+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df[y&z],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],temp_cmap='Oranges')
fig.savefig(save_path,bbox_inches='tight',dpi=500)
xlim = plt.xlim()  # we want each comparison plot to look the same so we get the axes from the first and give it to the second
ylim = plt.ylim()
plt.clf()
#-------------------------------------------
y1,z1='D','pGAL1-WHI5'
temp_xv,temp_yv,save_ext='$V_{b,ell}$','$V_{d,ell}$','vb_vd'
temp_xlab, temp_ylab = 'Volume at birth (fL)', 'Volume at division (fL)'
bin_range=(15,45)
y,z,save_path=df.type==y1,df.genotype==z1,bp+y1+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df[y&z],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],temp_cmap='Blues', temp_xlim=xlim, temp_ylim=ylim,gridsize=[20,8])
fig.savefig(save_path,bbox_inches='tight',dpi=500)
plt.clf()

bp = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/correlations_perturbation/'
#-------------------------------------------
y1,z1='D','pWHI5-WHI5'
temp_xv,temp_yv,save_ext='$V_{b,ell}$','$c1_{b,seg,norm}$','c1b_vb'
temp_xlab, temp_ylab = 'Volume at birth (fL)', '[Whi5] (a.u.)'
bin_range=(15,45)
y,z,save_path=df.type==y1,df.genotype==z1,bp+y1+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df[y&z],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],temp_cmap='Oranges')
fig.savefig(save_path,bbox_inches='tight',dpi=500)
xlim = plt.xlim()  # we want each comparison plot to look the same so we get the axes from the first and give it to the second
ylim = plt.ylim()
plt.clf()
#-------------------------------------------
y1,z1='D','pGAL1-WHI5'
temp_xv,temp_yv,save_ext='$V_{b,ell}$','$c1_{b,seg,norm}$','c1b_vb'
temp_xlab, temp_ylab = 'Volume at birth (fL)', '[Whi5] (a.u.)'
bin_range=(15,45)
y,z,save_path=df.type==y1,df.genotype==z1,bp+y1+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df[y&z],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],temp_cmap='Blues', temp_xlim=xlim, temp_ylim=ylim,gridsize=[27,25])
fig.savefig(save_path,bbox_inches='tight',dpi=500)
plt.clf()
#-------------------------------------------
y1,z1='D','pWHI5-WHI5'
temp_xv,temp_yv,save_ext='$V_{b,ell}$','$c2_{b,seg,norm}$','c2b_vb'
temp_xlab, temp_ylab = 'Volume at birth (fL)', '[mCherry] (a.u.)'
bin_range=(15,45)
y,z,save_path=df.type==y1,df.genotype==z1,bp+y1+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df[y&z],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],temp_cmap='Oranges')
fig.savefig(save_path,bbox_inches='tight',dpi=500)
xlim = plt.xlim()  # we want each comparison plot to look the same so we get the axes from the first and give it to the second
ylim = plt.ylim()
plt.clf()
#-------------------------------------------
y1,z1='D','pGAL1-WHI5'
temp_xv,temp_yv,save_ext='$V_{b,ell}$','$c2_{b,seg,norm}$','c2b_vb'
temp_xlab, temp_ylab = 'Volume at birth (fL)', '[mCherry] (a.u.)'
bin_range=(15,45)
y,z,save_path=df.type==y1,df.genotype==z1,bp+y1+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df[y&z],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],temp_cmap='Blues', temp_xlim=xlim, temp_ylim=ylim,gridsize=[18,4])
fig.savefig(save_path,bbox_inches='tight',dpi=500)
plt.clf()