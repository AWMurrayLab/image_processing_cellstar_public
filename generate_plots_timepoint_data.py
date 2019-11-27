import pandas as pd
import numpy as np
import cPickle as pickle
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from scikits import bootstrap as boot
import scipy
import custom_image_toolkit as c

strain_db={'yFB29':r'pGAL1-WHI5-mVenNB', 'yFB30':r'pGAL1-WHI5-mVenNB, $\Delta$bck2',
           'yFB41':r'pWHI5-WHI5-mVenNB', 'yFB43':r'pWHI5-WHI5-mVenNB', 'yFB25':r'WT', 'yFB86':r'pGAL1-CLN3',
          'yFB45':r'pWHI5-WHI5-mVenNB, $\Delta$bck2', 'yFB46':r'pWHI5-WHI5-mVenNB, $\Delta$bck2',
           'yFB78':r'pGAL1-WHI5-mVenNB','yFB79':r'pWHI5-WHI5-mVenNB',
          'yFB93':r'WT MATa', 'yFB94':r'WT MAT$\alpha$', 'yFB95':r'WT MAT$\alpha$ -leu',
          'yFB96':r'$\Delta$whi5', 'yFB97':r'$\Delta$cln3', 'yFB98':'$\Delta$bck2',
          'yFB99':r'$\Delta$swe1', 'yFB100':r'$\Delta$cln3, $\Delta$whi5',
           'yFB101':r'$\Delta$whi5, $\Delta$bck2',
           'yFB102':r'$\Delta$whi5, $\Delta$cln3, $\Delta$bck2',
           'yFB103':r'$\Delta$whi5, $\Delta$bck2, $\Delta$swe1',
           'yFB104':r'$\Delta$whi5, $\Delta$cln3, $\Delta$bck2, $\Delta$swe1',
          'yFB108':r'$\Delta$whi5, $\Delta$cln3, $\Delta$swe1'}

expt_ids = ['/190403_timepoint', '/190417_timepoint', '/190607_timepoint', '/190322_timepoint']
for ind in range(len(expt_ids)):
    expt_id = expt_ids[ind]
    pickle_in = open("./expt_ids" + expt_id + '.pkl', "rb")
    temp_df = pickle.load(pickle_in)
    temp_df['expt_id'] = expt_id
    if ind == 0:
        df1 = temp_df
    else:
        df1 = df1.append(temp_df)

df1['pixel_thresh_fluor_vals_av'] = df1['pixel_thresh_fluor_vals'] / df1['pixel_thresh_vol']
df1['pixel_thresh_fluor_vals_av_c2'] = df1['pixel_thresh_fluor_vals_c2'] / df1['pixel_thresh_vol']
df1['nucl_cyt_fluor_ratio'] = df1['nuclear_fluor_int'] / df1['cytoplasmic_fluor_int']
df1['nucl_cyt_fluor_ratio_c2'] = df1['nuclear_fluor_int_c2'] / df1['cytoplasmic_fluor_int_c2']
df1['nucl_cyt_vol_ratio'] = df1['nuclear_vol'] / (df1['pixel_thresh_vol'] - df1['nuclear_vol'])
df1['zproj_fluor_vals_conc'] = df1['zproj_fluor_vals'] / df1['ellipse_volume']
df1['strain_num'] = [df1.iloc[i0].Strain[:5] for i0 in range(len(df1))]
df1['gal_conc'] = [df1.iloc[i0].Strain[6:] for i0 in range(len(df1))]
df1['genotype'] = [strain_db[df1.iloc[i0].strain_num] for i0 in range(len(df1))]
df1['Condition'] = [df1.iloc[i0].genotype[:10] + ', ' + df1.iloc[i0].gal_conc for i0 in range(len(df1))]
df1['Condition1'] = [df1.iloc[i0].genotype + ', ' + df1.iloc[i0].gal_conc for i0 in range(len(df1))]


# this normalizes the fluorescence values relative to WT for each condition so that it is easier to see the
# differences

x=df1.expt_id =='/190322_timepoint'
y=df1.nuclear_whi5==1
z=df1.strain_num=='yFB79'
df1.loc[x&y,'pixel_thresh_fluor_vals_av'] = df1[x&y].pixel_thresh_fluor_vals_av/np.mean(df1[x&y&z].pixel_thresh_fluor_vals_av)
df1.loc[x&y,'pixel_thresh_fluor_vals_av_c2'] = df1[x&y].pixel_thresh_fluor_vals_av_c2/np.mean(df1[x&y&z].pixel_thresh_fluor_vals_av_c2)
x=df1.expt_id !='/190322_timepoint'
df1.loc[x&y,'pixel_thresh_fluor_vals_av'] = df1.loc[x&y].pixel_thresh_fluor_vals_av/np.mean(df1[x&y&z].pixel_thresh_fluor_vals_av)
df1.loc[x&y,'pixel_thresh_fluor_vals_av_c2'] = df1[x&y].pixel_thresh_fluor_vals_av_c2/np.mean(df1[x&y&z].pixel_thresh_fluor_vals_av_c2)

bp = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/correlations_perturbation/'
titles={'pGAL1-WHI5-mVenNB':r'$P_{GAL1}$ perturbed cells','pWHI5-WHI5-mVenNB':r'$P_{WHI5}$ unperturbed cells'}
grid=[]
#-------------------------------------------
y1, w1,z1=1,'125uMGal','pWHI5-WHI5-mVenNB'
x1 = df1.pixel_thresh_fluor_vals_av<4.0
x2 = df1.ellipse_volume<200
w11 = df1.gal_conc == w1
temp_xv,temp_yv,save_ext='ellipse_volume','pixel_thresh_fluor_vals_av','c1_vg1_timepoint'
temp_xlab, temp_ylab = 'G1 volume (fL)', '[Whi5] (a.u.)'
bin_range=(20,120)
y,z,save_path=df1.nuclear_whi5==y1,df1.genotype==z1,bp+str(y1)+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df1[y&z&x1&x2&w11],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],gridsize=30,temp_cmap='Oranges')
plt.xlim(xmin=0,xmax=140)
plt.ylim(ymin=0.2,ymax=2.0)
fig.savefig(save_path,bbox_inches='tight',dpi=500)
xlim = plt.xlim()  # we want each comparison plot to look the same so we get the axes from the first and give it to the second
ylim = plt.ylim()
plt.clf()
#-------------------------------------------
y1, w1,z1=1,'125uMGal','pGAL1-WHI5-mVenNB'
x1 = df1.pixel_thresh_fluor_vals_av<4.0
x2 = df1.ellipse_volume<200
w11 = df1.gal_conc == w1
temp_xv,temp_yv,save_ext='ellipse_volume','pixel_thresh_fluor_vals_av','c1_vg1_timepoint'
temp_xlab, temp_ylab = 'G1 volume (fL)', '[Whi5] (a.u.)'
bin_range=(20,120)
y,z,save_path=df1.nuclear_whi5==y1,df1.genotype==z1,bp+str(y1)+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df1[y&z&x1&x2&w11],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],gridsize=[40,30],temp_cmap='Blues', temp_xlim=xlim, temp_ylim=ylim)
plt.xlim(xmin=0,xmax=140)
plt.ylim(ymin=0.2,ymax=2.0)
fig.savefig(save_path,bbox_inches='tight',dpi=500)
plt.clf()
#-------------------------------------------
y1, w1,z1=1,'125uMGal','pWHI5-WHI5-mVenNB'
x1 = df1.pixel_thresh_fluor_vals_av<4.0
x2 = df1.ellipse_volume<200
w11 = df1.gal_conc == w1
temp_xv,temp_yv,save_ext='ellipse_volume','pixel_thresh_fluor_vals_av_c2','c2_vg1_timepoint'
temp_xlab, temp_ylab = 'G1 volume (fL)', '[mCherry] (a.u.)'
bin_range=(20,120)
y,z,save_path=df1.nuclear_whi5==y1,df1.genotype==z1,bp+str(y1)+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df1[y&z&x1&x2&w11],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],gridsize=30,temp_cmap='Oranges')
plt.xlim(xmin=0,xmax=140)
plt.ylim(ymin=0.2,ymax=2.0)
fig.savefig(save_path,bbox_inches='tight',dpi=500)
xlim = plt.xlim()  # we want each comparison plot to look the same so we get the axes from the first and give it to the second
ylim = plt.ylim()
plt.clf()
#-------------------------------------------
y1, w1,z1=1,'125uMGal','pGAL1-WHI5-mVenNB'
x1 = df1.pixel_thresh_fluor_vals_av<4.0
x2 = df1.ellipse_volume<200
w11 = df1.gal_conc == w1
temp_xv,temp_yv,save_ext='ellipse_volume','pixel_thresh_fluor_vals_av_c2','c2_vg1_timepoint'
temp_xlab, temp_ylab = 'G1 volume (fL)', '[mCherry] (a.u.)'
bin_range=(20,120)
y,z,save_path=df1.nuclear_whi5==y1,df1.genotype==z1,bp+str(y1)+'_'+z1[:5]+'_'+save_ext+'.png'
fig=c.plotting_heatmap(temp_xv,temp_yv,df1[y&z&x1&x2&w11],bin_range,xlab=temp_xlab,ylab=temp_ylab,temp_title=titles[z1],gridsize=[45,15],temp_cmap='Blues', temp_xlim=xlim, temp_ylim=ylim)
plt.xlim(xmin=0,xmax=140)
plt.ylim(ymin=0.2,ymax=2.0)
fig.savefig(save_path,bbox_inches='tight',dpi=500)
plt.clf()
#-------------------------------------------