import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import manual_image_toolkit as M
import os
import scipy
from scipy import stats

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

pixel_size = {'60X': 0.267, '100X': 0.16}
drange = 65535.0
# this script takes a list of cell tracks with daughter lineages assigned (e.g. the output from assign_lineages_1.py)
# and creates a collection of reliable cell cycles from this by iterating through this list.
color = cm.tab10(np.linspace(0, 1, 10))

# importing the cell cycle variables

# pACT1-mKate2 experiment on 180910
scale = pixel_size['60X']
base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180910_pACT1_mKate2/timelapse'  # path for the desktop
# storage of data.
# base_path, expt_path = '/mnt/d/Lab/image_analysis', '/180831_yFB29_800uMGal'  # path for the laptop storage of images
script_path = '/mnt/c/Users/felix/Documents/GitHub/image_processing_cellstar'
timestep = 5.0
fluor_c2 = True
zscale = 0.7

# # yFB29 experiment on 180831
# scale = pixel_size['60X']
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180831_yFB29_800uMGal/timelapse'
# timestep = 10.0
# fluor_c2 = False

directory = base_path + expt_path + '/plots'
if not os.path.exists(directory):
    os.makedirs(directory)
with open(base_path+expt_path+'/cell_cycles_filtered.pkl', 'rb') as input:
    filt_cc = pickle.load(input)

fig = plt.figure(figsize=[7, 7])
xv = [obj.vb for obj in filt_cc]
yv = [(obj.zproj_fl[-1]+obj.zproj_fl_bud[-1]-obj.zproj_fl[0])/((len(obj.frames)-1)*timestep) for obj in filt_cc]
temp = scipy.stats.pearsonr(np.asarray(xv)*scale**3, yv)
plt.plot(np.asarray(xv)*scale**3, yv, marker='.', label=' Z sum Fluor. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(filt_cc)), linestyle='None')
# plt.hexbin(np.asarray(xv)*scale**3, yv)
plt.xlabel('$V_b$ ($\mu m^3$)')
plt.ylabel('Fluorescence rate of production')
plt.legend()
plt.title('Fluorescence rate of change')
fig.savefig(directory+'/df_zproj_vb_point.png', bbox_inches='tight', dpi=fig.dpi)
del fig


fig = plt.figure(figsize=[7, 7])
xv = [obj.pixel_thresh_vol[0]*scale**2*zscale for obj in filt_cc]
yv = [(obj.pixel_thresh_fluor_vals[-1]+obj.segment_fl_bud[-1]-obj.pixel_thresh_fluor_vals[0])/((len(obj.frames)-1)*timestep) for obj in filt_cc]
temp = scipy.stats.pearsonr(np.asarray(xv), yv)
plt.plot(np.asarray(xv), yv, marker='.', label=' Pixel thresh Fluor. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(filt_cc)), linestyle='None')
# plt.hexbin(np.asarray(xv)*scale**3, yv)
plt.xlabel('$V_b$ pixel threshold ($\mu m^3$)')
plt.ylabel('Fluorescence rate of production')
plt.legend()
plt.title('Fluorescence rate of change')
fig.savefig(directory+'/df_thresh_vb_thresh_point.png', bbox_inches='tight', dpi=fig.dpi)
del fig


fig = plt.figure(figsize=[7, 7])
xv = [obj.pixel_thresh_vol[0]*scale**2*zscale for obj in filt_cc]
yv = [(obj.pixel_thresh_fluor_vals_c2[-1]+obj.segment_fl_bud_c2[-1]-obj.pixel_thresh_fluor_vals_c2[0])/((len(obj.frames)-1)*timestep) for obj in filt_cc]
temp = scipy.stats.pearsonr(np.asarray(xv), yv)
plt.plot(np.asarray(xv), yv, marker='.', label=' Pixel thresh Fluor. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(filt_cc)), linestyle='None')
# plt.hexbin(np.asarray(xv)*scale**3, yv)
plt.xlabel('$V_b$ pixel threshold ($\mu m^3$)')
plt.ylabel('Fluorescence rate of production')
plt.legend()
plt.title('Fluorescence rate of change')
fig.savefig(directory+'/df_thresh_vb_thresh_point_c2.png', bbox_inches='tight', dpi=fig.dpi)
del fig


fig = plt.figure(figsize=[7, 7])
xv = [obj.pixel_thresh_vol[0]*scale**2*zscale for obj in filt_cc]
yv = [(obj.pixel_thresh_vol[-1]+obj.segment_vol_bud[-1]-obj.pixel_thresh_vol[0])*scale**2*zscale/((len(obj.frames)-1)*timestep) for obj in filt_cc]
temp = scipy.stats.pearsonr(np.asarray(xv), yv)
plt.plot(np.asarray(xv), yv, marker='.', label=' Volume. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(filt_cc)), linestyle='None')
# plt.hexbin(np.asarray(xv)*scale**3, yv, label=' Volume. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 1), len(filt_cc)))
plt.xlabel('$V_b$ ($\mu m^3$)')
plt.ylabel('Volume rate of production')
plt.legend()
plt.title('Volume rate of change')
fig.savefig(directory+'/dv_vb_point_thresh.png', bbox_inches='tight', dpi=fig.dpi)
del fig


fig = plt.figure(figsize=[7, 7])
xv = np.asarray([obj.pixel_thresh_vol[0]*scale**2*zscale for obj in filt_cc])
yv = np.asarray([(obj.pixel_thresh_vol[-1]+obj.segment_vol_bud[-1]-obj.pixel_thresh_vol[0])*scale**2*zscale/((len(obj.frames)-1)*timestep) for obj in filt_cc])
temp = scipy.stats.pearsonr(xv, yv/xv)
plt.plot(xv, yv/xv, marker='.', label=' Volume. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(filt_cc)), linestyle='None')
# plt.hexbin(np.asarray(xv)*scale**3, yv, label=' Volume. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 1), len(filt_cc)))
plt.xlabel('$V_b$ ($\mu m^3$)')
plt.ylabel('Volume rate of production/$V_b$')
plt.legend()
plt.title('Volume rate of change')
fig.savefig(directory+'/dv_specific_vb_point_thresh.png', bbox_inches='tight', dpi=fig.dpi)
del fig


fig = plt.figure(figsize=[7, 7])
xv = [obj.vb for obj in filt_cc]
yv = [(obj.vd+obj.vbud[-1]-obj.vb)/((len(obj.frames)-1)*timestep) for obj in filt_cc]
temp = scipy.stats.pearsonr(np.asarray(xv)*scale**3, yv)
plt.plot(np.asarray(xv)*scale**3, yv, marker='.', label=' Volume. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(filt_cc)), linestyle='None')
# plt.hexbin(np.asarray(xv)*scale**3, yv, label=' Volume. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 1), len(filt_cc)))
plt.xlabel('$V_b$ ($\mu m^3$)')
plt.ylabel('Volume rate of production')
plt.legend()
plt.title('Volume rate of change')
fig.savefig(directory+'/dv_vb_point.png', bbox_inches='tight', dpi=fig.dpi)
del fig


fig = plt.figure(figsize=[7, 7])
xv = np.asarray([obj.vb for obj in filt_cc])*scale**3
yv = np.asarray([(obj.vd+obj.vbud[-1]-obj.vb)/((len(obj.frames)-1)*timestep) for obj in filt_cc])*scale**3
temp = scipy.stats.pearsonr(xv, yv/xv)
plt.plot(xv, yv/xv, marker='.', label=' Volume. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(filt_cc)), linestyle='None')
# plt.hexbin(np.asarray(xv)*scale**3, yv, label=' Volume. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 1), len(filt_cc)))
plt.xlabel('$V_b$ ($\mu m^3$)')
plt.ylabel('Volume rate of production')
plt.legend()
plt.title('Volume rate of change')
fig.savefig(directory+'/dv_specific_vb_point.png', bbox_inches='tight', dpi=fig.dpi)
del fig


fig = plt.figure(figsize=[7, 7])
xv = np.asarray([obj.vb for obj in filt_cc])*scale**3
yv = np.asarray([obj.pixel_thresh_vol[0]*scale**2*zscale for obj in filt_cc])
temp = scipy.stats.pearsonr(xv, yv)
plt.plot(xv, yv, marker='.', label=' Volume. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(filt_cc)), linestyle='None')
plt.xlabel('$V_b$ ($\mu m^3$)')
plt.ylabel('$V_b$ Threshold estimate ($\mu m^3$)')
plt.legend()
plt.title('Volume measurement')
fig.savefig(directory+'/vb_two_meas_corr.png', bbox_inches='tight', dpi=fig.dpi)
del fig


fig = plt.figure(figsize=[7, 7])
xv = np.asarray([obj.zproj_fl_c2[0] for obj in filt_cc])
yv = np.asarray([obj.pixel_thresh_fluor_vals_c2[0] for obj in filt_cc])
temp = scipy.stats.pearsonr(xv, yv)
plt.plot(xv, yv, marker='.', label='PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(filt_cc)), linestyle='None')
plt.xlabel('$F_b$ z projected')
plt.ylabel('$F_b$ Fluorescence Threshold estimate')
plt.legend()
plt.title('Fluorescence measurement')
fig.savefig(directory+'/fb_two_meas_corr.png', bbox_inches='tight', dpi=fig.dpi)
del fig


if fluor_c2:
    fig = plt.figure(figsize=[7, 7])
    xv = [obj.vb for obj in filt_cc]
    yv = [(obj.zproj_fl_c2[-1] + obj.zproj_fl_bud_c2[-1] - obj.zproj_fl_c2[0]) / ((len(obj.frames) - 1) * timestep) for obj in
          filt_cc]
    temp = scipy.stats.pearsonr(np.asarray(xv) * scale ** 3, yv)
    plt.plot(np.asarray(xv) * scale ** 3, yv, marker='.',
             label=' Z sum Fluor. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3), np.round(temp[1], 2),
                                                                           len(filt_cc)), linestyle='None')
    # plt.hexbin(np.asarray(xv) * scale ** 3, yv,
    #            label=' Z sum Fluor. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3), np.round(temp[1], 1),
    #                                                                          len(filt_cc)))
    plt.xlabel('$V_b$ ($\mu m^3$)')
    plt.ylabel('Fluorescence rate of production')
    plt.legend()
    plt.title('Fluorescence rate of change')
    fig.savefig(directory + '/df_zproj_c2_vb_point.png', bbox_inches='tight', dpi=fig.dpi)
    del fig


fig = plt.figure(figsize=[7, 7])
celltype = ['Mothers', 'Daughters']
import seaborn as sns
for i0 in range(2):
    xv = [len(obj.frames) for obj in filt_cc if obj.celltype==i0]
    sns.distplot(np.asarray(xv)*timestep, label=celltype[i0]+', num cells={0}'.format(len(xv)))
xv = [len(obj.frames)*timestep for obj in filt_cc]
sns.distplot(np.asarray(xv), label='Population, num cells={0}'.format(len(xv)))

plt.xlabel('$t_{div}$ (minutes)')
plt.legend()
plt.title('Interdivision time by celltype')
fig.savefig(directory+'/tdiv_dists.png', bbox_inches='tight', dpi=fig.dpi)
del fig
