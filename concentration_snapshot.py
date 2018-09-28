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

# importing the cells
#
# # yFB29 experiment on 180831
# scale = pixel_size['60X']
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180831_yFB29_800uMGal/timelapse'
# timestep = 10.0
# fluor_c2 = False
# num_scenes = 7  # num_scenes should be the number of scenes to analyze + 1
# num_frames = np.asarray([56, 56, 56, 56, 56, 56])
# temp_frame_num = 55

# pACT1-mKate2 experiment on 180910
scale = pixel_size['60X']
base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180910_pACT1_mKate2/timelapse'
# path for the desktop storage of data.
# base_path, expt_path = '/mnt/d/Lab/image_analysis', '/180831_yFB29_800uMGal'  # path for the laptop storage of images
# script_path = '/mnt/c/Users/felix/Documents/GitHub/image_processing_cellstar'
timestep = 5.0
fluor_c2 = True
num_scenes = 8  # num_scenes should be the number of scenes to analyze + 1
num_frames = np.asarray([55, 70, 65, 66, 51, 66, 66])
temp_frame_num = 50

# beginning analysis
relevant_cells = []
for scene in range(1, num_scenes):
    print 'Scene Number: {0}'.format(scene)
    data_index = [base_path + expt_path, scene]
    directory = base_path + expt_path + '/scene_{0}/outputs'.format(scene)
    with open(directory + '/cells_scene_{0}_v3.pkl'.format(scene), 'rb') as input:
        c = pickle.load(input)
    relevant_cells += [obj for obj in c if temp_frame_num in obj.frames]
f_c1 = np.asarray([obj.zproj_fluor_vals[obj.frames.index(temp_frame_num)] for obj in relevant_cells])
v = np.asarray([obj.ellipse_volume[obj.frames.index(temp_frame_num)] for obj in relevant_cells])

# generating figures

directory = base_path + expt_path + '/plots'
if not os.path.exists(directory):
    os.makedirs(directory)

fig = plt.figure(figsize=[7, 7])
xv = v
yv = f_c1/v
temp = scipy.stats.pearsonr(np.asarray(xv)*scale**3, yv)
plt.plot(np.asarray(xv)*scale**3, yv, marker='.', label=' F1 conc. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(xv)), linestyle='None', alpha=0.2)
# plt.hexbin(np.asarray(xv)*scale**3, yv)
plt.xlabel('$V_b$ ($\mu m^3$)')
plt.ylabel('Channel 1 fluor concentration')
plt.legend()
plt.title('Snapshot fluorescence concentration')
fig.savefig(directory+'/snap_f1conc_vol.png', bbox_inches='tight', dpi=fig.dpi)
del fig

fig = plt.figure(figsize=[7, 7])
xv = v
yv = f_c1
temp = scipy.stats.pearsonr(np.asarray(xv)*scale**3, yv)
plt.plot(np.asarray(xv)*scale**3, yv, marker='.', label=' F1 int. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(xv)), linestyle='None', alpha=0.2)
# plt.hexbin(np.asarray(xv)*scale**3, yv)
plt.xlabel('$V_b$ ($\mu m^3$)')
plt.ylabel('Channel 1 fluor abundance')
plt.legend()
plt.title('Snapshot fluorescence abundance')
fig.savefig(directory+'/snap_f1ab_vol.png', bbox_inches='tight', dpi=fig.dpi)
del fig


# if there is a second fluorescence channel.
if fluor_c2:
    f_c2 = np.asarray([obj.zproj_fluor_vals_c2[obj.frames.index(temp_frame_num)] for obj in relevant_cells])
    fig = plt.figure(figsize=[7, 7])
    xv = v
    yv = f_c2 / v
    temp = scipy.stats.pearsonr(np.asarray(xv) * scale ** 3, yv)
    plt.plot(np.asarray(xv) * scale ** 3, yv, marker='.',
             label=' F2 conc. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3), np.round(temp[1], 2),
                                                                       len(xv)), linestyle='None', alpha=0.2)
    # plt.hexbin(np.asarray(xv)*scale**3, yv)
    plt.xlabel('$V_b$ ($\mu m^3$)')
    plt.ylabel('Channel 2 fluor concentration')
    plt.legend()
    plt.title('Snapshot fluorescence concentration')
    fig.savefig(directory + '/snap_f2conc_vol.png', bbox_inches='tight', dpi=fig.dpi)
    del fig

    fig = plt.figure(figsize=[7, 7])
    xv = v
    yv = f_c2
    temp = scipy.stats.pearsonr(np.asarray(xv)*scale**3, yv)
    plt.plot(np.asarray(xv)*scale**3, yv, marker='.', label=' F2 int. PCC={0}, pval={1}, num cells={2}'.format(np.round(temp[0], 3),np.round(temp[1], 2), len(xv)), linestyle='None', alpha=0.2)
    # plt.hexbin(np.asarray(xv)*scale**3, yv)
    plt.xlabel('$V_b$ ($\mu m^3$)')
    plt.ylabel('Channel 2 fluor abundance')
    plt.legend()
    plt.title('Snapshot fluorescence abundance')
    fig.savefig(directory+'/snap_f2ab_vol.png', bbox_inches='tight', dpi=fig.dpi)
    del fig