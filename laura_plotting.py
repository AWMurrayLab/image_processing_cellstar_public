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


# Expt with Laura 181015
pixel_size = {'60X': 0.267, '100X': 0.16}
image_filename, bf_filename, fl_filename = '/2018_10_15_yLB256_yLB365_yFJB71_cellasics_01_', \
                                           'w1Brightfield confocal_', \
                                           'w2594 laser 20_'  # note to change fluorescence path to match laser power
base_path, expt_path = '/scratch/lab/image_analysis_scratch', ['/181015_spinning_disk/timelapse']
scenes = [3, 5, 7, 9, 10, 11, 12, 13]
num_scenes = len(scenes)
bkgd_scene = 22  # the last scene is the background
num_frames = 22*np.ones(bkgd_scene, dtype=int)
drange = 65535.0  # image fluorescence maximum
bf_base_name = '/2018_10_15_yLB256_yLB365_yFJB71_cellasics_01_w1Brightfield confocal'
date = '/2018_10_15'
bkgd_details = None
timestep = 15.0
scale = pixel_size['60X']
frame_max = 21

c = []
for i0 in range(len(expt_path)):
    with open(base_path + expt_path[i0] + '/cells_compiled.pkl','rb') as input:
        c.append(pickle.load(input))

directory = base_path + expt_path[0] + '/plots'
if not os.path.exists(directory):
    os.makedirs(directory)
print len(c)
fig=plt.figure(figsize=[5,5])
x = np.asarray([obj.ellipse_volume[-1]*scale**3 for obj in c[0] if np.sum(obj.nuclear_whi5)==0 and frame_max in obj.frames])
y = np.asarray([obj.zproj_fluor_vals[-1] for obj in c[0] if np.sum(obj.nuclear_whi5)==0 and frame_max in obj.frames])
vals = scipy.stats.pearsonr(x, y)
plt.plot(x, y, '.', label='PCC={0}, pval={1}, num cells={2}'.format(np.around(vals[0], 3), np.around(vals[1], 3), len(x)), alpha=0.5)
plt.xlabel('Volume')
plt.ylabel('Integrated Fluorescence')
plt.legend()
plt.title('Fluorescence Abundance vs. volume')
fig.savefig(directory+'/Integrated_fluorescence.png', bbox_inches='tight', dpi=fig.dpi)
del fig

fig=plt.figure(figsize=[5,5])
x = np.asarray([obj.ellipse_volume[-1]*scale**3 for obj in c[0] if np.sum(obj.nuclear_whi5)==0 and frame_max in obj.frames])
y = np.asarray([obj.zproj_fluor_vals[-1] for obj in c[0] if np.sum(obj.nuclear_whi5)==0 and frame_max in obj.frames])/x
vals = scipy.stats.pearsonr(x, y)
plt.plot(x, y, '.', label='PCC={0}, pval={1}, num cells={2}'.format(np.around(vals[0], 3), np.around(vals[1], 3), len(x)), alpha=0.5)
plt.xlabel('Volume')
plt.ylabel('Fluorescence concentration')
plt.legend()
plt.title('Fluorescence concentration vs. volume')
fig.savefig(directory+'/concentrated_fluorescence.png', bbox_inches='tight', dpi=fig.dpi)
del fig

fig=plt.figure(figsize=[5,5])
x = np.asarray([obj.ellipse_volume[0]*scale**3 for obj in c[0] if np.sum(obj.nuclear_whi5)==0 and len(obj.frames)>3])
yv=[]
for obj in c[0]:
    if np.sum(obj.nuclear_whi5) == 0 and len(obj.frames)>3:
        tempy = obj.zproj_fluor_vals[:4]
        tempx = obj.frames[:4]
        tempv = scipy.stats.linregress(np.asarray(tempx)*timestep, tempy)
        yv.append(tempv[0])
y = np.asarray(yv)/x
vals = scipy.stats.pearsonr(x, y)
plt.plot(x, y, '.', label='PCC={0}, pval={1}, num cells={2}'.format(np.around(vals[0], 3), np.around(vals[1], 3), len(x)), alpha=0.5)
plt.xlabel('Volume')
plt.ylabel('Specific rate of fluor production (normalized to volume)')
plt.legend()
plt.title('Specific rate of fluorescence production vs. volume')
fig.savefig(directory+'/fluor_prodn_rate.png', bbox_inches='tight', dpi=fig.dpi)
del fig


fig=plt.figure(figsize=[5,5])
x = np.asarray([obj.ellipse_volume[0]*scale**3 for obj in c[0] if np.sum(obj.nuclear_whi5)==0 and len(obj.frames)>5])
yv=[]
for obj in c[0]:
    if np.sum(obj.nuclear_whi5) == 0 and len(obj.frames)>5:
        tempy = obj.zproj_fluor_vals[3:6]
        tempx = obj.frames[3:6]
        tempv = scipy.stats.linregress(np.asarray(tempx)*timestep, tempy)
        yv.append(tempv[0])
y = np.asarray(yv)/x
vals = scipy.stats.pearsonr(x, y)
plt.plot(x, y, '.', label='PCC={0}, pval={1}, num cells={2}'.format(np.around(vals[0], 3), np.around(vals[1], 3), len(x)), alpha=0.5)
plt.xlabel('Volume')
plt.ylabel('Specific rate of fluor production (normalized to volume)')
plt.legend()
plt.title('Specific rate of fluorescence production vs. volume')
fig.savefig(directory+'/fluor_prodn_rate_excl_buds.png', bbox_inches='tight', dpi=fig.dpi)
del fig


fig=plt.figure(figsize=[5,5])
x = np.asarray([obj.ellipse_volume[0]*scale**3 for obj in c[0] if np.sum(obj.nuclear_whi5)==0 and len(obj.frames)>5])
yv=[]
for obj in c[0]:
    if np.sum(obj.nuclear_whi5) == 0 and len(obj.frames)>5:
        tempy = obj.zproj_fluor_vals[3:6]
        tempx = obj.frames[3:6]
        tempv = scipy.stats.linregress(np.asarray(tempx)*timestep, tempy)
        yv.append(tempv[0])
y = np.asarray(yv)
vals = scipy.stats.pearsonr(x, y)
plt.plot(x, y, '.', label='PCC={0}, pval={1}, num cells={2}'.format(np.around(vals[0], 3), np.around(vals[1], 3), len(x)), alpha=0.5)
plt.xlabel('Volume')
plt.ylabel('Absolute rate of fluor production (normalized to volume)')
plt.legend()
plt.title('Absolute rate of fluorescence production vs. volume')
fig.savefig(directory+'/absolute_fluor_prodn_rate_excl_buds.png', bbox_inches='tight', dpi=fig.dpi)
del fig


# for i0 in range(len(c)):
#     int_fluor = np.zeros([num_scenes, num_frames[i0]])
#     x = timestep * np.asarray(range(num_frames[i0]))
#     for scene in range(1, num_scenes+1):
#         temp_cells = [obj for obj in c[i0] if obj.frames[-1]==num_frames[i0] and not(obj.edge_cell) and obj.scene_num==scene]
#         # integrating the fluorescence at each timepoint
#         for i1 in range(1, num_frames[i0]+1):
#             temp_fl_scalar=0.0
#             for obj in temp_cells:
#                 if i1 in obj.frames:
#                     temp_fl_scalar += obj.zproj_fluor_vals[obj.frames.index(i1)]
#             int_fluor[scene-1, i1-1] = temp_fl_scalar
#     temp_norm = int_fluor[:, 0]
#     int_fluor = int_fluor/np.repeat(temp_norm[:, np.newaxis], num_frames[i0], axis=1)
#     plt.plot(x, np.mean(int_fluor, axis=0), label=labels[i0])
#     plt.fill_between(x, np.mean(int_fluor, axis=0)-np.std(int_fluor, axis=0), np.mean(int_fluor, axis=0)+
#                      np.std(int_fluor, axis=0), alpha=0.5)
# plt.axvline(x=timestep*(timepoint_added-1), label='Cycloheximide added')  # note that zero is timepoint # 1
# plt.xlabel('Time (minutes)')
# plt.ylabel('Relative Integrated Fluorescence')
# plt.legend()
# fig.savefig(directory+'/Integrated_fluorescence_full_trace.png', bbox_inches='tight', dpi=fig.dpi)
#
# fig=plt.figure(figsize=[5,5])
# for i0 in range(len(c)):
#     temp_cells = [obj for obj in c[i0] if obj.frames[-1]==num_frames[i0] and not(obj.edge_cell)]
#     int_fluor = np.zeros(num_frames[i0])
#     x = timestep*np.asarray(range(num_frames[i0]))
#     # integrating the fluorescence at each timepoint
#     for i1 in range(1, num_frames[i0]+1):
#         temp_fl_scalar=0.0
#         for obj in temp_cells:
#             if i1 in obj.frames:
#                 temp_fl_scalar += obj.zproj_fluor_vals[obj.frames.index(i1)]
#         int_fluor[i1-1]=temp_fl_scalar
#     plt.plot(x[:10], int_fluor[:10]/int_fluor[0], label=labels[i0])
# plt.xlabel('Time (minutes)')
# plt.ylabel('Relative Integrated Fluorescence')
# plt.legend()
# fig.savefig(directory+'/Integrated_fluorescence.png', bbox_inches='tight', dpi=fig.dpi)