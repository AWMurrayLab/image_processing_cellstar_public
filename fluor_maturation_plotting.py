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
pixel_size = {'60X': 0.267, '100X': 0.16}

# # yFB7 on 180731 with csm or cycloheximide
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', ['/180731_csm_fluor_mat/timelapse',
#                                                                '/180731_cycloheximide_fluor_mat/timelapse']
# num_frames = [40, 40]
# num_scenes = 6
# timestep = 2.0
# timepoint_added=4

# # yFB79 on 181204 with csm or cycloheximide
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', ['/181204_yFB79_fluor_maturation/CSM_timelapse',
#                                                                '/181204_yFB79_fluor_maturation/CHX_timelapse']
# num_frames = [46, 46]  # should be number of frames analyzed
# num_scenes = [5, 4]  # should be number of scenes analyzed
# timestep = 2.0
# timepoint_added = 5

# # yFB79 on 181212 with csm or cycloheximide (mCherry)
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', ['/181212_yFB79_CSM_dex_fluor_maturation/CSM_timelapse',
#                                                                '/181212_yFB79_CSM_dex_fluor_maturation/CHX_timelapse']
# num_frames = [30, 40]  # should be number of frames analyzed
# num_scenes = [5, 6]  # should be number of scenes analyzed
# timestep = 4.0
# timepoint_added = 4
# fluor = 'mCherry'

# yFB7 on 181213 with csm or cycloheximide (mVenNB)
base_path, expt_path = '/scratch/lab/image_analysis_scratch', ['/181213_yFB7_fluor_maturation/CSM_timelapse',
                                                               '/181213_yFB7_fluor_maturation/CHX_timelapse']
num_frames = [46, 46]  # should be number of frames analyzed
num_scenes = [5, 6]  # should be number of scenes analyzed
timestep = 2.0
timepoint_added = 5
fluor = 'mVenNB'

c = []
labels=['CSM', 'CSM + $125\mu M$ cyc.']
for i0 in range(len(expt_path)):
    with open(base_path + expt_path[i0] + '/cells_compiled.pkl','rb') as input:
        c.append(pickle.load(input))

directory = base_path + expt_path[0] + '/plots'
if not os.path.exists(directory):
    os.makedirs(directory)
fig=plt.figure(figsize=[5,5])


for i0 in range(len(c)):
    int_fluor = np.zeros([num_scenes[i0], num_frames[i0]])
    x = timestep * np.asarray(range(num_frames[i0]))
    for scene in range(1, num_scenes[i0]+1):
        temp_cells = [obj for obj in c[i0] if obj.frames[-1]==num_frames[i0] and not(obj.edge_cell) and obj.scene_num==scene]
        print i0, scene, len(temp_cells)
        # integrating the fluorescence at each timepoint
        for i1 in range(1, num_frames[i0]+1):
            temp_fl_scalar=0.0
            for obj in temp_cells:
                if i1 in obj.frames:
                    temp_fl_scalar += obj.zproj_fluor_vals[obj.frames.index(i1)]
            int_fluor[scene-1, i1-1] = temp_fl_scalar
    temp_norm = int_fluor[:, timepoint_added-1]
    int_fluor = int_fluor/np.repeat(temp_norm[:, np.newaxis], num_frames[i0], axis=1)
    plt.plot(x, np.mean(int_fluor, axis=0), label=labels[i0])
    plt.fill_between(x, np.mean(int_fluor, axis=0)-np.std(int_fluor, axis=0), np.mean(int_fluor, axis=0)+
                     np.std(int_fluor, axis=0), alpha=0.5)
plt.axvline(x=timestep*(timepoint_added-1), label='Cycloheximide added',c='k')  # note that zero is timepoint # 1
plt.xlabel('Time (minutes)')
plt.ylabel('Relative Integrated Fluorescence')
plt.legend()
fig.savefig(directory+'/Integrated_fluorescence_full_trace.png', bbox_inches='tight', dpi=fig.dpi)
del fig

fig=plt.figure(figsize=[5,5])
for i0 in range(len(c)):
    int_fluor = np.zeros([num_scenes[i0], num_frames[i0]])
    x = timestep * np.asarray(range(num_frames[i0]))
    for scene in range(1, num_scenes[i0]+1):
        temp_cells = [obj for obj in c[i0] if obj.frames[-1]==num_frames[i0] and not(obj.edge_cell) and obj.scene_num==scene]
        print i0, scene, len(temp_cells)
        # integrating the fluorescence at each timepoint
        for i1 in range(1, num_frames[i0]+1):
            temp_fl_scalar=0.0
            for obj in temp_cells:
                if i1 in obj.frames:
                    temp_fl_scalar += obj.zproj_fluor_vals[obj.frames.index(i1)]
            int_fluor[scene-1, i1-1] = temp_fl_scalar
    temp_norm = int_fluor[:, timepoint_added-1]
    int_fluor = int_fluor/np.repeat(temp_norm[:, np.newaxis], num_frames[i0], axis=1)
    plt.plot(x, np.mean(int_fluor, axis=0), label=labels[i0])
    plt.fill_between(x, np.mean(int_fluor, axis=0)-np.std(int_fluor, axis=0), np.mean(int_fluor, axis=0)+
                     np.std(int_fluor, axis=0), alpha=0.5)
plt.axvline(x=timestep*(timepoint_added-1), label='Cycloheximide added',c='k')  # note that zero is timepoint # 1
plt.xlabel('Time (minutes)')
plt.ylabel('Relative Integrated Fluorescence')
plt.legend()
fig.savefig(directory+'/Integrated_fluorescence_partial_trace.png', bbox_inches='tight', dpi=fig.dpi)
del fig

fig=plt.figure(figsize=[5,5])
for i0 in range(len(c)):
    int_fluor = np.zeros([num_scenes[i0], num_frames[i0]])
    x = timestep * np.asarray(range(num_frames[i0]))
    for scene in range(1, num_scenes[i0]+1):
        temp_cells = [obj for obj in c[i0] if obj.frames[-1]==num_frames[i0] and not(obj.edge_cell) and obj.scene_num==scene]
        print i0, scene, len(temp_cells)
        # integrating the fluorescence at each timepoint
        for i1 in range(1, num_frames[i0]+1):
            temp_fl_scalar=0.0
            for obj in temp_cells:
                if i1 in obj.frames:
                    temp_fl_scalar += obj.ellipse_volume[obj.frames.index(i1)]
            int_fluor[scene-1, i1-1] = temp_fl_scalar
    temp_norm = int_fluor[:, timepoint_added-1]
    int_fluor = int_fluor/np.repeat(temp_norm[:, np.newaxis], num_frames[i0], axis=1)
    plt.plot(x, np.mean(int_fluor, axis=0), label=labels[i0])
    plt.fill_between(x, np.mean(int_fluor, axis=0)-np.std(int_fluor, axis=0), np.mean(int_fluor, axis=0)+
                     np.std(int_fluor, axis=0), alpha=0.5)
plt.axvline(x=timestep*(timepoint_added-1), label='Cycloheximide added',c='k')  # note that zero is timepoint # 1
plt.xlabel('Time (minutes)')
plt.ylabel('Relative Integrated Volume')
plt.title(fluor+' Relative integrated volume curve')
plt.legend()
fig.savefig(directory+'/Integrated_volume_full_trace.png', bbox_inches='tight', dpi=fig.dpi)
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/microscopy_validation/'+fluor+'_volume_curve.png', bbox_inches='tight', dpi=fig.dpi)

fig=plt.figure(figsize=[5,5])
for i0 in range(len(c)):
    temp_cells = [obj for obj in c[i0] if obj.frames[-1]==num_frames[i0] and not(obj.edge_cell)]
    int_fluor = np.zeros(num_frames[i0])
    x = timestep*np.asarray(range(num_frames[i0]))
    # integrating the fluorescence at each timepoint
    for i1 in range(1, num_frames[i0]+1):
        temp_fl_scalar=0.0
        for obj in temp_cells:
            if i1 in obj.frames:
                temp_fl_scalar += obj.zproj_fluor_vals[obj.frames.index(i1)]
        int_fluor[i1-1]=temp_fl_scalar
    plt.plot(x[:10], int_fluor[:10]/int_fluor[0], label=labels[i0])
plt.xlabel('Time (minutes)')
plt.ylabel('Relative Integrated Fluorescence')
plt.legend()
fig.savefig(directory+'/Integrated_fluorescence.png', bbox_inches='tight', dpi=fig.dpi)