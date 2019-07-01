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

# USE FITTING1, NOT THIS VERSION

# yFB79 on 181212 with csm or cycloheximide (mCherry)
base_path, expt_path = '/scratch/lab/image_analysis_scratch', ['/181212_yFB79_CSM_dex_fluor_maturation/CSM_timelapse',
                                                               '/181212_yFB79_CSM_dex_fluor_maturation/CHX_timelapse']
num_frames = [30, 40]  # should be number of frames analyzed
num_scenes = [5, 6]  # should be number of scenes analyzed
timestep = 4.0
timepoint_added = 4

# # yFB7 on 181213 with csm or cycloheximide (mVenNB)
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', ['/181213_yFB7_fluor_maturation/CSM_timelapse',
#                                                                '/181213_yFB7_fluor_maturation/CHX_timelapse']
# num_frames = [46, 46]  # should be number of frames analyzed
# num_scenes = [5, 6]  # should be number of scenes analyzed
# timestep = 2.0
# timepoint_added = 5
time_fitting = 100/np.int(timestep)  # number of minutes in up to which the data should be fitted
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
        # integrating the fluorescence at each timepoint
        for i1 in range(1, num_frames[i0]+1):
            temp_fl_scalar=0.0
            for obj in temp_cells:
                if i1 in obj.frames:
                    temp_fl_scalar += obj.zproj_fluor_vals[obj.frames.index(i1)]
            int_fluor[scene-1, i1-1] = temp_fl_scalar
    temp_norm = int_fluor[:, timepoint_added-1]
    int_fluor = int_fluor/np.repeat(temp_norm[:, np.newaxis], num_frames[i0], axis=1)
    # fitting maturation curve
    if i0==1:  # only do this for the cycloheximide version
        # time_fitting = np.argmax(np.mean(int_fluor,axis=0))
        time_fitting1 = time_fitting-2
        # print time_fitting
        mat_vals = np.mean(int_fluor,axis=0)[timepoint_added:time_fitting1]  # finish the fitting 2 timepoints before the max value
        mat_vals_2 = np.log((np.mean(int_fluor,axis=0)[time_fitting]-mat_vals)/(np.amax(np.mean(int_fluor,axis=0))-mat_vals[0]))
        # mat_vals_2 = np.log((1.6 - mat_vals) / (1.6 - mat_vals[0]))
        x1 = timestep*np.asarray(range(time_fitting1-timepoint_added))
        x2=x1[:,np.newaxis]
        # print x2.shape, mat_vals_2.shape, mat_vals_2
        # plt.plot(x2, mat_vals_2, label='Fitting')
        # plt.show()
        temp_vals = scipy.stats.linregress(x1, mat_vals_2)
        # temp_vals = np.linalg.lstsq(x2,-np.log(mat_vals_2)) # this should give the slope
        # print temp_vals

        # exit()
        x3 = timestep*np.asarray(range(timepoint_added,time_fitting))
        y1 = 1+(np.mean(int_fluor,axis=0)[time_fitting]-mat_vals[0])*(1-np.exp(temp_vals[0]*(x3-timestep*timepoint_added)))
        plt.plot(x3, y1, label=r'fitting, $\tau_{1/2}=$'+str(np.around(-np.log(2)/temp_vals[0],1))+' minutes',color='b')
    plt.plot(x, np.mean(int_fluor, axis=0), label=labels[i0])
    plt.fill_between(x, np.mean(int_fluor, axis=0)-np.std(int_fluor, axis=0), np.mean(int_fluor, axis=0)+
                     np.std(int_fluor, axis=0), alpha=0.5)
plt.axvline(x=timestep*(timepoint_added-1), label='Cycloheximide added',c='k')  # note that zero is timepoint # 1
plt.xlabel('Time (minutes)')
plt.ylabel('Relative Integrated Fluorescence')
plt.legend()
fig.savefig(directory+'/Integrated_fluorescence_fitting.png', bbox_inches='tight', dpi=fig.dpi)
print
del fig
