import numpy as np
import scipy
from scipy import stats
import scipy.optimize
import scipy.io as sio
import weakref
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from scipy import stats
import pandas as pd
import cPickle as pickle
import custom_image_toolkit as C

# takes the output of populate_cells_all_scenes_3.py
base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt'
timestep = 10
pixel_size={'60X': 0.267, '100X':0.16}
# labeled is
labeling = {'labeled':'yFB29', 'unlabeled':'yFB46'}
# loading our data
c=[]
for scene in range(1, 4):
    with open(base_path+expt_path+'/timelapse/scene_{0}/outputs/cells_fl_lab_scene_{0}.pkl'.format(scene), 'rb') as input:
        temp=pickle.load(input)
    for obj in temp:
        c.append(obj)
    print len(c)


# plotting volume traces for all cells
fig=plt.figure(figsize=[6,6])
a=[temp for temp in c if temp.type==1]
for obj in a:
    plt.plot(timestep*(np.asarray(obj.frames)-obj.frames[0]), np.asarray(obj.ellipse_volume)*pixel_size['60X']**3, color='r', alpha=0.5)
b=[temp for temp in c if temp.type==0]
for obj in b:
    plt.plot(timestep*(np.asarray(obj.frames)-obj.frames[0]), np.asarray(obj.ellipse_volume)*pixel_size['60X']**3, color='b', alpha=0.5)
plt.xlabel('Time (minutes)')
plt.ylabel('Volume ($\mu m^3$)')
plt.title('Red=Labeled {0} ({1}), Blue=Unlabeled {2} ({3})'.format(labeling['labeled'], len(a), labeling['unlabeled'], len(b)))
fig.savefig(base_path+expt_path+'/figure_volume_traces.png',bbox_inches='tight')
del fig

# plotting fluorescence traces for all cells
fig=plt.figure(figsize=[6,6])
a=[temp for temp in c if temp.type==1]
for obj in a:
    plt.plot(timestep*(np.asarray(obj.frames)-obj.frames[0]), np.asarray(obj.int_fl), color='r', alpha=0.5)
b=[temp for temp in c if temp.type==0]
for obj in b:
    plt.plot(timestep*(np.asarray(obj.frames)-obj.frames[0]), np.asarray(obj.int_fl), color='b', alpha=0.5)
plt.xlabel('Time (minutes)')
plt.ylabel('Integrated Fluorescence (A.U.)')
plt.title('Red=Labeled {0} ({1}), Blue=Unlabeled {2} ({3})'.format(labeling['labeled'], len(a), labeling['unlabeled'], len(b)))
fig.savefig(base_path+expt_path+'/figure_fluor_traces.png',bbox_inches='tight')
del fig

# plotting volume traces for small cells only
fig=plt.figure(figsize=[6,6])
a=[temp for temp in c if temp.type==1 and temp.area[0]<100]
for obj in a:
    plt.plot(timestep*(np.asarray(obj.frames)-obj.frames[0]), np.asarray(obj.ellipse_volume)*pixel_size['60X']**3, color='r', alpha=0.5)
b=[temp for temp in c if temp.type==0 and temp.area[0]<100]
for obj in b:
    plt.plot(timestep*(np.asarray(obj.frames)-obj.frames[0]), np.asarray(obj.ellipse_volume)*pixel_size['60X']**3, color='b', alpha=0.5)
plt.xlabel('Time (minutes)')
plt.ylabel('Volume ($\mu m^3$)')
plt.title('Red=Labeled {0} ({1}), Blue=Unlabeled {2} ({3})'.format(labeling['labeled'], len(a), labeling['unlabeled'], len(b)))
fig.savefig(base_path+expt_path+'/figure_volume_traces_small_cells.png',bbox_inches='tight')
del fig

# plotting fluorescence traces for small cells only
fig=plt.figure(figsize=[6,6])
a=[temp for temp in c if temp.type==1 and temp.area[0]<100]
for obj in a:
    plt.plot(timestep*(np.asarray(obj.frames)-obj.frames[0]), np.asarray(obj.int_fl), color='r', alpha=0.5)
b=[temp for temp in c if temp.type==0 and temp.area[0]<100]
for obj in b:
    plt.plot(timestep*(np.asarray(obj.frames)-obj.frames[0]), np.asarray(obj.int_fl), color='b', alpha=0.5)
plt.xlabel('Time (minutes)')
plt.ylabel('Integrated Fluorescence (A.U.)')
plt.title('Red=Labeled {0} ({1}), Blue=Unlabeled {2} ({3})'.format(labeling['labeled'], len(a), labeling['unlabeled'], len(b)))
fig.savefig(base_path+expt_path+'/figure_fluor_traces_small_cells.png',bbox_inches='tight')
del fig