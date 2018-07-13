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

# This script takes the cell output from populate_cells_2.py and generates fluorescence plots.
scene=1
timestep = 10
pixel_size={'60X': 0.267, '100X':0.16}
with open('./scene_{0}/cells_fl_scene_{0}.pkl'.format(scene), 'rb') as input:
    c = pickle.load(input)

fig = plt.figure(figsize=[6,6])
gr = []
# print c[10].int_fl
# for obj in [obj for obj in c if ((obj.area[0]<50) & (len(obj.frames)>5))]:
for obj in [obj for obj in c if (len(obj.frames) > 5)]:
    x = timestep*(np.asarray(obj.frames)-obj.frames[0])
    y = np.asarray(obj.int_fl)
    plt.plot(x, y, label='cell {0}'.format(obj.index), alpha=0.5)
    gr.append(scipy.stats.linregress(x, y)[0])
    # print scipy.stats.linregress(x, y)[0]
# plt.legend()
plt.xlabel('Time (minutes)')
plt.ylabel('Fluorescence Signal')
fig.savefig('./scene_1/fluor_traces_all_cells.png',bbox_inches='tight')

fig = plt.figure(figsize=[6,6])
gr = []
# print c[10].int_fl
for obj in [obj for obj in c if ((obj.area[0]<50) & (len(obj.frames)>5))]:
# for obj in [obj for obj in c if (len(obj.frames) > 5)]:
    x = timestep*(np.asarray(obj.frames)-obj.frames[0])
    y = np.asarray(obj.int_fl)
    plt.plot(x, y, label='cell {0}'.format(obj.index), alpha=0.5)
    gr.append(scipy.stats.linregress(x, y)[0])
    # print scipy.stats.linregress(x, y)[0]
plt.legend()
plt.xlabel('Time (minutes)')
plt.ylabel('Fluorescence Signal')
fig.savefig('./scene_1/fluor_traces_small_cells_only.png',bbox_inches='tight')