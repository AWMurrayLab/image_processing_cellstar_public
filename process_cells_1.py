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

with open('./scene_1/cells_scene_1.pkl', 'rb') as input:
    c = pickle.load(input)

timestep = 10
pixel_size={'60X': 0.267, '100X':0.16}
# plotting volume traces
fig=plt.figure(figsize=[6,6])
for obj in c[:50]:
    plt.plot(timestep*(np.asarray(obj.frames)-obj.frames[0]), np.asarray(obj.ellipse_volume)*pixel_size['60X']**3, label='cell {0}'.format(obj.index), alpha=0.5)
plt.xlabel('Time (minutes)')
plt.ylabel('Volume ($\mu m^3$)')
fig.savefig('./scene_1/volume_traces.png',bbox_inches='tight')

fig=plt.figure(figsize=[6,6])
for obj in c[:50]:
    plt.plot(timestep*(np.asarray(obj.frames)-obj.frames[0]), np.asarray(obj.ellipse_volume)/obj.ellipse_volume[0], label='cell {0}'.format(obj.index), alpha=0.5)
plt.xlabel('Time (minutes)')
plt.ylabel('V(t)/V(0)')
fig.savefig('./scene_1/relative_volume_traces.png',bbox_inches='tight')

fig=plt.figure(figsize=[6,6])
for obj in c[:50]:
    plt.plot(timestep*(np.asarray(obj.frames)-obj.frames[0]), np.log(np.asarray(obj.ellipse_volume)/obj.ellipse_volume[0]), label='cell {0}'.format(obj.index), alpha=0.5)
plt.xlabel('Time (minutes)')
plt.ylabel('$\ln (V(t)/V(0))$')
fig.savefig('./scene_1/log_volume_traces.png',bbox_inches='tight')

fig=plt.figure(figsize=[6,6])
gr = []
for obj in [obj for obj in c if ((obj.area[0]<50) & (len(obj.frames)>5))]:
    x=timestep*(np.asarray(obj.frames)-obj.frames[0])
    y=np.log(np.asarray(obj.ellipse_volume)/obj.ellipse_volume[0])
    plt.plot(x, y, label='cell {0}'.format(obj.index), alpha=0.5)
    gr.append(scipy.stats.linregress(x, y)[0])
    # print scipy.stats.linregress(x, y)[0]
plt.legend()
plt.xlabel('Time (minutes)')
plt.ylabel('$\ln (V(t)/V(0))$')
fig.savefig('./scene_1/log_volume_traces_small_cells_only.png',bbox_inches='tight')

print 'gr', np.nanmean(gr), np.nanstd(gr), 'td', np.nanmean(np.log(2)/gr), np.nanstd(np.log(2)/gr)


fig=plt.figure(figsize=[6,6])
for obj in [obj for obj in c if ((obj.area[0]<50) & (len(obj.frames)>5))]:
    x=timestep*(np.asarray(obj.frames)-obj.frames[0])
    y=np.asarray(obj.area)*pixel_size['60X']**2
    plt.plot(x, y, label='cell {0}'.format(obj.index), alpha=0.5)
    gr.append(scipy.stats.linregress(x, y)[0])
    # print scipy.stats.linregress(x, y)[0]
plt.xlabel('Time (minutes)')
plt.ylabel('2D area ($\mu m^2$)')
fig.savefig('./scene_1/area_traces_small_cells_only.png',bbox_inches='tight')

fig=plt.figure(figsize=[10,10])
for obj in c[:100]:
    x=timestep*(np.asarray(obj.frames)-obj.frames[0])
    y=np.asarray(obj.area)*pixel_size['60X']**2
    plt.plot(x, y, label='cell {0}'.format(obj.index), alpha=0.5)
    # print obj.frames
# plt.legend()
    # gr.append(scipy.stats.linregress(x, y)[0])
    # print scipy.stats.linregress(x, y)[0]
plt.xlabel('Time (minutes)')
plt.ylabel('2D area ($\mu m^2$)')
fig.savefig('./scene_1/area_traces_all_cells.png',bbox_inches='tight')
