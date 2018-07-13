import numpy as np
import scipy
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

base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt'
color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))

scene=1
c = []  # this will track all the cells over this timelapse
for frame_num in range(1,31):
    filename = '180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w1Brightfield confocal_s{0}_t{1}_segmentation.mat'.format(str(scene),str(frame_num).zfill(2))
    path = base_path+expt_path+'/timelapse/segments/'+filename
    tracking_csv_file = pd.DataFrame.from_csv(base_path+expt_path+'/timelapse/segments/tracking/tracking.csv', index_col=None)
    c, im = C.single_frame(path, c, frame_num, tracking_csv_file, color_seq_val)
    # plt.show(im)
    im.savefig('./scene_{0}/images/frame_{1}.png'.format(scene,frame_num))
    del im
C.save_object(c, './scene_{0}/cells_scene_{1}.pkl'.format(scene, scene))
# del c
# with open('./scene_1/cells_scene_1.pkl', 'rb') as input:
#     c = pickle.load(input)
# for obj in c[:50]:
#     print obj.index, obj.frames


# Calculating volume: orientation, major axis length, minor axis length