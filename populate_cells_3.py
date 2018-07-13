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
from skimage import io


# This script takes the cell output from populate_cells_2.py and adds data showing whether cells are stained or not.
scene = 1
timestep = 10
frame_num=1
pixel_size = {'60X': 0.267, '100X':0.16}
with open('./scene_{0}/cells_fl_scene_{0}.pkl'.format(scene), 'rb') as input:
    c = pickle.load(input)
base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
scenes = [1,2,3]
filename = '180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w1Brightfield confocal_s{0}_t{1}_segmentation.mat'.format(str(scene),str(frame_num).zfill(2))



