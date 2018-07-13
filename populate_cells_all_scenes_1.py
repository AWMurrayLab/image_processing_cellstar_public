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
import os

base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
num_frames = [31, 51, 51, 51, 51, 51, 51]
dims = 512

# scene 1 = 31 frames
# scene 2-5 = 51 frames

for scene in range(1, 8):
    c = []  # this will track all the cells over this timelapse
    # making necessary directory tree
    directory = base_path+expt_path+'/scene_{0}/outputs'.format(scene)
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists(directory+'/images'):
        os.makedirs(directory+'/images')
    outlines = np.zeros([num_frames[scene-1], dims, dims])
    for frame_num in range(1, num_frames[scene-1]):
        filename = '180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w1Brightfield confocal_s{0}_t{1}_segmentation.mat'.format(str(scene),str(frame_num).zfill(2))
        path = base_path+expt_path+'/scene_{0}/segments/'.format(scene)+filename
        tracking_csv_file = pd.DataFrame.from_csv(base_path+expt_path+'/scene_{0}/segments/tracking/tracking.csv'.format(scene), index_col=None)
        c, im, fig, temp_outlines = C.single_frame(path, c, frame_num, tracking_csv_file, color_seq_val)
        outlines[frame_num-1, :, :] = temp_outlines[:, :]

        # plt.show(im)
        # extent = im.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fig.subplots_adjust(bottom=0)
        fig.subplots_adjust(top=1)
        fig.subplots_adjust(right=1)
        fig.subplots_adjust(left=0)
        # extent = mpl.transforms.Bbox(((0, 0), (5, 5)))
        fig.savefig(directory+'/images/frame_{1}.tif'.format(scene, frame_num))
        del im, fig
    np.save(directory+'/cell_outlines_scene_{0}'.format(scene), outlines)
    C.save_object(c, directory+'/cells_scene_{0}.pkl'.format(scene))
    del c
