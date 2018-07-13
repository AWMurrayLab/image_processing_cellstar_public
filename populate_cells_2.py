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

# this takes the cell output from script populate_cells_1.py and adds fluorescence data to it
base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt'
color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
pixel_size = {'60X': 0.267, '100X': 0.16}
scene = 1
z_scale, z_offset = 0.4/pixel_size['60X'], 0

with open('./scene_1/cells_scene_1.pkl', 'rb') as input:
    c = pickle.load(input)
frame_list = [obj.frames for obj in c]
update_list = []
for obj in c:
    obj.add_fluor_placeholders()
for frame_num in range(1, 31):
    temp = [(frame_num in temp1) for temp1 in frame_list]
    update_list.append([i for i, e in enumerate(temp) if e != 0])  # gives the list of indices that have to be addressed
    # at each frame
    # print update_list
    filename_fl = '/180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w2515 laser 20_s{0}_t{1}.TIF'.format(str(scene), str(frame_num))
    filename_fl_bkgd = '/180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w2515 laser 20_s10_t{0}.TIF'.format(str(frame_num))
    # format is z, y, x
    bkgd_im = io.imread(base_path+expt_path+'/timelapse'+filename_fl_bkgd)
    bkgd_im1 = np.zeros(bkgd_im.shape)
    temp1 = np.mean(bkgd_im, axis=0)  # taking the mean with respect to the z axis. Do it this way since there doesn't seem
    # to be any systematic bias in that direction.
    for i0 in range(bkgd_im.shape[0]):
        bkgd_im1[i0, :, :] = temp1[:, :]
    del temp1, bkgd_im
    c, mask = C.add_fluorescence_traces_v1(temp_path=base_path+expt_path+'/timelapse'+filename_fl, temp_cells=c, frame_list=update_list[-1],
                                  current_frame=frame_num, z_scaling=z_scale, z_offset=z_offset, bkgd=bkgd_im1)
    print 'done frame {0}'.format(frame_num)
C.save_object(c, './scene_{0}/cells_fl_scene_{1}.pkl'.format(scene, scene))
#     c = C.single_frame(path, c, frame_num, tracking_csv_file, color_seq_val)
#     # plt.show(im)
#     im.savefig('./scene_{0}/frame_{1}.png'.format(scene,frame_num))
#     del im
# C.save_object(c, './scene_{0}/cells_scene_{1}.pkl'.format(scene, scene))

# del c
# with open('./scene_1/cells_scene_1.pkl', 'rb') as input:
#     c = pickle.load(input)
# for obj in c[:50]:
#     print obj.index, obj.frames
# Calculating volume: orientation, major axis length, minor axis length
