import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import os
import scipy
from scipy import stats

# pACT1-mKate2 experiment on 180910

pixel_size = {'60X': 0.267, '100X': 0.16}
scale = pixel_size['60X']
base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180910_pACT1_mKate2/timelapse'
image_filename, bf_filename, fl_filename = '/180910_yFB11_12_mated_hap3_1_60X_5min_10lp_v1_', \
                                           'w1Brightfield confocal_', \
                                           'w2515 laser 10_'  # note to change fluorescence path to match laser power
fl_filename_c2 = 'w3594 laser 10_'  # secondary fluorescence channel
label_path = None  # if label_path is None then this experiment doesn't use labeling
manual_annotation = False  # if manual_annotation then we will use manual annotation to assign ambiguous pairs.
num_scenes = 8  # num_scenes should be the number of scenes to analyze + 1
num_frames = np.asarray([55, 70, 65, 66, 51, 66, 66])
# num_frames should be the number of frames + 1. Default is the same for
# each field of view.
num_frames_analyzed = 30  # number of analyzed timepoints for Whi5 localization
bkgd_scene = 8  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
analyzed_scene = 1  # which scene will be used to manually track Whi5 localization
threshold = 10000  # threshold for visualizing log or linear fluorescence data
drange = 65535.0  # image fluorescence maximum
size_thresh = False  # do we remove the largest cells?

C.populate_cells_all_scenes_2(temp_base_path=base_path, temp_expt_path=expt_path, temp_image_filename=image_filename, temp_fl_filename=fl_filename,
                                temp_num_scenes=num_scenes, temp_num_frames=num_frames, temp_bkgd_scene=bkgd_scene, temp_fl_filename_c2=fl_filename_c2)



# def segment_cell_volume(temp_base_path, temp_expt_path, temp_image_filename, temp_fl_filename_c2, temp_num_frames,
#                         temp_scene, temp_bkgd_scene):
#     with open(temp_base_path+temp_expt_path+ '/cells_scene_{0}.pkl', 'rb') as input:
#         c = pickle.load(input)
#     frame_list = [obj.frames for obj in c]
#             update_list = []
#             for obj in c:
#                 obj.add_fluor_placeholders()
#             for frame_num in range(1, temp_num_frames[temp_scene - 1]):
#                 temp = [(frame_num in temp1) for temp1 in frame_list]
#                 update_list.append(
#                     [i for i, e in enumerate(temp) if e != 0])  # gives the list of indices that have to be addressed
#                 # at each frame
#                 filename_fl = temp_image_filename + temp_fl_filename_c2 + 's{0}_t{1}.TIF'.format(str(temp_scene),
#                                                                                                               str(frame_num))
#                 filename_fl_bkgd = temp_image_filename + temp_fl_filename_c2 + \
#                                    's{0}_t{1}.TIF'.format(str(temp_bkgd_scene), str(frame_num))
#                 # format is z, y, x
#                 bkgd_im = io.imread(temp_base_path + temp_expt_path + filename_fl_bkgd)
#                 bkgd_im1 = np.zeros(bkgd_im.shape)
#                 temp1 = np.mean(bkgd_im, axis=0)
#                 # taking the mean with respect to the z axis. Do it this way since there doesn't seem
#                 # to be any systematic bias in that direction.
#                 for i0 in range(bkgd_im.shape[0]):
#                     bkgd_im1[i0, :, :] = temp1[:, :]
#                 del temp1, bkgd_im