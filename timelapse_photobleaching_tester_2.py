import numpy as np
import custom_image_toolkit as C
import pandas as pd
import os
import cPickle as pickle
import skimage
from skimage import io
import scipy
import scipy.io as sio
import matplotlib.pyplot as plt

# expt_ids = ['/181917_yFB79_60X_Raff_125uMGal', '/190607_yFB78_60X_Raff_125uMGal', '/181114_yFB78_Raff_125Gal']
# expt_ids = ['/190607_yFB110_60X_Raff_125uMGal']
expt_ids = ['/190629_yFB110_60X_Raff_125uMGal']
# '/181207_yFB79_60X_Raff_125uMGal'  # appears to have a corrupted .mat file for Scene 1 Frame 45
# expt_ids = ['/181114_yFB78_Raff_125Gal']
for i0 in range(len(expt_ids)):
    pickle_in = open("./expt_ids"+expt_ids[i0]+'.pickle',"rb")
    ep = pickle.load(pickle_in)  # this gives us the experimental parameters, so that we can load everything in an
    # automated fashion
    print ep['num_frames']
    max_len = np.amax(ep['num_frames'])
    # av_brightness = []
    av_brightness = np.empty([max_len, len(ep['num_frames']), 2, 2])  # gives us our array to fill with data from each frame
    av_brightness[:]=np.nan
    bkgd_brightness = np.empty([max_len,  len(ep['num_frames'])+1, 2, 2])
    bkgd_brightness[:]=np.nan
    # 2 fluorescence channels, and also records the average and median fluorescence
    temp_dir = ep['base_path'] + ep['expt_path']

    for i2 in range(1, max_len):  # only iterating through to the number of scenes stated -1
        temp_bkgd1 = io.imread(
            temp_dir + ep['image_filename'] + ep['fl_filename'] + 's{0}_t{1}.TIF'.format(ep['bkgd_scene'], i2))
        bkgd_brightness[i2-1,-1, 0, 0] = np.mean(temp_bkgd1)
        bkgd_brightness[i2-1, -1, 0, 1] = np.median(temp_bkgd1)
        im_shape = temp_bkgd1.shape
        temp_bkgd2 = io.imread(
            temp_dir + ep['image_filename'] + ep['fl_filename_c2'] + 's{0}_t{1}.TIF'.format(ep['bkgd_scene'], i2))
        bkgd_brightness[i2-1, -1, 1, 0] = np.mean(temp_bkgd2)
        bkgd_brightness[i2-1, -1, 1, 1] = np.median(temp_bkgd2)
        for i1 in range(1, len(ep['num_frames']) + 1):  # iterating through the different frames
            # binary mask for cells
            if i2<ep['num_frames'][i1 - 1]:  # only do this part if there are this many frames in that scene
                temp_filename = ep['image_filename']+ep['bf_filename']+'s{0}_t{1}_segmentation.mat'.format(
                    str(i1), str(i2).zfill(2))
                temp_seg = sio.loadmat(temp_dir+'/scene_{0}/segments/'.format(i1) + temp_filename)['segments']>0
                temp_seg1 = np.tile(temp_seg, [im_shape[0], 1, 1])

                # importing fluorescence channel 2
                temp_im1 = io.imread(
                    temp_dir + ep['image_filename'] + ep['fl_filename_c2'] + 's{0}_t{1}.TIF'.format(i1, i2))
                bkgd_vals = temp_im1[np.nonzero(temp_seg1==0)]  # gives us the background values
                bkgd_med = np.median(bkgd_vals)
                bkgd_std = bkgd_med-np.percentile(bkgd_vals, 16)
                bkgd_brightness[i2-1,i1-1,1,0] = np.mean(bkgd_vals)
                bkgd_brightness[i2 - 1, i1 - 1, 1, 1] = bkgd_med
                temp_seg2 = temp_im1>bkgd_med+4.0*bkgd_std  # segmenting the actual cell volume
                temp_seg2 *= temp_seg1  # this gives us the cell segments in 3D

                # subtracting the background signal
                temp_im1 = temp_im1 - bkgd_med
                av_brightness[i2 - 1, i1 - 1, 1, 0] = np.mean(temp_im1[np.nonzero(temp_seg2)])
                av_brightness[i2 - 1, i1 - 1, 1, 1] = np.median(temp_im1[np.nonzero(temp_seg2)])

                # importing fluorescence channel 1
                temp_im1 = io.imread(temp_dir+ep['image_filename']+ep['fl_filename']+'s{0}_t{1}.TIF'.format(i1, i2))
                # subtracting the background signal
                bkgd_vals = temp_im1[np.nonzero(temp_seg1 == 0)]  # gives us the background values
                bkgd_med = np.median(bkgd_vals)
                bkgd_brightness[i2 - 1, i1 - 1, 0, 0] = np.mean(bkgd_vals)
                bkgd_brightness[i2 - 1, i1 - 1, 0, 1] = bkgd_med
                temp_im1 = temp_im1 - bkgd_med

                av_brightness[i2-1, i1-1, 0, 0] = np.mean(temp_im1[np.nonzero(temp_seg2)])
                av_brightness[i2-1, i1-1, 0, 1] = np.median(temp_im1[np.nonzero(temp_seg2)])

            print expt_ids[i0], 'Finished Scene {0}, Frame {1}'.format(i1, i2)
    np.save(temp_dir+'/photobleaching.npy', av_brightness)
    np.save(temp_dir+'/bkgd_brightness.npy', bkgd_brightness)