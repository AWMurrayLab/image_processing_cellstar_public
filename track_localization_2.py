import numpy as np
import matplotlib.pyplot as plt
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import pandas as pd

# This script takes the .npy outputs from track_localization_1 and assigns these points to individual cells at the
# timepoint in question for the cell outputs from populate_cells_all_scenes_3.py

base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
scene = 1
directory = base_path+expt_path+'/scene_{0}/outputs'.format(scene)
fluor_name = '/180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w2515 laser 20_'
drange = 65535.0
# opening cell file
with open(directory + '/cells_fl_lab_scene_{0}.pkl'.format(scene), 'rb') as input:
    c = pickle.load(input)
    frame_list = [obj.frames for obj in c]
# loading cell outlines to confirm assignments
outlines = np.load(directory + '/cell_outlines_scene_{0}.npy'.format(scene))
for frame_num in range(1, 31):
    temp = np.load(directory+'/fl_loc_centres/scene_{0}_frame_{1}.npy'.format(scene, frame_num))  # tracked centroids
    # in format x, y
    coords = zip(temp[:, 0], temp[:, 1])
    # print coords, temp
    temp2 = [(frame_num in temp1) for temp1 in frame_list]
    update_list = [i for i, e in enumerate(temp2) if e != 0]  # gives the indices of each object in the current frame
    # coords now stores the coordinates of cell center points for the first frame
    c, assignments = C.assign_labels_3(c, update_list, coords, frame_num)
    C.save_object(c, directory + '/cells_scene_{0}_v1.pkl'.format(scene))

    # organizing figure data
    temp_im = io.imread(
        base_path + expt_path + fluor_name + 's{0}_t{1}.TIF'.format(str(scene), str(frame_num))) / drange
    # print temp_im.dtype
    # exit()
    temp_im1 = np.amax(temp_im, axis=0) / np.amax(temp_im)  # scaling for visualization purposes.
    # print np.amax(temp_im1)
    # exit()
    temp_im1 *= outlines[frame_num - 1, :, :] == 0
    temp_im1 += outlines[frame_num - 1, :, :].astype('uint16')
    # print c[ind].frames
    temp_cell_coords = [c[ind].position[frame_num-c[ind].frames[0]] for ind in update_list if c[ind].nuclear_whi5[
        frame_num-c[ind].frames[0]] == 1]  # store the centroids of the correct cells.
    # plotting figures
    fig = plt.figure(figsize=[10, 10], frameon=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(temp_im1)
    temp_x, temp_y = zip(*temp_cell_coords)
    plt.plot(temp_x, temp_y, 'x', color='r',linestyle='None')
    fig.subplots_adjust(bottom=0)
    fig.subplots_adjust(top=1)
    fig.subplots_adjust(right=1)
    fig.subplots_adjust(left=0)
    fig.savefig(directory+'/images/whi5_assignments_frame_{0}.tif'.format(frame_num))



