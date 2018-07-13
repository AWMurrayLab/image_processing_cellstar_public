import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import manual_image_toolkit as M
import os
import scipy
from scipy import stats

# this script takes a list of cell tracks with daughter lineages assigned (e.g. the output from assign_lineages_1.py)
# and creates a collection of reliable cell cycles from this by iterating through this list.
color=cm.tab10(np.linspace(0,1,10))
print color[0]
base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
with open(base_path+expt_path+'/cell_cycles_compiled.pkl', 'rb') as input:
    cc = pickle.load(input)
# initial definitions
fluor_name = '/180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w2515 laser 20_'
drange = 65535.0

plot_cycles=False
frame_nums = [31, 50, 50, 50, 50]

if plot_cycles:
# testing the correct breakdown of single cell cycles
    for scene in range(1,6):
        directory = base_path + expt_path + '/scene_{0}/outputs'.format(scene)
        scene_inds = [obj for obj in cc if obj.data_origin[1]==scene]  # picking out cells that are only from this scene
        frame_list = [obj.frames for obj in scene_inds]
        for frame_num in range(1, frame_nums[
                    scene - 1]):  # we go through frame by frame, starting with frame # 2 since we need a previous
            temp2 = [(frame_num in temp1) for temp1 in frame_list]
            update_list = [i for i, e in enumerate(temp2) if e != 0]  # cells to plot at each timepoint
            # loading image
            temp_im = io.imread(
                base_path + expt_path + fluor_name + 's{0}_t{1}.TIF'.format(str(scene), str(frame_num))) / drange
            temp_im1 = np.amax(temp_im, axis=0) / np.amax(temp_im)  # scaling for visualization purposes.

            fig = plt.figure(figsize=[5.12, 5.12], frameon=False)
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)
            ax.imshow(temp_im1, cmap='viridis')
            ax.set_xlim(left=0, right=512)
            ax.set_ylim(bottom=0, top=512)
            for i0 in update_list:
                temp_vals = scene_inds[i0].ellipse_fit[frame_num-scene_inds[i0].frames[0]]
                # print temp_vals, i0
                temp_coords = C.generate_ellipse_coords(temp_vals)
                temp_x, temp_y = zip(*temp_coords)
                plt.plot(temp_x, temp_y, color=color[scene_inds[i0].index % 10])

            fig.subplots_adjust(bottom=0)
            fig.subplots_adjust(top=1)
            fig.subplots_adjust(right=1)
            fig.subplots_adjust(left=0)
            fig.savefig(directory + '/images/cell_cycle_assignments_frame_{0}.tif'.format(frame_num))

# Testing the assignment of buds
plot_buds = True

filt_cc = [obj for obj in cc if obj.celltype != 2 and obj.error == False and obj.complete]


if plot_buds:
# testing the correct breakdown of single cell cycles
    for scene in range(1, 6):
        directory = base_path + expt_path + '/scene_{0}/outputs'.format(scene)
        scene_inds = [obj for obj in filt_cc if obj.data_origin[1] == scene]
        # picking out cells that are only from this scene and that are not buds
        frame_list = [obj.frames for obj in scene_inds]
        for frame_num in range(1, frame_nums[
                    scene - 1]):  # we go through frame by frame, starting with frame # 2 since we need a previous
            temp2 = [(frame_num in temp1) for temp1 in frame_list]
            update_list = [i for i, e in enumerate(temp2) if e != 0]  # cells to plot at each timepoint
            # loading image
            temp_im = io.imread(
                base_path + expt_path + fluor_name + 's{0}_t{1}.TIF'.format(str(scene), str(frame_num))) / drange
            temp_im1 = np.amax(temp_im, axis=0) / np.amax(temp_im)  # scaling for visualization purposes.

            fig = plt.figure(figsize=[5.12, 5.12], frameon=False)
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)
            ax.imshow(temp_im1, cmap='viridis')
            ax.set_xlim(left=0, right=512)
            ax.set_ylim(bottom=0, top=512)
            for i0 in update_list:
                temp_vals = scene_inds[i0].ellipse_fit[frame_num-scene_inds[i0].frames[0]]
                # print temp_vals, i0
                temp_coords = C.generate_ellipse_coords(temp_vals)
                temp_x, temp_y = zip(*temp_coords)
                plt.plot(temp_x, temp_y, color=color[scene_inds[i0].index % 10])
                if not(scene_inds[i0].daughter is None):
                    temp_vals = scene_inds[i0].ellipse_fit_bud[frame_num - scene_inds[i0].frames[0]]
                    # print temp_vals, i0
                    if temp_vals!=0:
                        temp_coords = C.generate_ellipse_coords(temp_vals)
                        temp_x, temp_y = zip(*temp_coords)
                        plt.plot(temp_x, temp_y, color=color[scene_inds[i0].index % 10])
            fig.subplots_adjust(bottom=0)
            fig.subplots_adjust(top=1)
            fig.subplots_adjust(right=1)
            fig.subplots_adjust(left=0)
            fig.savefig(directory + '/images/bud_assignments_frame_{0}.tif'.format(frame_num))
