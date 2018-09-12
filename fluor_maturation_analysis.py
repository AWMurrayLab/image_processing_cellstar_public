import numpy as np
import custom_image_toolkit as C
import pandas as pd
import os

# This script takes the output from a fluorescence maturation experiment, and returns the integrated fluorescence of
# all cells without contact with the walls.

# # ######################################### experiment_specific_information

# yFB7 on 180731 with csm
pixel_size = {'60X': 0.267, '100X': 0.16}
base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180731_csm_fluor_mat/timelapse'
image_filename, bf_filename, fl_filename = '/180731_yFB7_60X_5lp_timelapse_2min_csm_', \
                                           'w1Brightfield confocal_', \
                                           'w2515 laser 5_'  # note to change fluorescence path to match laser power
num_scenes = 7  # num_scenes should be the number of scenes to analyze + 1
num_frames = 41*np.ones(num_scenes, dtype=int)  # num_frames should be the number of frames + 1. Default is the same for
# each field of view.
bkgd_scene = num_scenes  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
drange = 65535.0  # image fluorescence maximum
bkgd_details = None  # alternative to account for the fact that some timelapses may not have a dedicated background
# scene

# # yFB7 on 180731 with cycloheximide
# pixel_size = {'60X': 0.267, '100X': 0.16}
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180731_cycloheximide_fluor_mat/timelapse'
# image_filename, bf_filename, fl_filename = '/180731_yFB7_60X_5lp_timelapse_2min_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 5_'  # note to change fluorescence path to match laser power
# num_scenes = 7  # num_scenes should be the number of scenes to analyze + 1
# num_frames = 41*np.ones(num_scenes, dtype=int)  # num_frames should be the number of frames + 1. Default is the same for
# # each field of view.
# bkgd_scene = None  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
# drange = 65535.0  # image fluorescence maximum
# bkgd_details = [6, [209, 399], [60, 250]]  # alternative to account for the fact that some timelapses may not have a dedicated background
# # scene


#  Analysis of experimental data.

# Step 0:
# testing whether the progress report already exists
if os.path.exists(base_path+expt_path+'/progress_report'):
    print base_path+expt_path+'/progress_report'
    temp_df = pd.read_csv(base_path+expt_path+'/progress_report', sep='\t', index_col=0)
else:
    d = {'run': [1, 1, 0], 'script': ['populate_cells_all_scenes_1', 'populate_cells_all_scenes_2_v2',
                                      'label_boundaries']}
    temp_df = pd.DataFrame(data=d)


# Step 1:
# producing the necessary lists of cells from the segmentation and tracking data. Assumes a file structure produced by
# make_image_dirs.py
if temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_1'].index[0], 'run'] == 0:
    # if this has not already been run
    C.populate_cells_all_scenes_1(base_path, expt_path, image_filename, bf_filename, num_scenes, num_frames)
    temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_1'].index[0], 'run'] = 1
    # update the progress report
    print 'Brightfield cell data assembled'
else:  # if this has already been run
    print 'Using previously assembled brightfield cell data'
    # temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_1'].index[0], 'run'] = 5

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')

# Step 2:
# integrating fluorescence data into the brightfield cell data generated above. Note that this allows you to not have a
# background scene
if temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_2_v2'].index[0], 'run'] == 0:
    # if this has not already been run
    C.populate_cells_all_scenes_2_v2(base_path, expt_path, image_filename, fl_filename, num_scenes, num_frames,
                                     bkgd_scene, temp_bkgd_details=bkgd_details)
    temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_2_v2'].index[0], 'run'] = 1
    # update the progress report
    print 'Fluorescence cell data assembled'
else:  # if this has already been run
    print 'Using previously assembled fluorescence cell data'
    # temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_2'].index[0], 'run'] = 2

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')

# Step 3:
# Mark any cells which intersect with the boundaries at any point so that they can be excluded from analysis
if temp_df.loc[temp_df[temp_df.script == 'label_boundaries'].index[0], 'run'] == 0:
    # if this has not already been run
    C.label_boundaries(base_path, expt_path, image_filename, bf_filename, num_scenes, num_frames)
    temp_df.loc[temp_df[temp_df.script == 'label_boundaries'].index[0], 'run'] = 1
    # update the progress report
    print 'Boundary cells labeled appropriately'
else:  # if this has already been run
    print 'Boundary cells already labeled'
    # temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_2'].index[0], 'run'] = 2

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')
