import numpy as np
import custom_image_toolkit as C
import pandas as pd
import os

# This script takes the output from a fluorescence maturation experiment, and returns the integrated fluorescence of
# all cells without contact with the walls.

# # ######################################### experiment_specific_information

# # yFB7 on 180731 with csm
# pixel_size = {'60X': 0.267, '100X': 0.16}
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180731_csm_fluor_mat/timelapse'
# image_filename, bf_filename, fl_filename = '/180731_yFB7_60X_5lp_timelapse_2min_csm_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 5_'  # note to change fluorescence path to match laser power
# num_scenes = 7  # num_scenes should be the number of scenes to analyze + 1
# num_frames = 41*np.ones(num_scenes, dtype=int)  # num_frames should be the number of frames + 1. Default is the same for
# # each field of view.
# bkgd_scene = num_scenes  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
# drange = 65535.0  # image fluorescence maximum
# bkgd_details = None  # alternative to account for the fact that some timelapses may not have a dedicated background
# # scene

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


# # Expt with Laura 181015
# pixel_size = {'60X': 0.267, '100X': 0.16}
# image_filename, bf_filename, fl_filename = '/2018_10_15_yLB256_yLB365_yFJB71_cellasics_01_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2594 laser 20_'  # note to change fluorescence path to match laser power
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181015_spinning_disk/timelapse'
# scenes = [3, 5, 7, 9, 10, 11, 12, 13]
# num_scenes = len(scenes)
# bkgd_scene = 22  # the last scene is the background
# num_frames = 22*np.ones(bkgd_scene, dtype=int)
# drange = 65535.0  # image fluorescence maximum
# bf_base_name = '/2018_10_15_yLB256_yLB365_yFJB71_cellasics_01_w1Brightfield confocal'
# date = '/2018_10_15'
# bkgd_details=None
# analyzed_scene = 9  # which scene will be used to manually track whether this cell can be segregated from the bulk
# num_frames_analyzed=21
# threshold = 10000  # threshold for visualizing log or linear fluorescence data
# drange = 65535.0  # image fluorescence maximum
# label_path=None


# Expt 181204 yFB79
#CHX
# image_filename, bf_filename, fl_filename = '/181204_yFB79_60X_timelapse_CHX_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2594 laser 5_'  # note to change fluorescence path to match laser power
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181204_yFB79_fluor_maturation/CHX_timelapse'
# num_scenes = 5
# bkgd_scene = 5  # the last scene is the background
# num_frames = 47*np.ones(bkgd_scene, dtype=int)
# bf_base_name = '/181204_yFB79_60X_timelapse_CHX_w1Brightfield confocal'
# CSM
# image_filename, bf_filename, fl_filename = '/181204_yFB79_60X_timelapse_good_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2594 laser 5_'  # note to change fluorescence path to match laser power
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181204_yFB79_fluor_maturation/CSM_timelapse'
# num_scenes = 6
# bkgd_scene = 6  # the last scene is the background
# num_frames = 47*np.ones(bkgd_scene, dtype=int)
# bf_base_name = '/181204_yFB79_60X_timelapse_good_w1Brightfield confocal'

# Expt 181212 timelapse yFB79 fluor maturation
# # CSM
# image_filename, bf_filename, fl_filename = '/181212_yFB79_60X_dex_fluor_maturation_timelapse_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2594 laser 10_'  # note to change fluorescence path to match laser power
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181212_yFB79_CSM_dex_fluor_maturation/CSM_timelapse'
# num_scenes = 7
# bkgd_scene = 7  # the last scene is the background
# num_frames = 31*np.ones(bkgd_scene, dtype=int)
# bf_base_name = '/181212_yFB79_60X_dex_fluor_maturation_timelapse_w1Brightfield confocal'
# date = '/181212'
# CSM
# image_filename, bf_filename, fl_filename = '/181212_yFB79_60X_dex_fluor_maturation_timelapse1_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2594 laser 10_'  # note to change fluorescence path to match laser power
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181212_yFB79_CSM_dex_fluor_maturation/CHX_timelapse'
# num_scenes = 7
# bkgd_scene = 7  # the last scene is the background
# num_frames = 41*np.ones(bkgd_scene, dtype=int)
# bf_base_name = '/181212_yFB79_60X_dex_fluor_maturation_timelapse1_w1Brightfield confocal'
# date = '/181212'

# Expt 181213 timelapse yFB7 fluor maturation
# # CSM
image_filename, bf_filename, fl_filename = '/181212_yFB7_60X_dex_fluor_maturation_timelapse_', \
                                           'w1Brightfield confocal_', \
                                           'w2515 laser 5_'  # note to change fluorescence path to match laser power
base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181213_yFB7_fluor_maturation/CSM_timelapse'
num_scenes = 6
bkgd_scene = 6  # the last scene is the background
num_frames = 47*np.ones(bkgd_scene, dtype=int)
# bf_base_name = '/181212_yFB7_60X_dex_fluor_maturation_timelapse_w1Brightfield confocal'
# # CHX
# image_filename, bf_filename, fl_filename = '/181212_yFB7_60X_dex_fluor_maturation_CHX_timelapse_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 5_'  # note to change fluorescence path to match laser power
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181213_yFB7_fluor_maturation/CHX_timelapse'
# num_scenes = 7
# bkgd_scene = 7  # the last scene is the background
# num_frames = 47*np.ones(bkgd_scene, dtype=int)
# bf_base_name = '/181212_yFB7_60X_dex_fluor_maturation_CHX_timelapse_w1Brightfield confocal'
# Generic
pixel_size = {'60X': 0.267, '100X': 0.16}
drange = 65535.0  # image fluorescence maximum
bkgd_details=None
threshold = 10000  # threshold for visualizing log or linear fluorescence data
drange = 65535.0  # image fluorescence maximum
label_path=None



#  Analysis of experimental data.

# Step 0:
# testing whether the progress report already exists
if os.path.exists(base_path+expt_path+'/progress_report'):
    print base_path+expt_path+'/progress_report'
    temp_df = pd.read_csv(base_path+expt_path+'/progress_report', sep='\t', index_col=0)
else:
    d = {'run': [0, 0, 0], 'script': ['populate_cells_all_scenes_1', 'populate_cells_all_scenes_2_v2',
                                      'label_boundaries']}
    temp_df = pd.DataFrame(data=d)

# print type(scenes[0])
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

# Step 4:
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