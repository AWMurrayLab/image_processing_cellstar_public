import numpy as np
import custom_image_toolkit as C
import pandas as pd
import os

# Expt with Laura 181015
pixel_size = {'60X': 0.267, '100X': 0.16}
image_filename, bf_filename, fl_filename = '/2018_10_15_yLB256_yLB365_yFJB71_cellasics_01_', \
                                           'w1Brightfield confocal_', \
                                           'w2594 laser 20_'  # note to change fluorescence path to match laser power
base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181015_spinning_disk/timelapse'
scenes = [3, 5, 7, 9, 10, 11, 12, 13]
num_scenes = len(scenes)
bkgd_scene = 22  # the last scene is the background
num_frames = 22*np.ones(bkgd_scene, dtype=int)
drange = 65535.0  # image fluorescence maximum
bf_base_name = '/2018_10_15_yLB256_yLB365_yFJB71_cellasics_01_w1Brightfield confocal'
date = '/2018_10_15'
bkgd_details=None
analyzed_scene = 9  # which scene will be used to manually track whether this cell can be segregated from the bulk
num_frames_analyzed=21
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
    d = {'run': [1, 1, 1, 0, 0], 'script': ['populate_cells_all_scenes_1', 'populate_cells_all_scenes_2_v2',
                                      'track_localization_manual_annotation', 'analyze_whi5_distribution',
                                            'label_boundaries']}
    temp_df = pd.DataFrame(data=d)

# print type(scenes[0])
# Step 1:
# producing the necessary lists of cells from the segmentation and tracking data. Assumes a file structure produced by
# make_image_dirs.py
if temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_1'].index[0], 'run'] == 0:
    # if this has not already been run
    C.populate_cells_all_scenes_1(base_path, expt_path, image_filename, bf_filename, num_scenes, num_frames,
                                  temp_sel_scenes=scenes)
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
                                     bkgd_scene, temp_bkgd_details=bkgd_details, temp_sel_scenes=scenes)
    temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_2_v2'].index[0], 'run'] = 1
    # update the progress report
    print 'Fluorescence cell data assembled'
else:  # if this has already been run
    print 'Using previously assembled fluorescence cell data'
    # temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_2'].index[0], 'run'] = 2

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')

# Step 3
# Tracking fluorescence localization with manual annotation of datasets. Assumes you have run track_localization_final.py with
# the above parameters (default is running this for scene 1). Note that cells will be labeled as having Whi5 'nuclear
# localized'
if os.path.exists(base_path+expt_path+'/scene_{0}/outputs/fl_loc_cells_scene_{0}.npy'.format(analyzed_scene)):
    # If we have manually annotated datasets
    if temp_df.loc[temp_df[temp_df.script == 'track_localization_manual_annotation'].index[0], 'run'] == 0:
        # if this has not already been run
        C.track_localization_manual_annotation(base_path, expt_path, image_filename, fl_filename, num_frames,
                                               analyzed_scene, temp_threshold=threshold, temp_drange=drange,
                                               temp_analyzed_frames=num_frames_analyzed, temp_label_path=label_path)
        temp_df.loc[temp_df[temp_df.script == 'track_localization_manual_annotation'].index[0], 'run'] = 1
        # update the progress report
        print 'Manually annotated Whi5 localization data assembled'
    else:  # if this has already been run
        print 'Using previously assembled Manually annotated Whi5 localization'
else:
    print 'Unable to track Whi5 distribution: You must manually annotate a dataset first. Exited and saved progress.'
    temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')
    exit()


# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')

# Step 4:
# Analyzing Whi5 localization via by a learned method based on the manual annotation of datasets.
if temp_df.loc[temp_df[temp_df.script == 'analyze_whi5_distribution'].index[0], 'run'] == 0:
    # if this has not already been run
    C.analyze_whi5_distribution(base_path, expt_path,image_filename, fl_filename, num_frames, analyzed_scene,
                              num_scenes, threshold, drange, temp_analyzed_frames=num_frames_analyzed,
                                temp_label_path=label_path, temp_sel_scenes=scenes)
    temp_df.loc[temp_df[temp_df.script == 'analyze_whi5_distribution'].index[0], 'run'] = 1
    # update the progress report
    print 'Whi5 localization data analyzed'
else:  # if this has already been run
    print 'Using previously analyzed Whi5 localization data'

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')


# Step 5:
# Mark any cells which intersect with the boundaries at any point so that they can be excluded from analysis
if temp_df.loc[temp_df[temp_df.script == 'label_boundaries'].index[0], 'run'] == 0:
    # if this has not already been run
    C.label_boundaries(base_path, expt_path, image_filename, bf_filename, num_scenes, num_frames, temp_sel_scenes=scenes)
    temp_df.loc[temp_df[temp_df.script == 'label_boundaries'].index[0], 'run'] = 1
    # update the progress report
    print 'Boundary cells labeled appropriately'
else:  # if this has already been run
    print 'Boundary cells already labeled'
    # temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_2'].index[0], 'run'] = 2

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')