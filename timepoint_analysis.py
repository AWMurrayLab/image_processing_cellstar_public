import numpy as np
import scipy
import matplotlib.pyplot as plt
import pandas as pd
import os
import custom_image_toolkit as C
import cPickle as pickle

# to be run after timepoint_make_image_dirs.py

# # ######################################### experiment_specific_information
# # timepoint experiment 3/22/19
# date = '/190322'
# expt_id = date+'_timepoint'
# completed_scenes = [1, 1, 1, 1]  # which scenes I have finished image segmentation on
# base_path, expt_path1 = '/scratch/lab/image_analysis_scratch', '/190322_yFB78_yFB79_CSM_Raff_Gal'
# expt_conds = ['/yFB78_125uMGal', '/yFB78_800uMGal', '/yFB79_125uMGal', '/yFB79_800uMGal']
# expt_paths = [expt_path1+cond for cond in expt_conds]
# bkgd_scenes = [1, 1, 1, 101]  # the first scene is not always the bkgd
# num_scenes = [101, 101, 102, 101]  # exact number of scenes to be analyzed
# # number of frames including first background. Note some later frames have a second background
# image_filenames = ['/190322_60X_yFB78_1XCSM_2Raff_125Gal_scene_',
#                 '/190322_60X_yFB78_1XCSM_2Raff_800Gal _scene_',
#                 '/190322_60X_yFB79_1XCSM_2XCSMpad_2Raff_125Gal_scene_',
#                 '/190322_60X_yFB79_1XCSM_2Raff_800Gal_better_timing_scene_']
# bf_filename = '_w1Brightfield confocal'
# fl_filename = '_w2515 laser 30'
# fl_filename_c2 = '_w3594 laser 30'
# date = '/190322'
# pixel_size = {'60X': 0.267, '100X': 0.16}
# zstep = 0.7  # distance between z steps
# drange = 65535.0  # image fluorescence maximum
# prog_vec = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
# num_frames_analyzed=30
# threshold = 10000
# thresh_blobs=[0.007,0.05,0.007,0.007]
#
# # timepoint experiment 4/3/19
# date = '/190403'
# expt_id = date+'_timepoint'
# completed_scenes = [1, 1]  # which scenes I have finished image segmentation on
# base_path, expt_path1 = '/scratch/lab/image_analysis_scratch', '/190403_yFB78_yFB79_CSM_Raff_timepoint'
# expt_conds = ['/yFB78_125uMGal', '/yFB79_125uMGal']
# expt_paths = [expt_path1+cond for cond in expt_conds]
# bkgd_scenes = [1, 1]  # the first scene is not always the bkgd
# num_scenes = [101, 100]  # exact number of scenes to be analyzed
# # number of frames including first background. Note some later frames have a second background
# image_filenames = ['/190403_yFB78_2XCSM_2Raff_125uMGal_60X_timepoint',
#                 '/190403_yFB79_2XCSM_2Raff_125uMGal_60X_timepoint']
# bf_filename = '_w1Brightfield confocal'
# fl_filename = '_w2515 laser 30'
# fl_filename_c2 = '_w3594 laser 30'
# date = '/190403'
# pixel_size = {'60X': 0.267, '100X': 0.16}
# zstep = 0.7  # distance between z steps
# drange = 65535.0  # image fluorescence maximum
# prog_vec = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]  # one column for each condition
# num_frames_analyzed = 30
# threshold = 10000
# thresh_blobs = [0.003, 0.003]

# # timepoint experiment 4/17/19
# date = '/190417'
# expt_id = date+'_timepoint'
# completed_scenes = [1, 1]  # which scenes I have finished image segmentation on
# base_path, expt_path1 = '/scratch/lab/image_analysis_scratch', '/190417_yFB78_79_CSM_Raff_timepoint'
# expt_conds = ['/yFB78_125uMGal', '/yFB79_125uMGal']
# expt_paths = [expt_path1+cond for cond in expt_conds]
# bkgd_scenes = [1, 1]  # the first scene is not always the bkgd
# num_scenes = [100, 100]  # exact number of scenes to be analyzed
# # number of frames including first background. Note some later frames have a second background
# image_filenames = ['/190417_yFB78_60X_timepoint_2XCSM_2Raff_125uMGal',
#                 '/190417_yFB79_60X_timepoint_2XCSM_2Raff_125uMGal']
# bf_filename = '_w1Brightfield confocal'
# fl_filename = '_w2515 laser 30'
# fl_filename_c2 = '_w3594 laser 30'
# date = '/190417'
# pixel_size = {'60X': 0.267, '100X': 0.16}
# zstep = 0.7  # distance between z steps
# drange = 65535.0  # image fluorescence maximum
# prog_vec = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]  # one column for each condition
# num_frames_analyzed = 30
# threshold = 10000
# thresh_blobs = [0.003, 0.003]

# # timepoint experiment 800uM Gal 6/7/19
# date = '/190607'
# expt_id = date+'_timepoint'
# completed_scenes = [1]  # which scenes I have finished image segmentation on
# base_path, expt_path1 = '/scratch/lab/image_analysis_scratch', '/190607_yFB78_800uMGal_timepoint'
# expt_conds = ['/yFB78_800uMGal']
# expt_paths = [expt_path1+cond for cond in expt_conds]
# bkgd_scenes = [1]  # the first scene is not always the bkgd
# num_scenes = [101]  # exact number of scenes to be analyzed
# # number of frames including first background. Note some later frames have a second background
# image_filenames = ['/190607_yFB78_60X_2XCSM_2Raff_800uMGal']
# bf_filename = '_w1Brightfield confocal'
# fl_filename = '_w2515 laser 30'
# fl_filename_c2 = '_w3594 laser 30'
# date = '/190607'
# pixel_size = {'60X': 0.267, '100X': 0.16}
# zstep = 0.7  # distance between z steps
# drange = 65535.0  # image fluorescence maximum
# prog_vec = [[0, 0, 0, 0, 0]]  # one column for each condition
# num_frames_analyzed = 30
# threshold = 10000
# thresh_blobs = [0.05]

# timepoint experiment 800uM Gal 10/9/19
date = '/191009'
expt_id = date+'_timepoint'
completed_scenes = [1,1]  # which scenes I have finished image segmentation on
base_path, expt_path1 = '/scratch/lab/image_analysis_scratch', '/191009_yFB78_yFB79_timepoint'
expt_conds = ['/yFB78_800uMGal','/yFB79_800uMGal']
expt_paths = [expt_path1+cond for cond in expt_conds]
bkgd_scenes = [101,101]  # the first scene is not always the bkgd
num_scenes = [101,101]  # exact number of scenes to be analyzed
# number of frames including first background. Note some later frames have a second background
image_filenames = ['/191009_yFB78_800uMGal_60X','/191009_yFB79_800uMGal_60X1']
bf_filename = '_w1Brightfield confocal'
fl_filename = '_w2515 laser 30'
fl_filename_c2 = '_w3594 laser 30'
pixel_size = {'60X': 0.267, '100X': 0.16}
zstep = 0.7  # distance between z steps
drange = 65535.0  # image fluorescence maximum
prog_vec = [[0, 0, 0, 0, 0],[0, 0, 0, 0, 0]]  # one column for each condition
num_frames_analyzed = 30
threshold = 10000
thresh_blobs = [0.01,0.003]


# saving the experimental parameters so we can just load them in future.



expt_params = {'base_path':base_path, 'expt_path':expt_paths, 'image_filename':image_filenames,
               'bf_filename':bf_filename, 'fl_filename':fl_filename, 'fl_filename_c2':fl_filename_c2,
               'fl_filename_c2':fl_filename_c2, 'manual_annotation':False,
               'num_scenes':num_scenes, 'num_frames_analyzed':num_frames_analyzed,
               'bkgd_scenes':bkgd_scenes, 'expt_id':expt_id, 'drange':drange, 'zstep':zstep}
pickle_out = open('./expt_ids'+expt_id+".pickle","wb")
pickle.dump(expt_params, pickle_out)
pickle_out.close()

#  Analysis of experimental data.

# Step 0:
# testing whether the progress report already exists
if os.path.exists(base_path+expt_path1+'/progress_report'):
    # print base_path+expt_path1+'/progress_report'
    temp_df = pd.read_csv(base_path+expt_path1+'/progress_report', sep='\t', index_col=0)
else:
    d = {'script': ['timepoint_populate_cells_all_scenes_1', 'timepoint_populate_cells_all_scenes_2',
                                     'timepoint_track_localization_manual_annotation',
                                     'timepoint_analyze_whi5_distribution', 'timepoint_pandas_export']}
    for i0 in range(len(prog_vec)):  # giving the number of different conditions
        d['run_{0}'.format(i0)] = prog_vec[i0]
        # print prog_vec[i0]
    # print d
    temp_df = pd.DataFrame(data=d)

# Step 1:
# producing the necessary lists of cells from the segmentation and tracking data. Assumes a file structure produced by
# make_image_dirs.py
if (temp_df[temp_df.script == 'timepoint_populate_cells_all_scenes_1']==0).any().any():  # if there are any conditions
    # that have not been completed at this timepoint
    # if this has not already been run
    temp_df = C.timepoint_populate_cells_all_scenes_1(base_path, expt_paths, expt_conds, image_filenames, bf_filename, completed_scenes,
                                  num_scenes, temp_df)
    # temp_df.loc[temp_df[temp_df.script == 'timepoint_populate_cells_all_scenes_1'].index[0], 'run'] = 1
    # update the progress report
    print 'Brightfield cell data assembled'
else:  # if this has already been run
    print 'Using previously assembled brightfield cell data'
    # temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_1'].index[0], 'run'] = 5

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path1+'/progress_report', sep='\t')

# Step 2:
# integrating fluorescence data into the brightfield cell data generated above
if (temp_df[temp_df.script == 'timepoint_populate_cells_all_scenes_2']==0).any().any():  # if there are any conditions
    # if this has not already been run
    C.timepoint_populate_cells_all_scenes_2(base_path, expt_paths, expt_conds, image_filenames, bf_filename, completed_scenes,
                                  num_scenes, bkgd_scenes, fl_filename, fl_filename_c2, temp_thresh=thresh_blobs,
                                            temp_df=temp_df)
    # temp_df.loc[temp_df[temp_df.script == 'timepoint_populate_cells_all_scenes_2'].index[0], 'run'] = 1
    # update the progress report
    print 'Fluorescence cell data assembled'
else:  # if this has already been run
    print 'Using previously assembled fluorescence cell data'
    # temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_2'].index[0], 'run'] = 2

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path1+'/progress_report', sep='\t')

# Step 3:
# Tracking Whi5 localization with manual annotation of datasets. Assumes you have run track_localization_final.py with
# the above parameters (default is running this for scene 1).

if (temp_df[temp_df.script == 'timepoint_track_localization_manual_annotation']==0).any().any():  # if there are any conditions
    # if this has not already been run
    C.timepoint_track_localization_manual_annotation(base_path, expt_paths, expt_conds, image_filenames,
                                                     completed_scenes, fl_filename, temp_threshold=threshold,
                                                     temp_drange=drange, temp_df=temp_df)
    # temp_df.loc[temp_df[temp_df.script == 'timepoint_track_localization_manual_annotation'].index[0], 'run'] = 1
    # update the progress report
    print 'Manually annotated Whi5 localization data assembled'
else:  # if this has already been run
    print 'Using previously assembled Manually annotated Whi5 localization'

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path1+'/progress_report', sep='\t')

# Step 4:
# Analyzing Whi5 localization via by a learned method based on the manual annotation of datasets.
if (temp_df[temp_df.script == 'timepoint_analyze_whi5_distribution']==0).any().any():  # if there are any conditions
    # if this has not already been run
    C.timepoint_analyze_whi5_distribution(base_path, expt_paths, expt_conds, image_filenames, completed_scenes,
                                          fl_filename, num_scenes, threshold, drange, temp_df)
    # temp_df.loc[temp_df[temp_df.script == 'timepoint_analyze_whi5_distribution'].index[0], 'run'] = 1
    # update the progress report
    print 'Whi5 localization data analyzed'
else:  # if this has already been run
    print 'Using previously analyzed Whi5 localization data'

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path1+'/progress_report', sep='\t')


# Step 5:
# Exporting all relevant data as a pandas dataframe for easy analysis
if (temp_df[temp_df.script == 'timepoint_pandas_export']==0).any().any():  # if there are any conditions
    # if this has not already been run
    C.timepoint_pandas_export(base_path, expt_paths, expt_conds, expt_id, completed_scenes, temp_df)
    # temp_df.loc[temp_df[temp_df.script == 'timepoint_pandas_export'].index[0], 'run'] = 1
    # update the progress report
    print 'Data exported to pandas'
else:  # if this has already been run
    print 'Data already exported'

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path1+'/progress_report', sep='\t')

C.timepoint_plot_distribution(base_path, expt_paths, expt_conds, image_filenames, completed_scenes,
                                          fl_filename, num_scenes, threshold, drange)
