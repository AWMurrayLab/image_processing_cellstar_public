import numpy as np
import custom_image_toolkit as C
import pandas as pd
import os
import cPickle as pickle

# # ######################################### experiment_specific_information

# # yFB29 and 41 experiment on 180531
# base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
# image_filename, bf_filename, fl_filename = '/180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_', 'w1Brightfield confocal_', \
#                                            'w2515 laser 20_'  # note to change fluorescence path to match laser power
# num_scenes = 10  # num_scenes should be the number of scenes to analyze + 1
# num_frames = 41*np.ones(num_scenes, dtype=int)  # num_frames should be the number of frames + 1. Default is the same
# # for each field of view.
# bkgd_scene = num_scenes  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.

# # yFB29 and 41 experiment on 180725
# pixel_size = {'60X': 0.267, '100X': 0.16}
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180725_timelapse/timelapse'
# image_filename, bf_filename, fl_filename = '/180725_yFB29dyed_yFB41_60X_10lp_10min_timelapse_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 10_'  # note to change fluorescence path to match laser power
# label_path = '/180725_timelapse/initial_dyed_test'
# num_scenes = 11  # num_scenes should be the number of scenes to analyze + 1
# num_frames = 41*np.ones(num_scenes, dtype=int)  # num_frames should be the number of frames + 1. Default is the same for
# # each field of view.
# bkgd_scene = num_scenes  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
# analyzed_scene = 1  # which scene will be used to manually track Whi5 localization
# threshold = 10000  # threshold for visualizing log or linear fluorescence data
# drange = 65535.0  # image fluorescence maximum


# # yFB30 and 45 experiment on 180725
# pixel_size = {'60X': 0.267, '100X': 0.16}
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180802_staining_mix_expt/timelapse'
# image_filename, bf_filename, fl_filename = '/180802_yFB30dyed_yFB45_60X_10min__', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 10_'  # note to change fluorescence path to match laser power
# label_path = '/180802_staining_mix_expt/initial_staining'
# num_scenes = 10  # num_scenes should be the number of scenes to analyze + 1
# num_frames = 74*np.ones(num_scenes, dtype=int)  # num_frames should be the number of frames + 1. Default is the same for
# # each field of view.
# num_frames_analyzed = 40  # number of analyzed timepoints for Whi5 localization
# bkgd_scene = num_scenes  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
# analyzed_scene = 1  # which scene will be used to manually track Whi5 localization
# threshold = 10000  # threshold for visualizing log or linear fluorescence data
# drange = 65535.0  # image fluorescence maximum


# # yFB29 experiment on 180823
# pixel_size = {'60X': 0.267, '100X': 0.16}
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180823_yFB29_800uMGal_timelapse/timelapse'
# image_filename, bf_filename, fl_filename = '/180823_yFB29_60X_800uMGal_10lp_timelapse_settings_test_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 5_'  # note to change fluorescence path to match laser power
# label_path = None  # if label_path is None then this experiment doesn't use labeling
# num_scenes = 7  # num_scenes should be the number of scenes to analyze + 1
# num_frames = np.asarray([61, 51, 61, 51, 51, 51])
# # num_frames should be the number of frames + 1. Default is the same for
# # each field of view.
# num_frames_analyzed = 40  # number of analyzed timepoints for Whi5 localization
# bkgd_scene = 12  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
# analyzed_scene = 1  # which scene will be used to manually track Whi5 localization
# threshold = 10000  # threshold for visualizing log or linear fluorescence data
# drange = 65535.0  # image fluorescence maximum
# prog_vec = [1, 0, 0, 0, 0, 0, 0, 0]


# yFB29 experiment on 180831
# pixel_size = {'60X': 0.267, '100X': 0.16}
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180831_yFB29_800uMGal/timelapse'
# image_filename, bf_filename, fl_filename = '/180831_yFB29_800uMGal_60X_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 10_'  # note to change fluorescence path to match laser power
# scale = pixel_size['60X']
# height = 5.0  # minimum height of the cellasics chamber
# fl_filename_c2 = None  # secondary fluorescence channel
# label_path = None  # if label_path is None then this experiment doesn't use labeling
# manual_annotation = False  # if manual_annotation then we will use manual annotation to assign ambiguous pairs.
# num_scenes = 7  # num_scenes should be the number of scenes to analyze + 1
# num_frames = np.asarray([56, 56, 56, 56, 56, 56])
# # num_frames should be the number of frames + 1. Default is the same for
# # each field of view.
# num_frames_analyzed = 40  # number of analyzed timepoints for Whi5 localization
# bkgd_scene = 12  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
# analyzed_scene = 1  # which scene will be used to manually track Whi5 localization
# threshold = 10000  # threshold for visualizing log or linear fluorescence data
# drange = 65535.0  # image fluorescence maximum
# size_thresh = True  # do we remove the largest cells?
# prog_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0]

# # pACT1-mKate2 experiment on 180910
#
# expt_id = '/180910_yFB71'
# pixel_size = {'60X': 0.267, '100X': 0.16}
# scale = pixel_size['60X']
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180910_pACT1_mKate2/timelapse'
# image_filename, bf_filename, fl_filename = '/180910_yFB11_12_mated_hap3_1_60X_5min_10lp_v1_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 10_'  # note to change fluorescence path to match laser power
# fl_filename_c2 = 'w3594 laser 10_'  # secondary fluorescence channel
# fluor_c2 = True
# label_path = None  # if label_path is None then this experiment doesn't use labeling
# manual_annotation = False  # if manual_annotation then we will use manual annotation to assign ambiguous pairs.
# num_scenes = 8  # num_scenes should be the number of scenes to analyze + 1
# num_frames = np.asarray([55, 70, 65, 66, 51, 66, 66])
# # num_frames should be the number of frames + 1. Default is the same for
# # each field of view.
# num_frames_analyzed = 30  # number of analyzed timepoints for Whi5 localization
# bkgd_scene = 8  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
# analyzed_scene = 1  # which scene will be used to manually track Whi5 localization
# threshold = 10000  # threshold for visualizing log or linear fluorescence data
# drange = 65535.0  # image fluorescence maximum
# prog_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0]
# size_thresh = False  # do we remove the largest cells?
# zstep = 0.7  # distance between z steps
# tstep = 5.0

# # pACT1-mCherry experiment on 181114
#
# expt_id = '/181114_yFB78_Raff_125Gal'
# pixel_size = {'60X': 0.267, '100X': 0.16}
# scale = pixel_size['60X']
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181114_yFB78_Raff_125Gal/timelapse'
# image_filename, bf_filename, fl_filename = '/181114_yFB78_60X_2Raff_125Gal__', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 10_'  # note to change fluorescence path to match laser power
# fl_filename_c2 = 'w3594 laser 10_'  # secondary fluorescence channel
# fluor_c2 = True
# label_path = None  # if label_path is None then this experiment doesn't use labeling
# manual_annotation = False  # if manual_annotation then we will use manual annotation to assign ambiguous pairs.
# num_scenes = 13  # num_scenes should be the number of scenes to analyze + 1
# # num_frames = np.asarray([66, 66, 66, 60, 65, 56, 56, 62, 62, 56, 56, 56])
# num_frames = np.asarray([51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51])
# # num_frames should be the number of frames + 1. Default is the same for
# # each field of view.
# num_frames_analyzed = 30  # number of analyzed timepoints for Whi5 localization
# bkgd_scene = 13  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
# analyzed_scene = 1  # which scene will be used to manually track Whi5 localization
# threshold = 10000  # threshold for visualizing log or linear fluorescence data
# drange = 65535.0  # image fluorescence maximum
# prog_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0]
# size_thresh = False  # do we remove the largest cells?
# zstep = 0.7  # distance between z steps
# tstep = 10.0

# # yFB79 expt 181207
#
# expt_id = '/181207_yFB79_60X_Raff_125uMGal'
# pixel_size = {'60X': 0.267, '100X': 0.16}
# scale = pixel_size['60X']
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181207_yFB79_60X_Raff_125uMGal/timelapse'
# image_filename, bf_filename, fl_filename = '/181207_yFB79_60X_Raff_125uMGal_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 10_'  # note to change fluorescence path to match laser power
# fl_filename_c2 = 'w3594 laser 10_'  # secondary fluorescence channel
# fluor_c2 = True
# label_path = None  # if label_path is None then this experiment doesn't use labeling
# manual_annotation = False  # if manual_annotation then we will use manual annotation to assign ambiguous pairs.
# num_scenes = 15  # num_scenes should be the number of scenes to analyze + 1
# # num_frames = np.asarray([66, 66, 66, 60, 65, 56, 56, 62, 62, 56, 56, 56])
# num_frames = np.asarray([45, 55, 61, 56, 61, 61, 61, 56, 61, 51, 51, 51, 61, 61])
# # num_frames should be the number of frames + 1. Default is the same for
# # each field of view.
# num_frames_analyzed = 30  # number of analyzed timepoints for Whi5 localization
# bkgd_scene = 15  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
# analyzed_scene = 1  # which scene will be used to manually track Whi5 localization
# threshold = 10000  # threshold for visualizing log or linear fluorescence data
# drange = 65535.0  # image fluorescence maximum
# prog_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0]
# size_thresh = False  # do we remove the largest cells?
# zstep = 0.7  # distance between z steps
# tstep = 10.0


# yFB79 expt 181207 ML replicate

expt_id = '/181207_yFB79_60X_Raff_125uMGal'
pixel_size = {'60X': 0.267, '100X': 0.16}
scale = pixel_size['60X']
base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181207_yFB79_60X_Raff_125uMGal_ML_test/timelapse'
image_filename, bf_filename, fl_filename = '/181207_yFB79_60X_Raff_125uMGal_', \
                                           'w1Brightfield confocal_', \
                                           'w2515 laser 10_'  # note to change fluorescence path to match laser power
fl_filename_c2 = 'w3594 laser 10_'  # secondary fluorescence channel
fluor_c2 = True
label_path = None  # if label_path is None then this experiment doesn't use labeling
manual_annotation = False  # if manual_annotation then we will use manual annotation to assign ambiguous pairs.
num_scenes = 15  # num_scenes should be the number of scenes to analyze + 1
# num_frames = np.asarray([66, 66, 66, 60, 65, 56, 56, 62, 62, 56, 56, 56])
num_frames = np.asarray([45, 55, 61, 56, 61, 61, 61, 56, 61, 51, 51, 51, 61, 61])
# num_frames should be the number of frames + 1. Default is the same for
# each field of view.
num_frames_analyzed = 30  # number of analyzed timepoints for Whi5 localization
bkgd_scene = 15  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
analyzed_scene = 1  # which scene will be used to manually track Whi5 localization
threshold = 10000  # threshold for visualizing log or linear fluorescence data
drange = 65535.0  # image fluorescence maximum
prog_vec = [0, 0, 0, 0, 0, 0, 0, 0, 0]
size_thresh = False  # do we remove the largest cells?
zstep = 0.7  # distance between z steps
tstep = 10.0

# saving the experimental parameters so we can just load them in future.

expt_params = {'scale':scale, 'base_path':base_path, 'expt_path':expt_path, 'image_filename':image_filename,
               'bf_filename':bf_filename, 'fl_filename':fl_filename, 'fluor_c2':fluor_c2,
               'fl_filename_c2':fl_filename_c2, 'label_path':label_path, 'manual_annotation':False,
               'num_scenes':num_scenes, 'num_frames':num_frames, 'num_frames_analyzed':num_frames_analyzed,
               'bkgd_scene':bkgd_scene, 'analyzed_scene':analyzed_scene, 'expt_id':expt_id, 'threshold':threshold,
               'drange':drange, 'size_thresh':False, 'zstep':zstep, 'tstep':tstep}
pickle_out = open('./expt_ids'+expt_id+".pickle","wb")
pickle.dump(expt_params, pickle_out)
pickle_out.close()

#  Analysis of experimental data.

# Step 0:
# testing whether the progress report already exists
if os.path.exists(base_path+expt_path+'/progress_report'):
    print base_path+expt_path+'/progress_report'
    temp_df = pd.read_csv(base_path+expt_path+'/progress_report', sep='\t', index_col=0)
else:
    d = {'run': prog_vec, 'script': ['populate_cells_all_scenes_1', 'populate_cells_all_scenes_2',
                                   'populate_cells_all_scenes_3', 'track_localization_manual_annotation',
                                         'analyze_whi5_distribution', 'assign_lineages', 'create_cycles',
                                                     'validate_cycles', 'filter_cycles']}
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
# integrating fluorescence data into the brightfield cell data generated above
if temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_2'].index[0], 'run'] == 0:
    # if this has not already been run
    C.populate_cells_all_scenes_2(base_path, expt_path, image_filename, fl_filename, num_scenes, num_frames, bkgd_scene,
                                  temp_fl_filename_c2=fl_filename_c2)
    temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_2'].index[0], 'run'] = 1
    # update the progress report
    print 'Fluorescence cell data assembled'
else:  # if this has already been run
    print 'Using previously assembled fluorescence cell data'
    # temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_2'].index[0], 'run'] = 2

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')

# Step 3:
# Adding labels for whether cells were dyed or not.
if temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_3'].index[0], 'run'] == 0 and not(label_path is None):
    # if this has not already been run
    C.populate_cells_all_scenes_3(base_path, expt_path, label_path, num_scenes)
    temp_df.loc[temp_df[temp_df.script == 'populate_cells_all_scenes_3'].index[0], 'run'] = 1
    # update the progress report
    print 'Dyed labeling cell data assembled'
else:  # if this has already been run
    print 'Using previously assembled dye labeling cell data'

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')

# Step 4:
# Tracking Whi5 localization with manual annotation of datasets. Assumes you have run track_localization_final.py with
# the above parameters (default is running this for scene 1).

if os.path.exists(base_path+expt_path+'/whi5_analysis/completed_samples.npy'):
    temp = np.load(base_path + expt_path + '/whi5_analysis/completed_samples.npy')
    if np.sum(temp == 0) == 0:
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
        print 'Unable to track Whi5 distribution: You must complete manually annotating a dataset first with track_localization_final.py. Exited and saved progress.'
        temp_df.to_csv(base_path + expt_path + '/progress_report', sep='\t')
        exit()
else:
    print 'Unable to track Whi5 distribution: You must manually annotate a dataset first with track_localization_final.py. Exited and saved progress.'
    temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')
    exit()

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')

# Step 5:
# Analyzing Whi5 localization via by a learned method based on the manual annotation of datasets.
if temp_df.loc[temp_df[temp_df.script == 'analyze_whi5_distribution'].index[0], 'run'] == 0:
    # if this has not already been run
    C.analyze_whi5_distribution(base_path, expt_path,image_filename, fl_filename, num_frames, analyzed_scene,
                              num_scenes, threshold, drange, temp_analyzed_frames=num_frames_analyzed,
                                temp_label_path=label_path)
    temp_df.loc[temp_df[temp_df.script == 'analyze_whi5_distribution'].index[0], 'run'] = 1
    # update the progress report
    print 'Whi5 localization data analyzed'
else:  # if this has already been run
    print 'Using previously analyzed Whi5 localization data'

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')

# Step 6:
# Assigning lineages based on Whi5 localization data
if temp_df.loc[temp_df[temp_df.script == 'assign_lineages'].index[0], 'run'] == 0:
    # if this has not already been run
    C.assign_lineages(base_path, expt_path,image_filename, fl_filename, num_frames, analyzed_scene,
                              num_scenes, threshold, drange, manual_annotation)
    temp_df.loc[temp_df[temp_df.script == 'assign_lineages'].index[0], 'run'] = 1
    # update the progress report
    print 'Lineage data analyzed'
else:  # if this has already been run
    print 'Using previously analyzed lineage data'

# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')

# Step 7:
# Creating individual cell cycles based on lineage data
if temp_df.loc[temp_df[temp_df.script == 'create_cycles'].index[0], 'run'] == 0:
    # if this has not already been run
    C.create_cycles_full(base_path, expt_path, num_scenes)
    temp_df.loc[temp_df[temp_df.script == 'create_cycles'].index[0], 'run'] = 1
    # update the progress report
    print 'Cycle data created'
else:  # if this has already been run
    print 'Using previously created cell cycle data'


# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')


# Step 8:
# filtering poor quality cell cycles
if temp_df.loc[temp_df[temp_df.script == 'filter_cycles'].index[0], 'run'] == 0:
    # if this has not already been run
    C.filter_cycles(base_path, expt_path, temp_scale=scale, temp_size_thresh=size_thresh)
    temp_df.loc[temp_df[temp_df.script == 'filter_cycles'].index[0], 'run'] = 1
    # update the progress report
    print 'Cycle data filtered'
else:  # if this has already been run
    print 'Using previously filtered cell cycle data'


# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')


# Step 9:
# Creating plots to validate cell cycle data
if temp_df.loc[temp_df[temp_df.script == 'validate_cycles'].index[0], 'run'] == 0:
    # if this has not already been run
    C.validate_cycles(base_path, expt_path, image_filename, bf_filename, num_frames, num_scenes)
    temp_df.loc[temp_df[temp_df.script == 'validate_cycles'].index[0], 'run'] = 1
    # update the progress report
    print 'Cycle data validation plots created'
else:  # if this has already been run
    print 'Cycle data validation plots already created'


# Saving the progress report
print 'Saving progress report'
temp_df.to_csv(base_path+expt_path+'/progress_report', sep='\t')


# After this runs to completion we can run plot_cycles with the same parameters as above.