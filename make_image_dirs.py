import os
import shutil


# # Expt 180725
# base_path, expt_path = '/scratch/lab/image_analysis', '/180725_timelapse/timelapse'
# num_scenes = 11
# bkgd_scene = 11  # the last scene is the background
# num_frames = 73
# bf_base_name = '/180725_yFB29dyed_yFB41_60X_10lp_10min_timelapse_w1Brightfield confocal'
# date = '/180725'

# # Expt 180802
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180802_staining_mix_expt/timelapse'
# num_scenes = 10
# bkgd_scene = 10  # the last scene is the background
# num_frames = 73
# bf_base_name = '/180802_yFB30dyed_yFB45_60X_10min__w1Brightfield confocal'
# date = '/180802'

# # Expt 180731 csm
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180731_csm_fluor_mat/timelapse'
# num_scenes = 7
# bkgd_scene = 7  # the last scene is the background
# num_frames = 40
# bf_base_name = '/180731_yFB7_60X_5lp_timelapse_2min_csm_w1Brightfield confocal'
# date = '/180731'

# # Expt 180731 cycloheximide
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180731_cycloheximide_fluor_mat/timelapse'
# num_scenes = 7
# bkgd_scene = 7  # the last scene is the background
# num_frames = 40
# bf_base_name = '/180731_yFB7_60X_5lp_timelapse_2min_w1Brightfield confocal'
# date = '/180731'


# # Expt 180823 overexpression Whi5
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180823_yFB29_800uMGal_timelapse/timelapse'
# num_scenes = 12
# bkgd_scene = 12  # the last scene is the background
# num_frames = 73
# bf_base_name = '/180823_yFB29_60X_800uMGal_10lp_timelapse_settings_test_w1Brightfield confocal'
# date = '/180823'


# # Expt 180831 overexpression Whi5
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180831_yFB29_800uMGal/timelapse'
# num_scenes = 12
# bkgd_scene = 12  # the last scene is the background
# num_frames = 55
# bf_base_name = '/180831_yFB29_800uMGal_60X_w1Brightfield confocal'
# date = '/180823'


# Expt 180831 overexpression Whi5
base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180910_pACT1_mKate2/timelapse'
num_scenes = 8
bkgd_scene = 8  # the last scene is the background
num_frames = 73
bf_base_name = '/180910_yFB11_12_mated_hap3_1_60X_5min_10lp_v1_w1Brightfield confocal'
date = '/180910'


bkgd_extension = '_s{0}_t{1}.TIF'.format(bkgd_scene, 1)

for scene in range(1, num_scenes):
    directory = base_path+expt_path +'/scene_{0}'.format(scene)
    if not os.path.exists(directory):
        os.makedirs(directory)
    for frame_num in range(1, num_frames+1):
        extension = '_s{0}_t{1}.TIF'.format(scene, frame_num)
        extension1 = '_s{0}_t{1}.TIF'.format(scene, str(frame_num).zfill(2))
        # print base_path+expt_path+bf_base_name+extension
        if os.path.exists(base_path+expt_path+bf_base_name+extension):
        #     print 'I made it'
        #     exit()
            shutil.move(base_path+expt_path+bf_base_name+extension, directory+bf_base_name+extension1)
        # moving the files into each new directory
    shutil.copyfile(base_path+expt_path+bf_base_name+bkgd_extension, directory+date+'_bkgd'+bkgd_extension)