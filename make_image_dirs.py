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


# # Expt 180831 overexpression Whi5
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180910_pACT1_mKate2/timelapse'
# num_scenes = 8
# bkgd_scene = 8  # the last scene is the background
# num_frames = 73
# bf_base_name = '/180910_yFB11_12_mated_hap3_1_60X_5min_10lp_v1_w1Brightfield confocal'
# date = '/180910'


# # Expt with Laura 181015
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181015_spinning_disk/timelapse'
# num_scenes = 22
# bkgd_scene = 22  # the last scene is the background
# num_frames = 69
# bf_base_name = '/2018_10_15_yLB256_yLB365_yFJB71_cellasics_01_w1Brightfield confocal'
# date = '/2018_10_15'

# # Diploid experiment 181004
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181004_diploid_pACT1_mKate2/timelapse'
# num_scenes = 11
# bkgd_scene = 11  # the last scene is the background
# num_frames = 91
# bf_base_name = '/181004_yFB71_72_diploid_pACT1_mKate2_60X_5min_10lp_D1_w1Brightfield confocal'
# date = '/181004'

# # Haploid expt 181115_yFB78
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181114_yFB78_Raff_125Gal/timelapse'
# num_scenes = 13
# bkgd_scene = 13  # the last scene is the background
# num_frames = 70
# bf_base_name = '/181114_yFB78_60X_2Raff_125Gal__w1Brightfield confocal'
# date = '/181114'

# # Fluorophore maturation expt 181204_yFB79
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181204_yFB79_fluor_maturation/CHX_timelapse'
# num_scenes = 5
# bkgd_scene = 5  # the last scene is the background
# num_frames = 46
# bf_base_name = '/181204_yFB79_60X_timelapse_CHX_w1Brightfield confocal'
# date = '/181204'

# # Fluorophore maturation expt 181204_yFB79
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181204_yFB79_fluor_maturation/CSM_timelapse'
# num_scenes = 6
# bkgd_scene = 6  # the last scene is the background
# num_frames = 46
# bf_base_name = '/181204_yFB79_60X_timelapse_good_w1Brightfield confocal'
# date = '/181204'

# # yFB79 timelapse expt 181207_yFB79
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181207_yFB79_60X_Raff_125uMGal/timelapse'
# num_scenes = 15
# bkgd_scene = 15  # the last scene is the background
# num_frames = 73
# bf_base_name = '/181207_yFB79_60X_Raff_125uMGal_w1Brightfield confocal'
# date = '/181207'

# # yFB79 timelapse expt 181207_yFB79
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181212_yFB79_CSM_dex_fluor_maturation/CSM_timelapse'
# num_scenes = 7
# bkgd_scene = 7  # the last scene is the background
# num_frames = 31
# bf_base_name = '/181212_yFB79_60X_dex_fluor_maturation_timelapse_w1Brightfield confocal'
# date = '/181212'

# # yFB79 timelapse expt 181207_yFB79
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181212_yFB79_CSM_dex_fluor_maturation/CHX_timelapse'
# num_scenes = 7
# bkgd_scene = 7  # the last scene is the background
# num_frames = 40
# bf_base_name = '/181212_yFB79_60X_dex_fluor_maturation_timelapse1_w1Brightfield confocal'
# date = '/181212'

# # yFB7 CHX timelapse expt 181213_yFB7
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181213_yFB7_fluor_maturation/CHX_timelapse'
# num_scenes = 7
# bkgd_scene = 7  # the last scene is the background
# num_frames = 46
# bf_base_name = '/181212_yFB7_60X_dex_fluor_maturation_CHX_timelapse_w1Brightfield confocal'
# date = '/181212'

# # yFB7 CSM timelapse expt 181213_yFB7
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181213_yFB7_fluor_maturation/CSM_timelapse'
# num_scenes = 6
# bkgd_scene = 6  # the last scene is the background
# num_frames = 46
# bf_base_name = '/181212_yFB7_60X_dex_fluor_maturation_timelapse_w1Brightfield confocal'
# date = '/181212'

# # yFB7 CSM timelapse expt 181213_yFB7
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/190417_yFB79_timelapse/timelapse'
# num_scenes = 19
# bkgd_scene = 19  # the last scene is the background
# num_frames = 70
# bf_base_name = '/190417_yFB79_60X_timelapse_2XCSM_2Raff_125uMGal_10min_w1Brightfield confocal'
# date = '/181212'

# # yFB77 CSM timelapse expt 181220_yFB77
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181220_yFB77_CSM_Raff_Gal/timelapse'
# num_scenes = 15
# bkgd_scene = 15  # the last scene is the background
# num_frames = 67
# bf_base_name = '/181220_yFB77_60X_Raff_Gal_timelapse_w1Brightfield confocal'
# date = '/181220'

# # yFB78 CSM timelapse expt 190607_yFB78
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/190607_yFB78_timelapse/timelapse'
# num_scenes = 13
# bkgd_scene = 13  # the last scene is the background
# num_frames = 60
# bf_base_name = '/190607_yFB78_60X_2XCSM_2Raff_125uMGal2_w1Brightfield confocal'
# date = '/190607'

# # yFB78 CSM timelapse expt 190606_yFB78
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/190606_yFB78_timelapse/timelapse'
# num_scenes = 13
# bkgd_scene = 13  # the last scene is the background
# num_frames = 74
# bf_base_name = '/190606_yFB78_60X_2XCSM_2Raff_125uMGal_timelapse_w1Brightfield confocal'
# date = '/190606'

# # yFB110 CSM timelapse expt 190629_yFB110
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/190629_yFB110_timelapse/timelapse'
# num_scenes = 18
# bkgd_scene = 18  # the last scene is the background
# num_frames = 60
# bf_base_name = '/190629_yFB110_60X_2XCSM_2Raff_125uMGal_w1Brightfield confocal'
# date = '/190629'

# # yFB78 CSM timelapse expt 190725_yFB78
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/190725_yFB78_timelapse/timelapse'
# num_scenes = 17
# bkgd_scene = 18  # the last scene is the background
# num_frames = 65
# bf_base_name = '/190725_yFB78_60X_2XCSM_2Raff_125uMGal_w1Brightfield confocal'
# date = '/190725'

# yFB79 CSM timelapse expt 190612_yFB79
base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/190612_yFB79_timelapse/timelapse'
num_scenes = 15
bkgd_scene = 15  # the last scene is the background
num_frames = 65
bf_base_name = '/190612_yFB79_60X_2XCSM_2Raff_125uMGal__w1Brightfield confocal'
date = '/190612'


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