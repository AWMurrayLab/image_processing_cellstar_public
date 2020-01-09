import os
import shutil


# # yFB79 timelapse expt 181207_yFB79
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181207_yFB79_60X_Raff_125uMGal/timelapse'
# num_scenes = 15
# bkgd_scene = 15  # the last scene is the background
# num_frames = 73
# bf_base_name = '/181207_yFB79_60X_Raff_125uMGal_w1Brightfield confocal'
# date = '/181207'

# # yFB79 CSM timelapse expt 190417_yFB79
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/190417_yFB79_timelapse/timelapse'
# num_scenes = 19
# bkgd_scene = 19  # the last scene is the background
# num_frames = 70
# bf_base_name = '/190417_yFB79_60X_timelapse_2XCSM_2Raff_125uMGal_10min_w1Brightfield confocal'
# date = '/190417'

# yFB79 CSM timelapse expt 190612_yFB79
base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/190612_yFB79_timelapse/timelapse'
num_scenes = 15
bkgd_scene = 15  # the last scene is the background
num_frames = 65
bf_base_name = '/190612_yFB79_60X_2XCSM_2Raff_125uMGal__w1Brightfield confocal'
date = '/190612'

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

# # yFB78 CSM timelapse expt 190725_yFB78
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/190725_yFB78_timelapse/timelapse'
# num_scenes = 17
# bkgd_scene = 18  # the last scene is the background
# num_frames = 65
# bf_base_name = '/190725_yFB78_60X_2XCSM_2Raff_125uMGal_w1Brightfield confocal'
# date = '/190725'

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
