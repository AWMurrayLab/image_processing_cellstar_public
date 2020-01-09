import os
import shutil

# # yFB78 and yFB79 CSM timepoint expt 190322
# base_path, expt_path1 = '/scratch/lab/image_analysis_scratch', '/190322_yFB78_yFB79_CSM_Raff_Gal'
# expt_conds = ['/yFB78_125uMGal', '/yFB78_800uMGal', '/yFB79_125uMGal', '/yFB79_800uMGal']
# bkgd_scene = [1, 1, 1, 101]  # the first scene is the bkgd
# num_scenes = [101, 102, 102, 101]
# # number of frames including first background. Note some later frames have a second background
# bf_base_name = ['/190322_60X_yFB78_1XCSM_2Raff_125Gal_scene__w1Brightfield confocal',
#                 '/190322_60X_yFB78_1XCSM_2Raff_800Gal _scene__w1Brightfield confocal',
#                 '/190322_60X_yFB79_1XCSM_2XCSMpad_2Raff_125Gal_scene__w1Brightfield confocal',
#                 '/190322_60X_yFB79_1XCSM_2Raff_800Gal_better_timing_scene__w1Brightfield confocal']
# date = '/190322'

# # yFB78 and yFB79 CSM timepoint expt 190403
# base_path, expt_path1 = '/scratch/lab/image_analysis_scratch', '/190403_yFB78_yFB79_CSM_Raff_timepoint'
# expt_conds = ['/yFB78_125uMGal', '/yFB79_125uMGal']
# bkgd_scene = [1, 1]  # the first scene is the bkgd
# num_scenes = [101, 100]
# # number of frames including first background. Note some later frames have a second background
# bf_base_name = ['/190403_yFB78_2XCSM_2Raff_125uMGal_60X_timepoint_w1Brightfield confocal',
#                 '/190403_yFB79_2XCSM_2Raff_125uMGal_60X_timepoint_w1Brightfield confocal']
# date = '/190403'
#
# # yFB78 and yFB79 CSM timepoint expt 190417
# base_path, expt_path1 = '/scratch/lab/image_analysis_scratch', '/190417_yFB78_79_CSM_Raff_timepoint'
# expt_conds = ['/yFB78_125uMGal', '/yFB79_125uMGal']
# bkgd_scene = [1, 1]  # the first scene is the bkgd
# num_scenes = [100, 100]
# # number of frames including first background. Note some later frames have a second background
# bf_base_name = ['/190417_yFB78_60X_timepoint_2XCSM_2Raff_125uMGal_w1Brightfield confocal',
#                 '/190417_yFB79_60X_timepoint_2XCSM_2Raff_125uMGal_w1Brightfield confocal']
# date = '/190417'

# # yFB78 CSM timepoint expt 190607
# base_path, expt_path1 = '/scratch/lab/image_analysis_scratch', '/190607_yFB78_800uMGal_timepoint'
# expt_conds = ['/yFB78_800uMGal']
# bkgd_scene = [1]  # the first scene is the bkgd
# num_scenes = [101]
# # number of frames including first background. Note some later frames have a second background
# bf_base_name = ['/190607_yFB78_60X_2XCSM_2Raff_800uMGal_w1Brightfield confocal']
# date = '/190417'

for i0 in range(len(expt_conds)):
    expt_path = expt_path1+expt_conds[i0]
    bkgd_extension = '_s{0}.TIF'.format(bkgd_scene[i0])
    directory = base_path+expt_path+'/bf_ims'
    print directory
    if not os.path.exists(directory):
        os.makedirs(directory)
        print 'hi'
    for scene in range(1, num_scenes[i0]+1):
        extension = '_s{0}.TIF'.format(scene)
        extension1 = '_s{0}.TIF'.format(str(scene).zfill(3))
        # print base_path+expt_path+bf_base_name+extension
        if os.path.exists(base_path+expt_path+bf_base_name[i0]+extension):
            shutil.copyfile(base_path+expt_path+bf_base_name[i0]+extension, directory+bf_base_name[i0]+extension1)
        # moving the files into each new directory
        shutil.copyfile(base_path+expt_path+bf_base_name[i0]+bkgd_extension, directory+date+'_0_bkgd.TIF')
        if os.path.exists(directory+date+'_0_bkgd'):
            os.remove(directory + date + '_0_bkgd')
        if os.path.exists(directory+date+'bkgd'):
            os.remove(directory + date + 'bkgd')
