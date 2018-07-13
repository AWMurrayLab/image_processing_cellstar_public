import os
from shutil import copyfile

base_path, date = '/scratch/lab/image_analysis', '/180619'
expt_path = date+'_dye_mix_expt'
directory1 = base_path+expt_path+'/timelapse'

num_scenes, num_timepoints = 11, 50
basename_bf = date+'_yFB43stain_yFB29unstain_60X_initial_staining_945pm_15min_10lp_w1Brightfield confocal'
basename_515 = date+'_yFB43stain_yFB29unstain_60X_initial_staining_945pm_15min_10lp_w2515 laser 10'

# first check that all files are labeled as expected
for scene in range(1, num_scenes):
    for frame_num in range(1, num_timepoints):
        init_ind = '_s{0}_t{1}.TIF'.format(scene, frame_num)
        if not os.path.exists(directory1 + basename_bf + init_ind):
            print 'bf ' + init_ind
        if not os.path.exists(directory1 + basename_515 + init_ind):
            print 'fl ' + init_ind
for scene in range(1, num_scenes+1):
    if scene == num_scenes:  # background frame
        if not os.path.exists(directory1 + '/bkgd'):
            os.makedirs(directory1 + '/bkgd')
        if not os.path.exists(directory1 + '/bkgd/fluorescence_data'):
            os.makedirs(directory1 + '/bkgd/fluorescence_data')
        for frame_num in range(1, num_timepoints):
            init_ind = '_s{0}_t{1}.TIF'.format(scene, frame_num)
            final_ind = '_s{0}_t{1}.TIF'.format(str(scene).zfill(2), str(frame_num).zfill(3))
            # renaming bf data
            if not os.path.exists(directory1 + '/bkgd' + basename_bf + final_ind):
                os.rename(directory1 + basename_bf + init_ind,
                          directory1 + '/bkgd' + basename_bf + final_ind)
            # renaming fl data
            if not os.path.exists(directory1 + '/bkgd/fluorescence_data' + basename_515 + final_ind):
                os.rename(directory1 + basename_515 + init_ind,
                          directory1 + '/bkgd/fluorescence_data' + basename_515 + final_ind)
                # print 'Scene {0}, Frame {1}'.format(scene, frame_num)
    else:
        if not os.path.exists(directory1 + '/scene_{0}'.format(scene)):
            os.makedirs(directory1 + '/scene_{0}'.format(scene))
            print
        if not os.path.exists(directory1 + '/scene_{0}/fluorescence_data'.format(scene)):
            os.makedirs(directory1 + '/scene_{0}/fluorescence_data'.format(scene))
        for frame_num in range(1, num_timepoints):
            init_ind = '_s{0}_t{1}.TIF'.format(scene, frame_num)
            final_ind = '_s{0}_t{1}.TIF'.format(str(scene).zfill(2), str(frame_num).zfill(3))
            # renaming bf data
            if not os.path.exists(directory1+'/scene_{0}'.format(scene) + basename_bf + final_ind):
                os.rename(directory1 + basename_bf + init_ind,
                          directory1+'/scene_{0}'.format(scene) + basename_bf + final_ind)
            # renaming fl data
            if not os.path.exists(directory1 + '/scene_{0}/fluorescence_data'.format(scene) + basename_515 + final_ind):
                os.rename(directory1 + basename_515 + init_ind,
                      directory1 + '/scene_{0}/fluorescence_data'.format(scene) + basename_515 + final_ind)
            # print 'Scene {0}, Frame {1}'.format(scene, frame_num)

# placing 1 background image in each bf scene so that we have something to compare them to
frame_num = 1
for scene in range(1, num_scenes):
    final_ind = '_s{0}_t{1}.TIF'.format(str(num_scenes).zfill(2), str(frame_num).zfill(3))
    initial_bkgd_name = basename_bf + final_ind
    final_bkgd_name = date+'_bkgd.TIF'
    output = directory1 + '/scene_{0}'.format(scene)+final_bkgd_name
    if not os.path.exists(output):
        copyfile(directory1+'/bkgd'+initial_bkgd_name, output)

