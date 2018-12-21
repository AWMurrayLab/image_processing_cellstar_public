import numpy as np
import matplotlib.pyplot as plt
from skimage import io
import os

# experiment 180531
# base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
# fluor_name = '/180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w2515 laser 20_'

# # experiment 180725
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
# threshold=10000  # this is the threshold we use to decide whether to plot the log or linear version of the image
# drange = 65535.0  # image fluorescence maximum

# # yFB30 and 45 experiment on 180725
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180802_staining_mix_expt/timelapse'
# image_filename, bf_filename, fl_filename = '/180802_yFB30dyed_yFB45_60X_10min__', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 10_'  # note to change fluorescence path to match laser power
# label_path = '/180802_staining_mix_expt/initial_staining'
# num_scenes = 10  # num_scenes should be the number of scenes to analyze + 1
# num_frames = 73*np.ones(num_scenes, dtype=int)  # We will only analyze 30 frames because otherwise it will be too
# # onerous
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
# num_frames = np.asarray([74, 50, 65, 50, 50, 55])
# # num_frames should be the number of frames + 1. Default is the same for
# # each field of view.
# num_frames_analyzed = 40  # number of analyzed timepoints for Whi5 localization
# bkgd_scene = 12  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
# analyzed_scene = 1  # which scene will be used to manually track Whi5 localization
# threshold = 10000  # threshold for visualizing log or linear fluorescence data
# drange = 65535.0  # image fluorescence maximum


# # yFB29 experiment on 180831
# pixel_size = {'60X': 0.267, '100X': 0.16}
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180831_yFB29_800uMGal/timelapse'
# image_filename, bf_filename, fl_filename = '/180831_yFB29_800uMGal_60X_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 10_'  # note to change fluorescence path to match laser power
# label_path = None  # if label_path is None then this experiment doesn't use labeling
# num_scenes = 5  # num_scenes should be the number of scenes to analyze + 1
# num_frames = np.asarray([56, 56, 56, 56])
# # num_frames should be the number of frames + 1. Default is the same for
# # each field of view.
# num_frames_analyzed = 40  # number of analyzed timepoints for Whi5 localization
# bkgd_scene = 12  # number of the bkgd_scene. Set equal to 1 greater than the scenes analyzed by default.
# analyzed_scene = 1  # which scene will be used to manually track Whi5 localization
# threshold = 10000  # threshold for visualizing log or linear fluorescence data
# drange = 65535.0  # image fluorescence maximum
# prog_vec = [0, 0, 0, 0, 0, 0, 0, 0]


# pACT1-mKate2 experiment on 180910

# pixel_size = {'60X': 0.267, '100X': 0.16}
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180910_pACT1_mKate2/timelapse'
# image_filename, bf_filename, fl_filename = '/180910_yFB11_12_mated_hap3_1_60X_5min_10lp_v1_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2515 laser 10_'  # note to change fluorescence path to match laser power
# fl_filename_c2 = 'w3594 laser 10_'  # secondary fluorescence channel
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
# prog_vec = [0, 0, 0, 0, 0, 0, 0, 0]


# # Expt with Laura 181015
# pixel_size = {'60X': 0.267, '100X': 0.16}
# image_filename, bf_filename, fl_filename = '/2018_10_15_yLB256_yLB365_yFJB71_cellasics_01_', \
#                                            'w1Brightfield confocal_', \
#                                            'w2594 laser 20_'  # note to change fluorescence path to match laser power
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181015_spinning_disk/timelapse'
# scenes = [3, 5, 7, 9, 10, 11, 12, 13]
# label_path=None
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
# num_scenes = 7  # num_scenes should be the number of scenes to analyze + 1
# # num_frames = np.asarray([66, 66, 66, 60, 65, 56, 56, 62, 62, 56, 56, 56])
# num_frames = np.asarray([51, 51, 51, 51, 51, 51])
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


# yFB79 expt 181207

expt_id = '/181207_yFB79_60X_Raff_125uMGal'
pixel_size = {'60X': 0.267, '100X': 0.16}
scale = pixel_size['60X']
base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181207_yFB79_60X_Raff_125uMGal/timelapse'
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


def onclick(event):
    ix, iy = event.xdata, event.ydata
    print 'x = %d, y = %d'%(
        ix, iy)
    x_pts.append(ix)
    y_pts.append(iy)
    line.set_xdata(x_pts)
    line.set_ydata(y_pts)
    fig.canvas.draw()
    global coords
    # coords.append(np.array([np.floor(ix), np.floor(iy)]).astype('int'))
    coords.append((np.floor(ix), np.floor(iy)))
    return coords


def keypress(event):
    # global val
    val = event.key

    if val == 'f':
        # close image
        fig.canvas.mpl_disconnect(cid)
        fig.canvas.mpl_disconnect(cid1)
        plt.close(fig)
        if not os.path.exists(directory+'/fl_loc_centres'):
            os.makedirs(directory+'/fl_loc_centres')
        x_vals, y_vals = zip(*coords)
        temp = np.zeros([len(x_vals),2])
        temp[:, 0] = x_vals[:]
        temp[:, 1] = y_vals[:]
        np.save(directory+'/fl_loc_centres/'+'scene_{0}_frame_{1}'.format(scene, frame_num), temp)
        print 'completed scene {0}, frame {1}, coordinates:'.format(scene, frame_num), temp
    elif val == 'b':
        temp = coords.pop()
        x_pts.pop()
        y_pts.pop()
        line.set_xdata(x_pts)
        line.set_ydatba(y_pts)
        fig.canvas.draw()
        print temp
    # elif val == 'b':
    #     temp = coords.pop()
    #     x_pts.pop()
    #     y_pts.pop()
    #     line.set_xdata(x_pts)
    #     line.set_ydata(y_pts)
    #     fig.canvas.draw()
    #     print temp
    # elif val == 'q':
    return coords


coords = []
global ax, temp_im, fig, scene, x_pts, y_pts, directory, frame_num
scene = analyzed_scene  # we want to get data for the first scene
directory = base_path+expt_path+'/scene_{0}/outputs'.format(scene)
outlines = np.load(directory+'/cell_outlines_scene_{0}.npy'.format(scene))
coord_list = []
# for frame_num in range(1, num_frames[scene-1]):
#     coord_list.append([])
# print np.amax(outlines)
# exit()
print 'Keep going until frame number', num_frames_analyzed
for frame_num in range(1, num_frames_analyzed):
    temp_im = io.imread(base_path+expt_path+image_filename+fl_filename+
                        's{0}_t{1}.TIF'.format(str(scene), str(frame_num)))
    temp_im1=temp_im/drange
    # print temp_im.dtype
    # exit()
    # temp_mask = np.amax(temp_im,axis=0)  # this will be the threshold to filter out random maxed out pixels
    # temp_im = np.amax(temp_im, axis=0)
    if np.sum(temp_im>threshold)>0:
        temp_im1 = np.log(np.amax(temp_im1, axis=0)/np.amax(temp_im1))
        # using a log scale in this case because this image contains the frustrating high energy pixels that sometimes
        # crop up
    else:
        temp_im1 = np.amax(temp_im1, axis=0) / np.amax(temp_im1)
    # print np.amax(temp_im1)
    # exit()
    temp_im1*= outlines[frame_num-1, :, :]==0
    temp_im1+= outlines[frame_num-1, :, :].astype('uint16')
    coords = []
    # filename = '/scene_{0}/outputs/images/frame_{1}.tif'.format(scene, frame)
    # temp_im = io.imread(base_path+expt_path+filename)
    fig = plt.figure(figsize=[10, 10], frameon=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(temp_im1)
    # variable for the points
    x_pts = []
    y_pts = []
    line, = ax.plot(x_pts, y_pts, marker="o", linestyle='None', color='r')
    # ax.imshow(outlines[frame_num-1, :, :])
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    cid1 = fig.canvas.mpl_connect('key_press_event', keypress)
    plt.show(fig)
    # print 'completed scene {0}, frame {1}, coordinates:'.format(scene, frame_num), coords
    coord_list.append(coords)
coord_vals = np.zeros([sum([len(obj) for obj in coord_list]), 3])
ind = 0
ind1 = 1
print coord_list
for obj in coord_list:
    print 'here'
    print obj
    if len(obj)>0:  # if we were unable to find any values in this frame then there should be nothing to unpack
        x,y = zip(*obj)
        coord_vals[ind:ind + len(obj), 0] = ind1
        coord_vals[ind:ind + len(obj), 1] = x[:]
        coord_vals[ind:ind + len(obj), 2] = y[:]
    ind1 += 1
    ind += len(obj)
print len(coord_list), sum([len(obj) for obj in coord_list]), coord_vals
np.save(directory+'/fl_loc_cells_scene_{0}.npy'.format(scene), coord_vals)
