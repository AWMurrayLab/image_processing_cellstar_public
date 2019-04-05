import numpy as np
import matplotlib.pyplot as plt
from skimage import io
import os




# yFB79 expt 181207 machine learning test
#
# expt_id = '/181207_yFB79_60X_Raff_125uMGal_ML_test'
# pixel_size = {'60X': 0.267, '100X': 0.16}
# scale = pixel_size['60X']
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/181207_yFB79_60X_Raff_125uMGal_ML_test/timelapse'
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
# print np.sum(num_frames)
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


# # ######################################### experiment_specific_information
date = '/190322'
expt_id = date+'_timepoint'
completed_scenes = [1, 0, 0, 1]  # which scenes I have finished image segmentation on
base_path, expt_path1 = '/scratch/lab/image_analysis_scratch', '/190322_yFB78_yFB79_CSM_Raff_Gal'
expt_conds = ['/yFB78_125uMGal', '/yFB78_800uMGal', '/yFB79_125uMGal', '/yFB79_800uMGal']
expt_paths = [expt_path1+cond for cond in expt_conds]
bkgd_scenes = [1, 1, 1, 101]  # the first scene is not always the bkgd
num_scenes = [101, 102, 102, 101]  # exact number of scenes to be analyzed
# number of frames including first background. Note some later frames have a second background
image_filenames = ['/190322_60X_yFB78_1XCSM_2Raff_125Gal_scene_',
                '/190322_60X_yFB78_1XCSM_2Raff_800Gal _scene_',
                '/190322_60X_yFB79_1XCSM_2XCSMpad_2Raff_125Gal_scene_',
                '/190322_60X_yFB79_1XCSM_2Raff_800Gal_better_timing_scene_']
bf_filename = '_w1Brightfield confocal'
fl_filename = '_w2515 laser 30'
fl_filename_c2 = '_w3594 laser 30'
date = '/190322'
pixel_size = {'60X': 0.267, '100X': 0.16}
zstep = 0.7  # distance between z steps
drange = 65535.0  # image fluorescence maximum
prog_vec = [0, 0]
num_frames_analyzed=30
threshold = 10000


def onclick(event):
    ix, iy = event.xdata, event.ydata
    print 'x = %d, y = %d' % (
        ix, iy)
    global coords
    if temp_editing == 0:  # if we are not in "editing mode"
        x_pts.append(ix)
        y_pts.append(iy)
        line.set_xdata(x_pts)
        line.set_ydata(y_pts)
        fig.canvas.draw()
        # coords.append(np.array([np.floor(ix), np.floor(iy)]).astype('int'))
        coords.append((np.floor(ix), np.floor(iy)))
        print 'Added point (x, y) =', coords[-1]
    else:  # if we are in "editing mode"
        temp_distances = np.linalg.norm(np.array(coords)-np.array([ix, iy]), axis=1)
        temp_min_ind = np.argmin(temp_distances)  # index of the coordinate corresponding to the point in question
        temp_coords = coords[temp_min_ind]
        del coords[temp_min_ind], x_pts[temp_min_ind], y_pts[temp_min_ind]
        line.set_xdata(x_pts)
        line.set_ydata(y_pts)
        fig.canvas.draw()
        print 'Removed point (x,y) =', temp_coords
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
        temp = np.zeros([len(x_vals), 2])
        temp[:, 0] = x_vals[:]
        temp[:, 1] = y_vals[:]
        # print temp
        np.save(directory+'/fl_loc_centres/'+'scene_{0}'.format(scene), temp)
    elif val == 'b':
        temp = coords.pop()
        x_pts.pop()
        y_pts.pop()
        line.set_xdata(x_pts)
        line.set_ydata(y_pts)
        fig.canvas.draw()
        print 'Removed last added point (x, y) =', temp
    elif val == 'h':
        print 'Left click to annotate a cell'
        print 'h = help'
        print 'q = quit'
        print 'f = move forward one frame'
        print 'b = remove the last seed'
        print 'e = enter editing mode'
    elif val == 'q':  # if we want to quit the analysis partway through
        print('Exiting')
        exit()
    elif val == 'e':  # if we want to enter "editing" mode
        global temp_editing
        temp_editing = np.mod(temp_editing + 1, 2)
        if temp_editing == 1:
            print 'Entered editing mode. Click on or near a point to remove it.'
        elif temp_editing == 0:
            print 'Exited editing mode. Click to add a cell.'
    return coords


for cond in range(len(expt_conds)):
    if completed_scenes[cond]:
        coords = []
        global ax, temp_im, fig, scene, x_pts, y_pts, directory, scene, directory1, temp_scenes1, temp_scenes, temp_editing
        temp_editing = 0
        # scene = analyzed_scene  # we want to get data for the first scene
        directory = base_path+expt_paths[cond]+'/outputs/whi5_analysis'
        if not os.path.exists(directory):
            os.makedirs(directory)
        if not os.path.exists(directory+'/images'):
            os.makedirs(directory+'/images')
        # selecting randomly which images we will consider from which scenes
        if os.path.exists(directory+'/samples.npy'):  # we pick up based on what has been done already
            temp_scenes = np.load(directory+'/samples.npy').tolist()
            temp_scenes1 = np.load(directory+'/completed_samples.npy')
            if len(temp_scenes) > num_frames_analyzed:
                raise ValueError('Cannot proceed since you are trying to analyze fewer frames than have already been analyzed')
        else:
            temp_scenes = []
            temp_scenes1 = np.zeros(num_frames_analyzed)
        temp_new = np.random.randint(low=0, high=num_scenes[cond]+1, size=1)  # randomly selecting which scenes we consider
        while len(temp_scenes) < num_frames_analyzed:  # we
            if not(temp_new[0] in temp_scenes) and not(temp_new[0] == bkgd_scenes[cond]):
                # making sure we don't have any repeats or sample the background
                temp_scenes.append(temp_new[0])
            temp_new = np.random.randint(low=0, high=num_scenes[cond] + 1,
                                         size=1)  # randomly selecting which scenes we consider
        print 'Keep going until scene number', num_frames_analyzed
        i0 = 0
        if np.sum(temp_scenes1 == 0) == 0:
            i1 = num_frames_analyzed # in this case we have nothing left to do
            print 'All frames already analyzed'
        else:
            i1 = np.nonzero(temp_scenes1 == 0)[0][0]
        for temp_ind in temp_scenes:
            print 'REACHED SCENE {0} OUT OF {1}'.format(temp_scenes.index(temp_ind), num_frames_analyzed)
            # picking which scene and frame we are considering here
            scene = temp_ind
            directory1 = base_path+expt_paths[cond]+'/outputs'
            print 'Condition '+expt_conds[cond]+', Scene {0}'.format(scene)
            outlines = np.load(directory1+'/cell_outlines.npy'.format(scene))
            temp_im = io.imread(base_path+expt_paths[cond]+image_filenames[cond]+fl_filename+
                                '_s{0}.TIF'.format(str(scene)))
            temp_im1 = temp_im/drange
            if np.sum(temp_im > threshold) > 0:
                temp_im1 = np.log(np.amax(temp_im1, axis=0) / np.amax(temp_im1))
                # using a log scale in this case because this image contains the frustrating high energy pixels that sometimes
                # arise
            else:
                temp_im1 = np.amax(temp_im1, axis=0) / np.amax(temp_im1)
            temp_im1 *= outlines[scene-1, :, :] == 0
            temp_im1 += outlines[scene-1, :, :].astype('uint16')

            if os.path.exists(directory+'/fl_loc_centres/'+'scene_{0}.npy'.format(scene)):  # if we already
                # have data on this frame then we do not need
                temp_vals = np.load(directory+'/fl_loc_centres/'+'scene_{0}.npy'.format(scene))
                coords = zip(*[temp_vals[:, 0], temp_vals[:, 1]])
                x_pts = list(temp_vals[:, 0])  # starting the points off correctly
                y_pts = list(temp_vals[:, 1])
                del temp_vals
            else:  # if we have no data on this frame
                coords = []
                x_pts = []
                y_pts = []
            # filename = '/scene_{0}/outputs/images/frame_{1}.tif'.format(scene, frame)
            # temp_im = io.imread(base_path+expt_path+filename)
            # print coords
            fig = plt.figure(figsize=[10, 10], frameon=False)
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)
            ax.imshow(temp_im1)
            # variable for the points
            line, = ax.plot(x_pts, y_pts, marker="o", linestyle='None', color='r')
            # ax.imshow(outlines[frame_num-1, :, :])
            cid = fig.canvas.mpl_connect('button_press_event', onclick)
            cid1 = fig.canvas.mpl_connect('key_press_event', keypress)
            plt.show(fig)
            temp_scenes1[i0] = 1  # indicating that we have completed that scene
            # print 'completed scene {0}, frame {1}, coordinates:'.format(scene, frame_num)
            # print 'Done {0} out of {1}'.format(i0, num_frames_analyzed)
            i0 += 1  # to track the next scene
            np.save(directory + '/completed_samples.npy', temp_scenes1)
            np.save(directory + '/samples.npy', temp_scenes)
