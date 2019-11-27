import numpy as np
import matplotlib.pyplot as plt
from skimage import io
import os
import cPickle as pickle


# # pACT1-mCherry experiment on 181114
# expt_id = '/181114_yFB78_Raff_125Gal'

# # yFB79 expt 181207 machine learning test
# expt_id = '/181207_yFB79_60X_Raff_125uMGal_ML_test'

# # yFB79 expt 181207
# expt_id = '/181207_yFB79_60X_Raff_125uMGal'

# # yFB79 expt 190417
# expt_id = '/181917_yFB79_60X_Raff_125uMGal'

# # yFB77 expt 181220
# expt_id = '/181220_yFB77_60X_Raff_125uMGal'

# yFB78 expt 190607
# expt_id = '/190607_yFB78_60X_Raff_125uMGal'

# # yFB110 expt 190629
# expt_id = '/190629_yFB110_60X_Raff_125uMGal'

# # yFB78 expt 190606
# expt_id = '/190606_yFB78_60X_Raff_125uMGal'

# # yFB78 expt 190725
# expt_id = '/190725_yFB78_60X_2Raff_125uMGal'

# yFB79 expt 190612_timelapse
expt_id = '/190612_yFB79_timelapse'

pickle_in = open("./expt_ids"+expt_id+'.pickle',"rb")
ep = pickle.load(pickle_in)  # this gives us the experimental parameters, so that we can load everything in an
base_path, expt_path = ep['base_path'], ep['expt_path']
num_frames_analyzed, num_frames = ep['num_frames_analyzed'], ep['num_frames']
image_filename, fl_filename = ep['image_filename'], ep['fl_filename']
drange, threshold = ep['drange'], ep['threshold']



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
        np.save(directory+'/fl_loc_centres/'+'scene_{0}_frame_{1}'.format(scene, frame_num), temp)
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


coords = []
global ax, temp_im, fig, scene, x_pts, y_pts, directory, frame_num, directory1, temp_scenes1, temp_scenes, temp_editing
temp_editing = 0
# scene = analyzed_scene  # we want to get data for the first scene
directory = base_path+expt_path+'/whi5_analysis'
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
temp_new = np.random.randint(low=0, high=np.sum(num_frames)+1, size=1)  # randomly selecting which scenes we consider
while len(temp_scenes) < num_frames_analyzed:  # we
    if not(temp_new[0] in temp_scenes):  # making sure we don't have any repeats
        temp_scenes.append(temp_new[0])
    temp_new = np.random.randint(low=0, high=np.sum(num_frames)+1, size=1)
indexing = [int(np.sum(num_frames[:i0])) for i0 in range(0, len(num_frames))]
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
    scene = np.searchsorted(indexing, temp_ind, side='right')
    frame_num = temp_ind - indexing[scene-1] + 1
    directory1 = base_path+expt_path+'/scene_{0}/outputs'.format(scene)
    print 'Scene {0}, frame {1}'.format(scene, frame_num)
    outlines = np.load(directory1+'/cell_outlines_scene_{0}.npy'.format(scene))
    temp_im = io.imread(base_path+expt_path+image_filename+fl_filename+
                        's{0}_t{1}.TIF'.format(str(scene), str(frame_num)))
    temp_im1 = temp_im/drange
    if np.sum(temp_im > threshold) > 0:
        temp_im1 = np.log(np.amax(temp_im1, axis=0)/np.amax(temp_im1))
        # using a log scale in this case because this image contains the frustrating high energy pixels that sometimes
        # arise
    else:
        temp_im1 = np.amax(temp_im1, axis=0) / np.amax(temp_im1)
    temp_im1 *= outlines[frame_num-1, :, :] == 0
    temp_im1 += outlines[frame_num-1, :, :].astype('uint16')
    if os.path.exists(directory+'/fl_loc_centres/'+'scene_{0}_frame_{1}.npy'.format(scene, frame_num)):  # if we already
        # have data on this frame then we do not need
        temp_vals = np.load(directory+'/fl_loc_centres/'+'scene_{0}_frame_{1}.npy'.format(scene, frame_num))
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
