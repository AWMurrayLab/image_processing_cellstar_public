import numpy as np
import scipy
import scipy.optimize
import scipy.io as sio
import weakref
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from scipy import stats
import pandas as pd
import cPickle as pickle
import custom_image_toolkit as C
from skimage import io

# This script takes the cell output from populate_cells_2.py and adds data showing whether cells are labeled as being
# stained or not.

timestep = 10
pixel_size = {'60X': 0.267, '100X': 0.16}
base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
label_path = '/180531_dye_mix_expt/initial_dyed_test'
color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
frame_num=1
for scene in range(1, 8):
    print 'scene = {0}'.format(scene)
    with open(base_path+expt_path+'/scene_{0}/outputs/cells_fl_scene_{0}.pkl'.format(scene), 'rb') as input:
        c = pickle.load(input)
    frame_list = [obj.frames for obj in c]
    temp = [(frame_num in temp1) for temp1 in frame_list]  # gives the indices of each object in the first frame
    update_list = [i for i, e in enumerate(temp) if e != 0]
    tracking_csv_file = pd.DataFrame.from_csv(base_path + label_path + '/scene_{0}.csv'.format(scene), index_col=None)
    coords=[]
    for i0 in range(len(tracking_csv_file['X'])):
        coords.append(np.array([tracking_csv_file['X'][i0], tracking_csv_file['Y'][i0]]))
    # coords now stores the coordinates of cell center points for the first frame
    c = C.assign_labels_2(c, update_list, coords)
    directory = base_path + expt_path + '/scene_{0}/outputs'.format(scene)
    C.save_object(c, directory+'/cells_fl_lab_scene_{0}.pkl'.format(scene))

