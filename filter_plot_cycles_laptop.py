import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import manual_image_toolkit as M
import os
import scipy
from scipy import stats

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

pixel_size = {'60X': 0.267, '100X': 0.16}
drange = 65535.0
# this script takes a list of cell tracks with daughter lineages assigned (e.g. the output from assign_lineages_1.py)
# and creates a collection of reliable cell cycles from this by iterating through this list.
color = cm.tab10(np.linspace(0, 1, 10))

# importing the cell cycle variables

# pACT1-mKate2 experiment on 180910
scale = pixel_size['60X']
# base_path, expt_path = '/scratch/lab/image_analysis_scratch', '/180910_pACT1_mKate2/timelapse'
base_path, expt_path = '/mnt/d/Lab/image_analysis', '/180831_yFB29_800uMGal'  # path for the laptop storage of images
script_path = '/mnt/c/Users/felix/Documents/GitHub/image_processing_cellstar'
timestep = 5.0
fluor_c2 = True

directory = base_path + expt_path + '/plots'
if not os.path.exists(directory):
    os.makedirs(directory)
with open(base_path+expt_path+'/cell_cycles_compiled.pkl', 'rb') as input:
    cc = pickle.load(input)
print 'I got here'
