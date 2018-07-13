import numpy as np
import matplotlib.pyplot as plt
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import manual_image_toolkit as M
import os
import scipy
from scipy import stats

base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
with open(base_path+expt_path+'/cell_cycles_compiled.pkl', 'rb') as input:
    cc = pickle.load(input)
# initial definitions

temp_cc = [obj for obj in cc if obj.complete and not(obj.daughter is None) and obj.celltype==1]

for obj in temp_cc:
    print obj.nuclear_whi5

