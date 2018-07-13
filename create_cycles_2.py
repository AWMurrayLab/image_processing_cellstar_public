import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import manual_image_toolkit as M
import os
import scipy
from scipy import stats

# this script takes a list of cell tracks with daughter lineages assigned (e.g. the output from assign_lineages_1.py)
# and creates a collection of reliable cell cycles from this by iterating through this list.

base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
cc = []
ind = 0  # this indexes the cell we are considering in the current timestep
# this will be repeated with each new iteration

for scene in range(1, 8):
    data_index = [base_path+expt_path, scene]
    directory = base_path+expt_path+'/scene_{0}/outputs'.format(scene)
    with open(directory + '/cells_scene_{0}_v3.pkl'.format(scene), 'rb') as input:
        c = pickle.load(input)
    for i0 in range(len(c)):
        c, cc, ind = C.create_cycles(c, i0, cc, ind, temp_data_origin=data_index)  # note that the unique index for each
        # cell cycle is unique across scenes also, but that the
    for i0 in range(len(c)):
        c, cc = C.stitch_cycles(c, cc, i0)
    # cc += cc
cc = C.inherit_lineage_properties(cc)  # inheriting the lineage properties of related cells
print len(cc), np.sum([1 for obj in cc if obj.complete and not(obj.daughter is None)])
print np.sum([1 for obj in cc if obj.complete and not(obj.daughter is None) and obj.celltype == 1])
print np.sum([1 for obj in cc if obj.complete and not(obj.daughter is None) and obj.celltype == 1 and not obj.error])
print np.sum([1 for obj in cc if obj.complete and not(obj.daughter is None) and obj.celltype == 0])
print np.sum([1 for obj in cc if obj.complete and not(obj.daughter is None) and obj.celltype == 0 and not obj.error])
cc = C.integrate_bud_data(cc)
print len(cc), np.sum([1 for obj in cc if obj.complete and not(obj.daughter is None)])
print np.sum([1 for obj in cc if obj.complete and not(obj.daughter is None) and obj.celltype == 1])
print np.sum([1 for obj in cc if obj.complete and not(obj.daughter is None) and obj.celltype == 1 and not obj.error])
print np.sum([1 for obj in cc if obj.complete and not(obj.daughter is None) and obj.celltype == 0])
print np.sum([1 for obj in cc if obj.complete and not(obj.daughter is None) and obj.celltype == 0 and not obj.error])
# print cc[-1].ellipse_fit[0]

# print [obj.index for obj in cc][len(cc)-20:len(cc)-10], range(len(cc))[len(cc)-20:len(cc)-10], ind
    # raise ValueError('something is not right')

C.save_object(cc, base_path + expt_path + '/cell_cycles_compiled.pkl'.format(scene))
# cc = C.filter_data(cc)  # this filters the dataset to


