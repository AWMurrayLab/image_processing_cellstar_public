import numpy as np
import matplotlib
import pandas as pd
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import os

pixel_size = {'60X': 0.267, '100X': 0.16}
drange = 65535.0

# # pACT1-mKate2 experiment on 180910
# expt_id = '/180910_yFB71'

# # pACT1-mCherry experiment on 181114
# expt_id = '/181114_yFB78_Raff_125Gal'

# yFB79 Raffinose experiment on 181207
expt_id = '/181207_yFB79_60X_Raff_125uMGal'

pickle_in = open("./expt_ids"+expt_id+'.pickle',"rb")
ep = pickle.load(pickle_in)
print ep
directory = ep['base_path'] + ep['expt_path'] + '/plots'
if not os.path.exists(directory):
    os.makedirs(directory)
with open(ep['base_path'] + ep['expt_path']+'/cell_cycles_filtered.pkl', 'rb') as input:
    cc = pickle.load(input)
cols = ['type', '$t_{G1}$', '$t_{budded}$', '$V_{b,ell}$', '$V_{b,seg}$', '$V_{s,ell}$', '$V_{s,seg}$',
        '$V_{d,ell}$', '$V_{d,seg}$', '$V_{bud,ell}$', '$V_{bud,seg}$',
        '$F1_{b,zproj}$', '$F1_{b,seg}$', '$F1_{s,zproj}$', '$F1_{s,seg}$', '$F1_{d,zproj}$', '$F1_{d,seg}$', '$F1_{bud,zproj}$',
        '$F1_{bud,seg}$']
celltype = ['M', 'D', 'Bud', 'NA']
# 0 indexes the start of the cell cycle, obj.start is the index of Start, -1 is the end of the cell cycle.
if ep['fluor_c2']:  # if there is a second fluorescence channel
    cols += ['$F2_{b,zproj}$', '$F2_{b,seg}$', '$F2_{s,zproj}$', '$F2_{s,seg}$', '$F2_{d,zproj}$', '$F2_{d,seg}$',
             '$F2_{bud,zproj}$', '$F2_{bud,seg}$', '$F1_{b,nuc,int}$', '$F1_{b,nuc,av}$', '$F1_{b,cyt,int}$',
             '$F1_{b,cyt,av}$', '$F2_{b,nuc,int}$', '$F2_{b,nuc,av}$', '$F2_{b,cyt,int}$', '$F2_{b,cyt,av}$' ]
    data_vals = [[celltype[obj.celltype], obj.start * ep['tstep'], (len(obj.frames) - obj.start) * ep['tstep'],
                  obj.vb * ep['scale'] ** 3, obj.pixel_thresh_vol[0] * ep['scale'] ** 2 * ep['zstep'],
                  obj.ellipse_volume[obj.start] * ep['scale'] ** 3,
                  obj.pixel_thresh_vol[obj.start] * ep['scale'] ** 2 * ep['zstep'],
                  obj.ellipse_volume[-1] * ep['scale'] ** 3, obj.pixel_thresh_vol[-1] * ep['scale'] ** 2 * ep['zstep'],
                  obj.vbud[-1] * ep['scale'] ** 3, obj.segment_vol_bud[-1] * ep['scale'] ** 2 * ep['zstep'],
                  obj.zproj_fl[0], obj.pixel_thresh_fluor_vals[0], obj.zproj_fl[obj.start],
                  obj.pixel_thresh_fluor_vals[obj.start], obj.zproj_fl[-1], obj.pixel_thresh_fluor_vals[-1],
                  obj.segment_fl_bud[-1], obj.zproj_fl_bud[-1],
                  obj.zproj_fl_c2[0], obj.pixel_thresh_fluor_vals_c2[0], obj.zproj_fl_c2[obj.start],
                  obj.pixel_thresh_fluor_vals_c2[obj.start], obj.zproj_fl_c2[-1], obj.pixel_thresh_fluor_vals_c2[-1],
                  obj.segment_fl_bud_c2[-1], obj.zproj_fl_bud_c2[-1], obj.nuclear_fluor_int[0]]
                    for obj in cc]
else:
    data_vals = [[celltype[obj.celltype], obj.start * ep['tstep'], (len(obj.frames) - obj.start) * ep['tstep'],
                  obj.vb * ep['scale'] ** 3, obj.pixel_thresh_vol[0] * ep['scale'] ** 2 * ep['zstep'],
                  obj.ellipse_volume[obj.start] * ep['scale'] ** 3,
                  obj.pixel_thresh_vol[obj.start] * ep['scale'] ** 2 * ep['zstep'],
                  obj.ellipse_volume[-1] * ep['scale'] ** 3, obj.pixel_thresh_vol[-1] * ep['scale'] ** 2 * ep['zstep'],
                  obj.vbud[-1] * ep['scale'] ** 3, obj.segment_vol_bud[-1] * ep['scale'] ** 2 * ep['zstep'],
                  obj.zproj_fl[0], obj.pixel_thresh_fluor_vals[0], obj.zproj_fl[obj.start],
                  obj.pixel_thresh_fluor_vals[obj.start], obj.zproj_fl[-1], obj.pixel_thresh_fluor_vals[-1],
                  obj.segment_fl_bud[-1], obj.zproj_fl_bud[-1]] for obj in cc]
df = pd.DataFrame(columns=cols)
for i in range(len(data_vals)):
    df.loc[i] = data_vals[i]
df.to_pickle('./expt_ids'+expt_id+'.pkl')
