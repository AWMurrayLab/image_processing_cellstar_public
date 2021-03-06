import numpy as np
import matplotlib
import pandas as pd
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import os
import scipy
from scipy import stats

pixel_size = {'60X': 0.267, '100X': 0.16}
drange = 65535.0

expt_ids = ['/190606_yFB78_60X_Raff_125uMGal', '/190607_yFB78_60X_Raff_125uMGal', '/190725_yFB78_60X_2Raff_125uMGal',
            '/181207_yFB79_60X_Raff_125uMGal', '/190417_yFB79_60X_Raff_125uMGal', '/190612_yFB79_timelapse']

for expt_id in expt_ids:
    pickle_in = open("./expt_ids"+expt_id+'.pickle',"rb")
    ep = pickle.load(pickle_in)
    print ep
    directory = ep['base_path'] + ep['expt_path'] + '/plots'
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(ep['base_path'] + ep['expt_path']+'/cell_cycles_filtered.pkl', 'rb') as input:
        cc = pickle.load(input)
    temp_cols = ['type', '$t_{G1}$', '$t_{budded}$', '$V_{b,ell}$', '$V_{b,seg}$', '$V_{s,ell}$', '$V_{s,seg}$',
            '$V_{d,ell}$', '$V_{d,seg}$', '$V_{bud,ell}$', '$V_{bud,seg}$',
            '$F1_{b,zproj}$', '$F1_{b,seg}$', '$F1_{s,zproj}$', '$F1_{s,seg}$', '$F1_{d,zproj}$', '$F1_{d,seg}$', '$F1_{bud,zproj}$',
            '$F1_{bud,seg}$', '$F2_{b,zproj}$', '$F2_{b,seg}$', '$F2_{s,zproj}$', '$F2_{s,seg}$', '$F2_{d,zproj}$', '$F2_{d,seg}$',
                 '$F2_{bud,zproj}$', '$F2_{bud,seg}$', '$F1_{b,nuc,int}$', '$c1_{b,nuc,int}$', '$F2_{b,nuc,int}$',
                 '$c2_{b,nuc,int}$', '$V_{b,nuc}$']
    cols = []
    celltype = ['M', 'D', 'Bud', 'NA']


    def nucl_vol(temp):
        if not(temp is None):
            out = len(temp) * ep['scale'] ** 2 * ep['zstep']
            # print temp
            # print ep['scale'] ** 2 * ep['zstep']
            # exit()
        else:
            out = np.nan
        return out
        #

    # 0 indexes the start of the cell cycle, obj.start is the index of Start, -1 is the end of the cell cycle.

    #
    cols+=['type']
    data_vals = [[celltype[obj.celltype]] for obj in cc]
    #
    cols+=['$t_{G1}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0]+=[obj.start * ep['tstep']]
    #
    cols+=['$t_{b}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0]+=[obj.tb * ep['tstep']]
    #
    cols += ['$t_{budded}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [(len(obj.frames) - obj.start) * ep['tstep']]
    #
    cols += ['$V_{b,ell}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [1.33*obj.vb * ep['scale'] ** 3]
    #
    cols += ['$V_{b,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.pixel_thresh_vol[0] * ep['scale'] ** 2 * ep['zstep']]
    #
    cols += ['$V_{s,ell}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [1.33*obj.ellipse_volume[obj.start] * ep['scale'] ** 3]
    #
    cols += ['$V_{s,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.pixel_thresh_vol[obj.start] * ep['scale'] ** 2 * ep['zstep']]
    #
    cols += ['$V_{d,ell}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [1.33*obj.ellipse_volume[-1] * ep['scale'] ** 3]
    #
    cols += ['$V_{d,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.segment_vol_bud[-1] * ep['scale'] ** 2 * ep['zstep']]
    #
    cols += ['$V_{bud,ell}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [1.33*obj.vbud[-1] * ep['scale'] ** 3]
    #
    cols += ['$V_{bud,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.segment_vol_bud[-1] * ep['scale'] ** 2 * ep['zstep']]
    #
    cols += ['$F1_{b,zproj}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.zproj_fl[0]]
    #
    cols += ['$F1_{b,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.pixel_thresh_fluor_vals[0]]
    #
    cols += ['$F1_{s,zproj}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.zproj_fl[obj.start]]
    #
    cols += ['$F1_{s,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.pixel_thresh_fluor_vals[obj.start]]
    #
    cols += ['$F1_{d,zproj}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.zproj_fl[-1]]
    #
    cols += ['$F1_{d,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.pixel_thresh_fluor_vals[-1]]
    #
    cols += ['$F1_{bud,zproj}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.zproj_fl_bud[-1]]
    #
    cols += ['$F1_{bud,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.segment_fl_bud[-1]]
    #
    cols += ['$F2_{b,zproj}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.zproj_fl_c2[0]]
    #
    cols += ['$F2_{b,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.pixel_thresh_fluor_vals_c2[0]]
    #
    cols += ['$F2_{s,zproj}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.zproj_fl_c2[obj.start]]
    #
    cols += ['$F2_{s,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.pixel_thresh_fluor_vals_c2[obj.start]]
    #
    cols += ['$F2_{d,zproj}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.zproj_fl_c2[-1]]
    #
    cols += ['$F2_{d,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.pixel_thresh_fluor_vals_c2[-1]]
    #
    cols += ['$F2_{bud,zproj}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.zproj_fl_bud_c2[-1]]
    #
    cols += ['$F2_{bud,seg}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.segment_fl_bud_c2[-1]]
    #
    cols += ['$F1_{b,nuc}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.nuclear_fluor_int[0]]
    #
    cols += ['$c1_{b,nuc}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.nuclear_fluor_av[0]]
    #
    cols += ['$F2_{b,nuc}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.nuclear_fluor_int_c2[0]]
    #
    cols += ['$c2_{b,nuc}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [obj.nuclear_fluor_av_c2[0]]
    #
    cols += ['$V_{b,nuc}$']
    for i0 in range(len(cc)):
        obj = cc[i0]
        data_vals[i0] += [nucl_vol(obj.nuclear_coords[0])]

    # #
    # cols += []
    # for i0 in range(len(cc)):
    #     obj = cc[i0]
    #     data_vals[i0] += []
    # #
    # cols += []
    # for i0 in range(len(cc)):
    #     obj = cc[i0]
    #     data_vals[i0] += []
    #
    df = pd.DataFrame(columns=cols)
    # print len(cols), len(data_vals[0])
    for i in range(len(data_vals)):
        # print i, data_vals[i]
        df.loc[i] = data_vals[i]

    # now we add values for the fitted slopes of dV/dt, dF1/dt, dF2/dt
    temp_name = '$dV_{fitted,ell}/dt$'
    yv = [(np.asarray(obj.ellipse_volume) + np.asarray(obj.vbud)) * ep['scale'] ** 3 for obj in cc]  # ell vol
    xv = [(np.asarray(range(len(obj.frames))) - obj.start) * ep['tstep'] for obj in cc]
    temp = np.zeros(len(cc))
    for ind in range(len(cc)):
        vals = scipy.stats.linregress(xv[ind], yv[ind])
        temp[ind] = vals[0]  # adding the linear regression
    df[temp_name] = temp

    temp_name = '$dV_{fitted,seg}/dt$'
    yv = [(np.asarray(obj.pixel_thresh_vol) + np.asarray(obj.segment_vol_bud)) * ep['scale'] ** 3 for obj in
          cc]  # seg vol
    temp = np.zeros(len(cc))
    for ind in range(len(cc)):
        vals = scipy.stats.linregress(xv[ind], yv[ind])
        temp[ind] = vals[0]  # adding the linear regression
    df[temp_name] = temp

    temp_name = '$dF1_{fitted,seg}/dt$'
    yv = [(np.asarray(obj.pixel_thresh_fluor_vals) + np.asarray(obj.segment_fl_bud)) for obj in cc]  # F1 seg
    temp = np.zeros(len(cc))
    for ind in range(len(cc)):
        vals = scipy.stats.linregress(xv[ind], yv[ind])
        temp[ind] = vals[0]  # adding the linear regression
    df[temp_name] = temp

    temp_name = '$dF1_{fitted,zproj}/dt$'
    yv = [(np.asarray(obj.zproj_fl) + np.asarray(obj.zproj_fl_bud)) for obj in cc]  # F1 zproj
    temp = np.zeros(len(cc))
    for ind in range(len(cc)):
        vals = scipy.stats.linregress(xv[ind], yv[ind])
        temp[ind] = vals[0]  # adding the linear regression
    df[temp_name] = temp

    temp_name = '$dF2_{fitted,seg}/dt$'
    yv = [(np.asarray(obj.pixel_thresh_fluor_vals_c2) + np.asarray(obj.segment_fl_bud_c2)) for obj in cc]  # F2 seg
    temp = np.zeros(len(cc))
    for ind in range(len(cc)):
        vals = scipy.stats.linregress(xv[ind], yv[ind])
        temp[ind] = vals[0]  # adding the linear regression
    df[temp_name] = temp

    temp_name = '$dF2_{fitted,zproj}/dt$'
    yv = [(np.asarray(obj.zproj_fl_c2) + np.asarray(obj.zproj_fl_bud_c2)) for obj in cc]  # F2 zproj
    temp = np.zeros(len(cc))
    for ind in range(len(cc)):
        vals = scipy.stats.linregress(xv[ind], yv[ind])
        temp[ind] = vals[0]  # adding the linear regression
    df[temp_name] = temp

    # now we add values for the fitted slopes of dV/dt, dF1/dt, dF2/dt
    temp_name = '$d ln(V_{ell})/dt$'
    yv = [(np.asarray(obj.ellipse_volume) + np.asarray(obj.vbud)) * ep['scale'] ** 3 for obj in cc]  # ell vol
    xv = [(np.asarray(range(len(obj.frames))) - obj.start) * ep['tstep'] for obj in cc]
    temp = np.zeros(len(cc))
    for ind in range(len(cc)):
        vals = scipy.stats.linregress(xv[ind], np.log(yv[ind]))
        temp[ind] = vals[0]  # adding the linear regression slope
    df[temp_name] = temp

    temp_name = '$dln(F1_{seg})/dt$'
    yv = [(np.asarray(obj.pixel_thresh_fluor_vals) + np.asarray(obj.segment_fl_bud)) for obj in cc]  # F1 seg
    temp = np.zeros(len(cc))
    for ind in range(len(cc)):
        vals = scipy.stats.linregress(xv[ind], np.log(yv[ind]))
        temp[ind] = vals[0]  # adding the linear regression
    df[temp_name] = temp

    temp_name = '$dln(F2_{seg})/dt$'
    yv = [(np.asarray(obj.pixel_thresh_fluor_vals_c2) + np.asarray(obj.segment_fl_bud_c2)) for obj in cc]  # F2 seg
    temp = np.zeros(len(cc))
    for ind in range(len(cc)):
        vals = scipy.stats.linregress(xv[ind], np.log(yv[ind]))
        temp[ind] = vals[0]  # adding the linear regression
    df[temp_name] = temp

    df.to_pickle('./expt_ids'+expt_id+'.pkl')
