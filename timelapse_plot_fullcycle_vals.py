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
# import seaborn as sns
import scipy
from scipy import stats

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

expt_ids = ['/190606_yFB78_60X_Raff_125uMGal', '/190607_yFB78_60X_Raff_125uMGal', '/190725_yFB78_60X_2Raff_125uMGal',
            '/181207_yFB79_60X_Raff_125uMGal', '/190417_yFB79_60X_Raff_125uMGal', '/190612_yFB79_timelapse']

for expt_id in expt_ids:
    pickle_in = open("./expt_ids"+expt_id+'.pickle',"rb")
    ep = pickle.load(pickle_in)
    directory = ep['base_path'] + ep['expt_path'] + '/plots'
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(ep['base_path'] + ep['expt_path']+'/cell_cycles_filtered.pkl', 'rb') as input:
        cc = pickle.load(input)
    celltype = ['M', 'D', 'Bud', 'NA']

    # volume vs. time
    yv = [(np.asarray(obj.ellipse_volume)+np.asarray(obj.vbud)) * ep['scale'] ** 3 for obj in cc]  # full volume
    xv = [(np.asarray(range(len(obj.frames)))-obj.start)*ep['tstep'] for obj in cc]
    inds = np.random.randint(low=0,high=len(yv),size=5)
    fig = plt.figure(figsize=[5, 5])
    for ind in inds:
        plt.plot(xv[ind],yv[ind],label='ind = {0}'.format(ind))
        vals=scipy.stats.linregress(xv[ind],yv[ind])
        plt.plot(xv[ind], vals[0]*xv[ind]+vals[1],label='linregress, ind={0}'.format(ind))
    plt.legend()
    plt.xlabel('Time relative to Start (mins)')
    plt.ylabel('Volume (fL)')
    fig.savefig('./plots/testing_cycle_vals/'+expt_id[1:] + '_Vol_time.png',bbox_inches='tight',dpi=fig.dpi)
    plt.clf()

    # volume pixel segmented vs. time
    yv = [(np.asarray(obj.pixel_thresh_vol)+np.asarray(obj.segment_vol_bud)) * ep['scale'] ** 3 for obj in cc]  # full volume
    xv = [(np.asarray(range(len(obj.frames)))-obj.start)*ep['tstep'] for obj in cc]
    fig = plt.figure(figsize=[5, 5])
    for ind in inds:
        plt.plot(xv[ind],yv[ind],label='ind = {0}'.format(ind))
        vals=scipy.stats.linregress(xv[ind],yv[ind])
        plt.plot(xv[ind], vals[0]*xv[ind]+vals[1],label='linregress, ind={0}'.format(ind))
    plt.legend()
    plt.xlabel('Time relative to Start (mins)')
    plt.ylabel('Volume (fL)')
    fig.savefig('./plots/testing_cycle_vals/'+expt_id[1:] + '_Vol_seg_time.png',bbox_inches='tight',dpi=fig.dpi)
    plt.clf()

    # C1 integrated fluorescence vs. time
    yv = [(np.asarray(obj.pixel_thresh_fluor_vals) + np.asarray(obj.segment_fl_bud)) for obj in cc]  # full volume
    fig = plt.figure(figsize=[5, 5])
    for ind in inds:
        plt.plot(xv[ind], yv[ind], label='ind = {0}'.format(ind))
        vals = scipy.stats.linregress(xv[ind], yv[ind])
        plt.plot(xv[ind], vals[0] * xv[ind] + vals[1], label='linregress, ind={0}'.format(ind))
    plt.legend()
    plt.xlabel('Time relative to Start (mins)')
    plt.ylabel('Integrated fluorescence (A.U.)')
    fig.savefig('./plots/testing_cycle_vals/'+expt_id[1:] + '_F1_time.png',bbox_inches='tight',dpi=fig.dpi)
    plt.clf()

    # C2 integrated fluorescence vs. time
    yv = [(np.asarray(obj.pixel_thresh_fluor_vals_c2) + np.asarray(obj.segment_fl_bud_c2)) for obj in cc]  # full volume
    fig = plt.figure(figsize=[5, 5])
    for ind in inds:
        plt.plot(xv[ind], yv[ind], label='ind = {0}'.format(ind))
        vals = scipy.stats.linregress(xv[ind], yv[ind])
        plt.plot(xv[ind], vals[0] * xv[ind] + vals[1], label='linregress, ind={0}'.format(ind))
    plt.legend()
    plt.xlabel('Time relative to Start (mins)')
    plt.ylabel('Integrated fluorescence (A.U.)')
    fig.savefig('./plots/testing_cycle_vals/'+expt_id[1:] + '_F2_time.png',bbox_inches='tight',dpi=fig.dpi)
    plt.clf()

    # C1 zproj fluorescence vs. time
    yv = [(np.asarray(obj.zproj_fl) + np.asarray(obj.zproj_fl_bud)) for obj in cc]  # full volume
    fig = plt.figure(figsize=[5, 5])
    for ind in inds:
        plt.plot(xv[ind], yv[ind], label='ind = {0}'.format(ind))
        vals = scipy.stats.linregress(xv[ind], yv[ind])
        plt.plot(xv[ind], vals[0] * xv[ind] + vals[1], label='linregress, ind={0}'.format(ind))
    plt.legend()
    plt.xlabel('Time relative to Start (mins)')
    plt.ylabel('Integrated fluorescence (A.U.)')
    fig.savefig('./plots/testing_cycle_vals/'+expt_id[1:] + '_F1_zproj_time.png',bbox_inches='tight',dpi=fig.dpi)
    plt.clf()

    # C2 integrated fluorescence vs. time
    yv = [(np.asarray(obj.zproj_fl_c2) + np.asarray(obj.zproj_fl_bud_c2)) for obj in cc]  # full volume
    fig = plt.figure(figsize=[5, 5])
    for ind in inds:
        plt.plot(xv[ind], yv[ind], label='ind = {0}'.format(ind))
        vals = scipy.stats.linregress(xv[ind], yv[ind])
        plt.plot(xv[ind], vals[0] * xv[ind] + vals[1], label='linregress, ind={0}'.format(ind))
    plt.legend()
    plt.xlabel('Time relative to Start (mins)')
    plt.ylabel('Integrated fluorescence (A.U.)')
    fig.savefig('./plots/testing_cycle_vals/'+expt_id[1:] + '_F2_zproj_time.png',bbox_inches='tight',dpi=fig.dpi)
    plt.clf()

    # Vol integrated fluorescence vs. time as a logarithm
    yv = [(np.asarray(obj.ellipse_volume) + np.asarray(obj.vbud)) * ep['scale'] ** 3 for obj in cc]  # full volume
    fig = plt.figure(figsize=[5, 5])
    for ind in inds:
        plt.plot(xv[ind], yv[ind], label='ind = {0}'.format(ind))
        vals = scipy.stats.linregress(xv[ind], np.log(yv[ind]))
        plt.plot(xv[ind], np.exp(vals[0] * xv[ind] + vals[1]), label='linregress, ind={0}'.format(ind))
    plt.legend()
    plt.xlabel('Time relative to Start (mins)')
    plt.ylabel('Integrated fluorescence (A.U.)')
    fig.savefig('./plots/testing_cycle_vals/'+expt_id[1:] + '_V_ell_log_time.png',bbox_inches='tight',dpi=fig.dpi)
    plt.clf()

    # C1 integrated fluorescence vs. time as a logarithm
    yv = [(np.asarray(obj.pixel_thresh_fluor_vals) + np.asarray(obj.segment_fl_bud)) for obj in cc]  # full volume
    fig = plt.figure(figsize=[5, 5])
    for ind in inds:
        plt.plot(xv[ind], yv[ind], label='ind = {0}'.format(ind))
        vals = scipy.stats.linregress(xv[ind], np.log(yv[ind]))
        plt.plot(xv[ind], np.exp(vals[0] * xv[ind] + vals[1]), label='linregress, ind={0}'.format(ind))
    plt.legend()
    plt.xlabel('Time relative to Start (mins)')
    plt.ylabel('Integrated fluorescence (A.U.)')
    fig.savefig('./plots/testing_cycle_vals/'+expt_id[1:] + '_F1_log_time.png',bbox_inches='tight',dpi=fig.dpi)
    plt.clf()

    # C2 integrated fluorescence vs. time as a logarithm
    yv = [(np.asarray(obj.pixel_thresh_fluor_vals_c2) + np.asarray(obj.segment_fl_bud_c2)) for obj in cc]  # full volume
    fig = plt.figure(figsize=[5, 5])
    for ind in inds:
        plt.plot(xv[ind], yv[ind], label='ind = {0}'.format(ind))
        vals = scipy.stats.linregress(xv[ind], np.log(yv[ind]))
        plt.plot(xv[ind], np.exp(vals[0] * xv[ind] + vals[1]), label='linregress, ind={0}'.format(ind))
    plt.legend()
    plt.xlabel('Time relative to Start (mins)')
    plt.ylabel('Integrated fluorescence (A.U.)')
    fig.savefig('./plots/testing_cycle_vals/'+expt_id[1:] + '_F2_log_time.png',bbox_inches='tight',dpi=fig.dpi)
    plt.clf()
