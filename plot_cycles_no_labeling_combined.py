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


# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)

pixel_size = {'60X': 0.267, '100X': 0.16}
drange = 65535.0
# this script takes a list of cell tracks with daughter lineages assigned (e.g. the output from assign_lineages_1.py)
# and creates a collection of reliable cell cycles from this by iterating through this list.
color = cm.tab10(np.linspace(0, 1, 10))
print color[0]

timestep = 10.0
fluor_c2 = True
zscale=0.7
scale = pixel_size['60X']

# yFB79 Raffinose experiment on 181207

base_path = '/scratch/lab/image_analysis_scratch'
expt_paths = [['/181207_yFB79_60X_Raff_125uMGal/timelapse','/190417_yFB79_timelapse/timelapse'],\
                ['/190606_yFB78_timelapse/timelapse', '/190725_yFB78_timelapse/timelapse']]
strains = ['$P_{WHI5}-WHI5$','$P_{GAL1}-WHI5$']
strain_labels = ['PWHI5','PGAL1']

celltypes = ['Mother', 'Daughter']
directory = '/home/felix/Dropbox/19_whi5_dilution_paper/plots/induction_validation'

for ind in range(len(strains)):
    cc = []
    for expt_path in expt_paths[ind]:
        with open(base_path+expt_path+'/cell_cycles_filtered.pkl', 'rb') as input:
            cc_temp = pickle.load(input)
        cc+=cc_temp
    print len(cc)
    fig = plt.figure(figsize=[5, 5])
    filt_cc = [obj for obj in cc if obj.complete and not (obj.error) and
               not (obj.daughter is None) and not np.sum(np.asarray(obj.pixel_thresh_vol)==0)>0 and obj.celltype==1]
    if len(filt_cc) > 0:
        xv = [len(obj.frames[:obj.start]) for obj in filt_cc]
        xv1 = timestep * np.arange(np.amax(xv))
        yv = np.full([len(xv), np.amax(xv)], np.nan)
        for i0 in range(len(filt_cc)):
            yv[i0, :xv[i0]] = 1.0 / np.asarray(filt_cc[i0].pixel_thresh_vol[:filt_cc[i0].start])
        yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
        # print np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
        thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
        print 'threshold', thresh
        # print np.nanmean(yv, axis=0)[:thresh]
        temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                          label=' dilution factor', linewidth=2, alpha=0.8)
        plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh]/np.sqrt(yv.shape[0]-np.sum(np.isnan(yv),axis=0))[:thresh],
                         np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh]/np.sqrt(yv.shape[0]-np.sum(np.isnan(yv),axis=0))[:thresh], alpha=0.2)
        # print np.nanmean(yv, axis=0)[:thresh]
        yv = np.full([len(xv), np.amax(xv)], np.nan)
        for i0 in range(len(filt_cc)):
            yv[i0, :xv[i0]] = np.asarray(filt_cc[i0].zproj_fl_c2[:filt_cc[i0].start]) / np.asarray(
                filt_cc[i0].pixel_thresh_vol[:filt_cc[i0].start])
        yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
        thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
        print 'threshold', thresh
        # print np.nanmean(yv, axis=0)[:thresh]
        temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                          label=r' mCherry', linewidth=2, alpha=0.8)
        plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh]/np.sqrt(yv.shape[0]-np.sum(np.isnan(yv),axis=0))[:thresh],
                         np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh]/np.sqrt(yv.shape[0]-np.sum(np.isnan(yv),axis=0))[:thresh], alpha=0.2)
        yv = np.full([len(xv), np.amax(xv)], np.nan)
        for i0 in range(len(filt_cc)):
            yv[i0, :xv[i0]] = np.asarray(filt_cc[i0].zproj_fl[:filt_cc[i0].start]) / np.asarray(
                filt_cc[i0].pixel_thresh_vol[:filt_cc[i0].start])
        yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
        thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
        print 'threshold', thresh
        # print np.nanmean(yv, axis=0)[:thresh]
        temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                          label=r' Whi5', linewidth=2, alpha=0.8)

        plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh]/np.sqrt(yv.shape[0]-np.sum(np.isnan(yv),axis=0))[:thresh],
                         np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh]/np.sqrt(yv.shape[0]-np.sum(np.isnan(yv),axis=0))[:thresh], alpha=0.2)

    plt.ylabel('Relative dilution')
    plt.xlabel('Time in G1 (mins)')
    plt.legend()
    plt.title(strains[ind]+'; Cell number = {0}'.format(len(xv)))
    fig.savefig(directory + '/G1_dilution_SEM_'+strain_labels[ind]+'.png', bbox_inches='tight', dpi=300)
    plt.clf()


for ind in range(len(strains)):
    cc = []
    for expt_path in expt_paths[ind]:
        with open(base_path+expt_path+'/cell_cycles_filtered.pkl', 'rb') as input:
            cc_temp = pickle.load(input)
        cc+=cc_temp
    print len(cc)
    fig = plt.figure(figsize=[5, 5])
    filt_cc = [obj for obj in cc if obj.complete and not (obj.error) and
               not (obj.daughter is None) and not np.sum(np.asarray(obj.pixel_thresh_vol)==0)>0 and obj.celltype==1]
    if len(filt_cc) > 0:
        xv = [len(obj.frames[:obj.start]) for obj in filt_cc]
        xv1 = timestep * np.arange(np.amax(xv))
        yv = np.full([len(xv), np.amax(xv)], np.nan)
        for i0 in range(len(filt_cc)):
            yv[i0, :xv[i0]] = 1.0 / np.asarray(filt_cc[i0].pixel_thresh_vol[:filt_cc[i0].start])
        yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
        # print np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
        thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
        print 'threshold', thresh
        # print np.nanmean(yv, axis=0)[:thresh]
        temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                          label=' dilution factor', linewidth=2, alpha=0.8)
        plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh],
                         np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh], alpha=0.2)
        # print np.nanmean(yv, axis=0)[:thresh]
        yv = np.full([len(xv), np.amax(xv)], np.nan)
        for i0 in range(len(filt_cc)):
            yv[i0, :xv[i0]] = np.asarray(filt_cc[i0].zproj_fl_c2[:filt_cc[i0].start]) / np.asarray(
                filt_cc[i0].pixel_thresh_vol[:filt_cc[i0].start])
        yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
        thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
        print 'threshold', thresh
        # print np.nanmean(yv, axis=0)[:thresh]
        temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                          label=r' mCherry', linewidth=2, alpha=0.8)
        plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh],
                         np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh], alpha=0.2)
        yv = np.full([len(xv), np.amax(xv)], np.nan)
        for i0 in range(len(filt_cc)):
            yv[i0, :xv[i0]] = np.asarray(filt_cc[i0].zproj_fl[:filt_cc[i0].start]) / np.asarray(
                filt_cc[i0].pixel_thresh_vol[:filt_cc[i0].start])
        yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
        thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
        print 'threshold', thresh
        # print np.nanmean(yv, axis=0)[:thresh]
        temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                          label=r' Whi5', linewidth=2, alpha=0.8)
        plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh],
                         np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh], alpha=0.2)

    plt.ylabel('Relative dilution')
    plt.xlabel('Time in G1 (mins)')
    plt.legend()
    plt.title(strains[ind]+'; Cell number = {0}'.format(len(xv)))
    fig.savefig(directory + '/G1_dilution_Std_dev_'+strain_labels[ind]+'.png', bbox_inches='tight', dpi=300)
    plt.clf()

