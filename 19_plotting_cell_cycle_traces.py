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
print color[0]

# list of experiments that we include in this analysis

# # pACT1-mCherry experiment on 181114, yFB78 expt 190607
# expt_ids = ['/181114_yFB78_Raff_125Gal']
# expt_ids = ['/181114_yFB78_Raff_125Gal', '/190607_yFB78_60X_Raff_125uMGal']  # not using this since the timestep is different
# timestep = 10.0
# genotype = 'pGAL1-WHI5-mVenNB'

# yFB79 Raffinose experiment on 181207, yFB79 expt 190417
expt_ids = ['/181207_yFB79_60X_Raff_125uMGal', '/181917_yFB79_60X_Raff_125uMGal']
timestep = 10.0
genotype = 'pWHI5-WHI5-mVenNB'

eps = []
cc = []
for i0 in range(len(expt_ids)):
    pickle_in = open("./expt_ids"+expt_ids[i0]+'.pickle',"rb")
    ep = pickle.load(pickle_in)
    print ep
    with open(ep['base_path'] + ep['expt_path']+'/cell_cycles_filtered.pkl', 'rb') as input:
        temp_cc = pickle.load(input)
        for obj in temp_cc:
            obj.expt_id = i0  # ensuring that we can use the correct settings for each experiment.
        cc+=temp_cc
    eps.append(ep)

# Fluor linear slope vs vb from birth to start
# Fluor pointwise slope from Start to division
# Concentration vs volume at division
# Thresholded volume growth rate

# initial definitions
celltypes = ['Mother', 'Daughter']

filt_cc = [obj for obj in cc if obj.complete and not (obj.error) and
           not (obj.daughter is None) and obj.celltype==1]
fig=plt.figure(figsize=[5,5])
if len(filt_cc) > 0:
    xv = [len(obj.frames[:obj.start]) for obj in filt_cc]
    xv1 = timestep * np.arange(np.amax(xv))
    yv = np.full([len(xv), np.amax(xv)], np.nan)
    for i0 in range(len(filt_cc)):
        yv[i0, :xv[i0]] = 1.0 / np.asarray(filt_cc[i0].pixel_thresh_vol[:filt_cc[i0].start])
    yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
    thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
    print 'threshold', thresh
    temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                      label=' dilution factor. Cell num={0}'.format(len(xv)), linewidth=2, alpha=0.8)
    plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh],
                     np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh], alpha=0.2)
    yv = np.full([len(xv), np.amax(xv)], np.nan)
    for i0 in range(len(filt_cc)):
        yv[i0, :xv[i0]] = np.asarray(filt_cc[i0].zproj_fl_c2[:filt_cc[i0].start]) / np.asarray(
            filt_cc[i0].pixel_thresh_vol[:filt_cc[i0].start])
    yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
    thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
    print 'threshold', thresh
    temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                      label=' pACT1-mCherry', linewidth=2, alpha=0.8)
    plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh],
                     np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh], alpha=0.2)
    for i0 in range(len(filt_cc)):
        yv[i0, :xv[i0]] = np.asarray(filt_cc[i0].zproj_fl[:filt_cc[i0].start]) / np.asarray(
            filt_cc[i0].pixel_thresh_vol[:filt_cc[i0].start])
    yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
    thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
    print 'threshold', thresh
    temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                      label=genotype, linewidth=2, alpha=0.8)
    plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh],
                     np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh], alpha=0.2)

plt.ylim(ymin=0)
plt.ylabel('Relative dilution')
plt.xlabel('Time post cell birth')
plt.legend()
plt.title('Dilution factor due to growth in G1')
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/microscopy_validation/G1_dilution_timecourse_'+genotype[:5]+'.png',
            bbox_inches='tight', dpi=fig.dpi)
plt.clf()

fig = plt.figure(figsize=[5, 5])
filt_cc = [obj for obj in cc if obj.complete and not (obj.error) and
           not (obj.daughter is None) and obj.celltype==1]
if len(filt_cc) > 0:
    xv = [len(obj.frames) for obj in filt_cc]
    xv1 = timestep * np.arange(np.amax(xv))
    yv = np.full([len(xv), np.amax(xv)], np.nan)
    for i0 in range(len(filt_cc)):
        yv[i0, :xv[i0]] = 1.0 / np.asarray(filt_cc[i0].pixel_thresh_vol)
    yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
    thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
    print 'threshold', thresh
    temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                      label=' dilution factor. Cell num={0}'.format(len(xv)), linewidth=2, alpha=0.8)
    plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh],
                     np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh], alpha=0.2)

    yv = np.full([len(xv), np.amax(xv)], np.nan)
    for i0 in range(len(filt_cc)):
        yv[i0, :xv[i0]] = (np.asarray(filt_cc[i0].zproj_fl_c2)+np.array(filt_cc[i0].zproj_fl_bud_c2)) / np.asarray(
            filt_cc[i0].pixel_thresh_vol)
    yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
    thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
    print 'threshold', thresh
    temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                      label=' pACT1-mCherry', linewidth=2, alpha=0.8)
    plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh],
                     np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh], alpha=0.2)
    for i0 in range(len(filt_cc)):
        yv[i0, :xv[i0]] = (np.asarray(filt_cc[i0].zproj_fl)+np.array(filt_cc[i0].zproj_fl_bud)) / np.asarray(
            filt_cc[i0].pixel_thresh_vol)
    yv /= np.repeat(np.reshape(yv[:, 0], [yv.shape[0], 1]), yv.shape[1], axis=1)
    thresh = np.nonzero(np.sum(~np.isnan(yv), axis=0) < 5)[0][0]
    print 'threshold', thresh
    temp_p = plt.plot(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh],
                      label=genotype, linewidth=2, alpha=0.8)
    plt.fill_between(xv1[:thresh], np.nanmean(yv, axis=0)[:thresh] - np.nanstd(yv, axis=0)[:thresh],
                     np.nanmean(yv, axis=0)[:thresh] + np.nanstd(yv, axis=0)[:thresh], alpha=0.2)

plt.ylabel('Relative dilution')
plt.xlabel('Time post cell birth')
plt.legend()
plt.title('Dilution factor due to growth')
fig.savefig('/home/felix/Dropbox/19_whi5_dilution_paper/plots/microscopy_validation/full_cycle_dilution_timecourse_'+genotype[:5]+'.png',
            bbox_inches='tight', dpi=fig.dpi)
plt.clf()
