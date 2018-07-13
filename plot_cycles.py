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

timestep=10.0
scale = 0.267
# this script takes a list of cell tracks with daughter lineages assigned (e.g. the output from assign_lineages_1.py)
# and creates a collection of reliable cell cycles from this by iterating through this list.
color=cm.tab10(np.linspace(0,1,10))
print color[0]
base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
directory = base_path + expt_path + '/plots'
if not os.path.exists(directory):
    os.makedirs(directory)
with open(base_path+expt_path+'/cell_cycles_compiled.pkl', 'rb') as input:
    cc = pickle.load(input)
# initial definitions
fluor_name = '/180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w2515 laser 20_'
drange = 65535.0
celltypes = ['Mother', 'Daughter']
pops = ['Pwhi5', 'Pgal']  # label = Pgal, unlabel = Pwhi5

fig = plt.figure(figsize=[5, 5])
for ind in range(2):
    filt_cc = [obj for obj in cc if obj.label_type == ind and obj.complete and not(obj.error) and
               not(obj.daughter is None)]
    print len(filt_cc)
    if len(filt_cc)>0:
        xv = [obj.ellipse_volume[obj.start] for obj in filt_cc]
        yv=[]
        for obj in filt_cc:
            # temp = obj.int_fl[1]+obj.
            yv1 = [np.asarray(obj.int_fl[obj.start :])+np.asarray(obj.int_fl_bud[obj.start :])]
            xvtemp = np.asarray(range(obj.start, len(obj.frames)))*timestep
            temp1 = scipy.stats.linregress(xvtemp, yv1)
            yv.append(temp1[0])
    # if len(filt_cc) > 0:
    temp = scipy.stats.pearsonr(xv, yv)
    plt.plot(xv, yv, marker='.', label=pops[ind]+' PCC={0}, num cells={1}'.format(np.round(temp[0], 3), len(filt_cc)), linestyle='None')
plt.xlabel('$V_b$ ($\mu m^3$)')
plt.ylabel('$\Delta F$')
plt.legend()
plt.title('Fluorescence rate of change')
fig.savefig(directory+'/df_vb.png', bbox_inches='tight', dpi=fig.dpi)

del fig
fig = plt.figure(figsize=[5, 5])
for ind in range(2):
    filt_cc = [obj for obj in cc if obj.label_type == ind and obj.complete and not(obj.error) and
               not(obj.daughter is None)]
    print len(filt_cc)
    if len(filt_cc)>0:
        xv = [obj.ellipse_volume[obj.start] for obj in filt_cc]
        yv=[]
        for obj in filt_cc:
            # temp = obj.int_fl[1]+obj.
            yv1 = [np.asarray(obj.zproj_fl[obj.start :])+np.asarray(obj.zproj_fl_bud[obj.start :])]
            xvtemp = np.asarray(range(obj.start, len(obj.frames)))*timestep
            temp1 = scipy.stats.linregress(xvtemp, yv1)
            yv.append(temp1[0])
    # if len(filt_cc) > 0:
    temp = scipy.stats.pearsonr(xv, yv)
    plt.plot(xv, yv, marker='.', label=pops[ind]+' Z sum Fluor. PCC={0}, num cells={1}'.format(np.round(temp[0], 3), len(filt_cc)), linestyle='None')
plt.xlabel('$V_b$ ($\mu m^3$)')
plt.ylabel('$\Delta F$')
plt.legend()
plt.title('Fluorescence rate of change')
fig.savefig(directory+'/df_zproj_vb.png', bbox_inches='tight', dpi=fig.dpi)
del fig


fig = plt.figure(figsize=[5, 5])
# for ind in range(2):
filt_cc = [obj for obj in cc if obj.complete and not(obj.error) and not(obj.daughter is None)]
print len(filt_cc)
if len(filt_cc)>0:
    xv = []
    yv = []
    for obj in filt_cc:
        # temp = obj.int_fl[1]+obj.
        xv1 = [obj.int_fl[ind] + obj.int_fl_bud[ind] for ind in range(len(obj.frames))]
        yv1 = [obj.zproj_fl[ind] + obj.zproj_fl_bud[ind] for ind in range(len(obj.frames))]
        xv+=xv1
        yv+=yv1
    # if len(filt_cc) > 0:
    temp = scipy.stats.pearsonr(xv, yv)
    plt.plot(xv, yv, marker='.', label='Fluor corr. PCC={0}'.format(np.round(temp[0], 3)), linestyle='None')
    xvals = np.linspace(np.mean(xv) - 2.5 * np.std(xv), np.mean(xv) + 5.0 * np.std(xv), 20)
    vals = scipy.stats.linregress(xv, yv)
    plt.plot(xvals, xvals, label='X=Y')
plt.xlabel('Volume integrated Fluorescence')
plt.ylabel('Z projected fluorescence')
plt.legend()
plt.title('Validation of fluorescence measurements')
fig.savefig(directory+'/fluor_validation.png', bbox_inches='tight', dpi=fig.dpi)


for ind in range(2):
    print celltypes[ind]
    filt_cc = [obj for obj in cc if obj.celltype == ind and obj.complete and not(obj.daughter is None) and not obj.error]
    filt_cc1 = [obj for obj in filt_cc if obj.label_type==1]
    filt_cc2 = [obj for obj in filt_cc if obj.label_type == 0]
    print 'length', len(filt_cc)
    if len(filt_cc)>=1:
        # Plots of single cell volume traces
        fig = plt.figure(figsize=[5, 5])
        # temp = []
        # amin = 0
        # amax = 0
        for obj in filt_cc:
            yv = np.asarray([(obj.ellipse_volume[i0]) * scale ** 3 for i0 in range(len(obj.frames))])
            xv = range(-obj.start, len(obj.frames)-obj.start)
            # amin = np.minimum(amin, -obj.start)
            # amax = np.max(amax, len(obj.frames)-obj.start)
        # xv = []\
        #     (amax-amin)
        # for obj in filt_cc:
        #     for ind1 in range(len(obj.frames)):
        #         if ind1 >= len(temp):
        #             temp.append([])
        #         temp[ind1].append(((obj.ellipse_volume[ind1]) / obj.ellipse_volume[0]))
        #     xmin = np.maximum(0.0, np.mean(xv) - 2.5 * np.std(xv))
        #     xvals = np.linspace(xmin, np.mean(xv) + 2.5 * np.std(xv), 20)
        #     plt.xlim(xmin=xmin)
            plt.plot(xv, yv / yv[0], marker=None, alpha=0.2, linewidth=1)
        # yv = [np.mean(vals) for vals in temp]
        # xv = np.asarray(range(len(yv))) * timestep
        # plt.plot(xv, yv, marker=None, label='Mean', linewidth=2, color='b', alpha=0.5)
        plt.xlabel('$t$ (minutes relative to Start)')
        plt.ylabel('$V(t)$ $(\mu m^3)$')
        plt.title(celltypes[ind] + r' Volume traces, cellnum = {0}'.format(len(filt_cc)))
        # temp = scipy.stats.linregress(xv, yv)
        # plt.plot(xvals, temp[0]*xvals+temp[1], label='slope = {0}'.format(str(np.round(temp[0],2))))
        # plt.legend()
        # print temp
        directory = base_path + expt_path + '/plots'
        if not os.path.exists(directory):
            os.makedirs(directory)
        fig.savefig(directory + '/' + celltypes[ind] + '_main_cell_volume_traces.png', bbox_inches='tight', dpi=fig.dpi)


        # Plots of single cell volume traces
        fig = plt.figure(figsize=[5, 5])
        temp = []
        for obj in filt_cc:
            xv, yv = np.asarray(range(len(obj.frames)))*timestep,\
                     np.asarray([(obj.ellipse_volume[i0]+obj.vbud[i0])*scale**3 for i0 in range(len(obj.frames))]),
            for ind1 in range(len(obj.frames)):
                if ind1 >= len(temp):
                    temp.append([])
                temp[ind1].append(((obj.ellipse_volume[ind1]+obj.vbud[ind1])/obj.ellipse_volume[0]))
            xmin = np.maximum(0.0, np.mean(xv)-2.5*np.std(xv))
            xvals = np.linspace(xmin, np.mean(xv)+2.5*np.std(xv), 20)
            # plt.xlim(xmin=xmin)
            plt.plot(xv, yv/yv[0], marker=None, alpha=0.2, linewidth=1)
        yv = [np.mean(vals) for vals in temp]
        xv = np.asarray(range(len(yv)))*timestep
        plt.plot(xv, yv, marker=None, label='Mean', linewidth=2, color='b', alpha=0.5)
        plt.xlabel('$t$ (minutes)')
        plt.ylabel('$V(t)$ $(\mu m^3)$')
        plt.title(celltypes[ind]+r' Volume traces, cellnum = {0}'.format(len(filt_cc)))
        # temp = scipy.stats.linregress(xv, yv)
        # plt.plot(xvals, temp[0]*xvals+temp[1], label='slope = {0}'.format(str(np.round(temp[0],2))))
        # plt.legend()
        # print temp
        directory = base_path+expt_path+'/plots'
        if not os.path.exists(directory):
            os.makedirs(directory)
        fig.savefig(directory+'/'+celltypes[ind]+'_volume_traces.png', bbox_inches='tight', dpi=fig.dpi)

        # Plots of single cell fluorescence traces
        fig = plt.figure(figsize=[5, 5])
        temp = []
        for obj in filt_cc:
            xv, yv = np.asarray(range(len(obj.frames)))*timestep,\
                     np.asarray([(obj.int_fl[i0]+obj.int_fl_bud[i0])*scale**3 for i0 in range(len(obj.frames))]),
            for ind1 in range(len(obj.frames)):
                if ind1 >= len(temp):
                    temp.append([])
                temp[ind1].append(((obj.int_fl[ind1]+obj.int_fl_bud[ind1])/obj.int_fl[0]))
            xmin = np.maximum(0.0, np.mean(xv)-2.5*np.std(xv))
            xvals = np.linspace(xmin, np.mean(xv)+2.5*np.std(xv), 20)
            # plt.xlim(xmin=xmin)
            plt.plot(xv, yv/yv[0], marker=None, alpha=0.2, linewidth=1)
        yv = [np.mean(vals) for vals in temp]
        xv = np.asarray(range(len(yv)))*timestep
        plt.plot(xv, yv, marker=None, label='Mean', linewidth=2, color='b', alpha=0.5)
        plt.xlabel('$t$ (minutes)')
        plt.ylabel('$F(t)$ (A.U.)')
        plt.title(celltypes[ind]+r' Fluorescence traces, cellnum = {0}'.format(len(filt_cc)))
        # temp = scipy.stats.linregress(xv, yv)
        # plt.plot(xvals, temp[0]*xvals+temp[1], label='slope = {0}'.format(str(np.round(temp[0],2))))
        # plt.legend()
        # print temp
        directory = base_path+expt_path+'/plots'
        if not os.path.exists(directory):
            os.makedirs(directory)
        fig.savefig(directory+'/'+celltypes[ind]+'_fluorescence_traces.png', bbox_inches='tight', dpi=fig.dpi)

        # print 'daughters exist!'
        fig = plt.figure(figsize=[5, 5])

        xv, yv = np.asarray([obj.vb for obj in filt_cc1])*scale**3, np.asarray([obj.vd+obj.vbud[-1] for obj in filt_cc1])*scale**3
        xmin = np.maximum(0.0, np.mean(xv)-2.5*np.std(xv))
        xvals = np.linspace(xmin, np.mean(xv)+2.5*np.std(xv), 20)
        # plt.xlim(xmin=xmin)
        plt.plot(xv, yv, '.', label=celltypes[ind]+' dyed data, cellnum={0}'.format(len(filt_cc1)), linestyle='None', color='b')
        temp = scipy.stats.linregress(xv, yv)
        plt.plot(xvals, temp[0] * xvals + temp[1], label='Dyed slope = {0}'.format(str(np.round(temp[0], 2))))

        xv, yv = np.asarray([obj.vb for obj in filt_cc2]) * scale ** 3, np.asarray(
            [obj.vd + obj.vbud[-1] for obj in filt_cc2]) * scale ** 3
        xmin = np.maximum(0.0, np.mean(xv) - 2.5 * np.std(xv))
        xvals = np.linspace(xmin, np.mean(xv) + 2.5 * np.std(xv), 20)
        # plt.xlim(xmin=xmin)
        plt.plot(xv, yv, '.', label=celltypes[ind] + ' undyed data, cellnum={0}'.format(len(filt_cc2)), linestyle='None', color='r')
        temp = scipy.stats.linregress(xv, yv)
        plt.plot(xvals, temp[0] * xvals + temp[1], label='Undyed slope = {0}'.format(str(np.round(temp[0], 2))))

        plt.xlabel('$V_b$ $(\mu m^3)$')
        plt.ylabel('$V_d$ $(\mu m^3)$')
        plt.title(celltypes[ind]+r' Volume correlations')

        plt.legend()
        # print temp
        directory = base_path+expt_path+'/plots'
        if not os.path.exists(directory):
            os.makedirs(directory)
        fig.savefig(directory+'/'+celltypes[ind]+'_vb_vd.png', bbox_inches='tight', dpi=fig.dpi)


        fig = plt.figure(figsize=[5, 5])
        xv, yv = np.asarray([obj.vb for obj in filt_cc1])*scale**3, np.asarray([len(obj.frames)*timestep for obj in filt_cc1])
        xmin = max(0.0, np.mean(xv) - 2.5 * np.std(xv))
        xvals = np.linspace(xmin, np.mean(xv) + 2.5 * np.std(xv), 20)
        temp = scipy.stats.pearsonr(xv, yv)
        plt.plot(xv, yv, '.', label=celltypes[ind]+' dyed data, PCC ={0}, Cellnum={1}'.format(np.round(temp[0], 3), len(filt_cc1)), linestyle='None')
        xv, yv = np.asarray([obj.vb for obj in filt_cc2]) * scale ** 3, np.asarray(
            [len(obj.frames) * timestep for obj in filt_cc2])
        xmin = max(0.0, np.mean(xv) - 2.5 * np.std(xv))
        xvals = np.linspace(xmin, np.mean(xv) + 2.5 * np.std(xv), 20)
        temp = scipy.stats.pearsonr(xv, yv)
        plt.plot(xv, yv, '.', label=celltypes[ind] + ' undyed data, PCC ={0}, Cellnum={1}'.format(np.round(temp[0], 3), len(filt_cc2)),
                 linestyle='None')
        plt.xlabel('$V_b$ $(\mu m^3)$')
        plt.ylabel('$t_{div}$ (minutes)')
        plt.title(celltypes[ind]+r' Division timing Correlations')
        # temp = scipy.stats.linregress(xv, yv)
        # plt.plot(xvals, temp[0]*xvals+temp[1], label='slope = {0}'.format(str(np.round(temp[0],2))))
        plt.legend()
        # print temp
        directory = base_path+expt_path+'/plots'
        if not os.path.exists(directory):
            os.makedirs(directory)
        fig.savefig(directory+'/'+celltypes[ind]+'_vb_tdiv.png', bbox_inches='tight', dpi=fig.dpi)

        fig = plt.figure(figsize=[5, 5])
        xv, yv = np.asarray([obj.vb for obj in filt_cc1])*scale**3, np.asarray([obj.start*timestep for obj in filt_cc1])
        xmin = max(0.0, np.mean(xv) - 2.5 * np.std(xv))
        xvals = np.linspace(xmin, np.mean(xv) + 2.5 * np.std(xv), 20)
        temp = scipy.stats.pearsonr(xv, yv)
        plt.plot(xv, yv, '.', label=celltypes[ind]+' data, PCC ={0}, Cellnum={1}'.format(np.round(temp[0], 3), len(filt_cc1)), linestyle='None')
        xv, yv = np.asarray([obj.vb for obj in filt_cc2]) * scale ** 3, np.asarray(
            [obj.start * timestep for obj in filt_cc2])
        xmin = max(0.0, np.mean(xv) - 2.5 * np.std(xv))
        xvals = np.linspace(xmin, np.mean(xv) + 2.5 * np.std(xv), 20)
        temp = scipy.stats.pearsonr(xv, yv)
        plt.plot(xv, yv, '.',
                 label=celltypes[ind] + ' data, PCC ={0}, Cellnum={1}'.format(np.round(temp[0], 3), len(filt_cc2)),
                 linestyle='None')
        plt.xlabel('$V_b$ $(\mu m^3)$')
        plt.ylabel('$t_{Start}$ (minutes)')
        plt.title(celltypes[ind]+r' Start timing Correlations')
        plt.legend()
        # print temp
        directory = base_path+expt_path+'/plots'
        if not os.path.exists(directory):
            os.makedirs(directory)
        fig.savefig(directory+'/'+celltypes[ind]+'_vb_tstart.png', bbox_inches='tight', dpi=fig.dpi)


        fig = plt.figure(figsize=[5, 5])
        xv, yv = np.asarray([obj.vb for obj in filt_cc])*scale**3, \
                 np.asarray([(obj.frames[-1]-obj.start)*timestep for obj in filt_cc])
        xmin = max(0.0, np.mean(xv) - 2.5 * np.std(xv))
        xvals = np.linspace(xmin, np.mean(xv) + 2.5 * np.std(xv), 20)
        # plt.xlim(xmin=xmin)
        temp = scipy.stats.pearsonr(xv, yv)
        plt.plot(xv, yv, '.', label=celltypes[ind]+' data, PCC ={0}'.format(np.round(temp[0], 3)), linestyle='None')
        plt.xlabel('$V_b$ $(\mu m^3)$')
        plt.ylabel('$t_{budded}$ (minutes)')
        plt.title(celltypes[ind]+r' Budded timing Correlations, cellnum = {0}'.format(len(filt_cc)))
        temp = scipy.stats.linregress(xv, yv)
        # plt.plot(xvals, temp[0]*xvals+temp[1], label='slope = {0}'.format(str(np.round(temp[0],2))))
        plt.legend()
        # print temp
        directory = base_path+expt_path+'/plots'
        if not os.path.exists(directory):
            os.makedirs(directory)
        fig.savefig(directory+'/'+celltypes[ind]+'_vb_tbud.png', bbox_inches='tight', dpi=fig.dpi)


        fig = plt.figure(figsize=[5, 5])
        xv, yv = np.asarray([obj.ellipse_volume[obj.start] for obj in filt_cc])*scale**3, \
                 np.asarray([(len(obj.frames)-obj.start)*timestep for obj in filt_cc])
        xmin = max(0.0, np.mean(xv) - 2.5 * np.std(xv))
        xvals = np.linspace(xmin, np.mean(xv) + 2.5 * np.std(xv), 20)
        # plt.xlim(xmin=xmin)
        temp = scipy.stats.pearsonr(xv, yv)
        plt.plot(xv, yv, '.', label=celltypes[ind]+' data, PCC ={0}'.format(np.round(temp[0], 3)), linestyle='None')
        plt.xlabel('$V_s$ $(\mu m^3)$')
        plt.ylabel('$t_{budded}$ (minutes)')
        plt.title(celltypes[ind]+r' Budded timing Correlations, cellnum = {0}'.format(len(filt_cc)))
        temp = scipy.stats.linregress(xv, yv)
        # plt.plot(xvals, temp[0]*xvals+temp[1], label='slope = {0}'.format(str(np.round(temp[0],2))))
        plt.legend()
        # print temp
        directory = base_path+expt_path+'/plots'
        if not os.path.exists(directory):
            os.makedirs(directory)
        fig.savefig(directory+'/'+celltypes[ind]+'_vs_tbud.png', bbox_inches='tight', dpi=fig.dpi)


        fig = plt.figure(figsize=[5, 5])
        xv, yv = np.asarray([obj.int_fl[0] for obj in filt_cc]), np.asarray([obj.int_fl[-1]+obj.int_fl_bud[-1] for obj in filt_cc])
        xmin = max(0.0, np.mean(xv) - 2.5 * np.std(xv))
        xvals = np.linspace(xmin, np.mean(xv) + 2.5 * np.std(xv), 20)
        # plt.xlim(xmin=xmin)
        plt.plot(xv, yv, '.', label=celltypes[ind]+' data', linestyle='None')
        plt.xlabel('$F_b$ (A.U.)')
        plt.ylabel('$F_{d}$ (A.U.)')
        plt.title(celltypes[ind]+r' Fluorescence Correlations, cellnum = {0}'.format(len(filt_cc)))
        temp = scipy.stats.linregress(xv, yv)
        plt.plot(xvals, temp[0]*xvals+temp[1], label='slope = {0}'.format(str(np.round(temp[0], 2))))
        plt.legend()
        # print temp
        directory = base_path+expt_path+'/plots'
        if not os.path.exists(directory):
            os.makedirs(directory)
        fig.savefig(directory+'/'+celltypes[ind]+'_fb_fdiv.png', bbox_inches='tight', dpi=fig.dpi)


print len([obj for obj in cc if obj.label_type == 1]), len([obj for obj in cc if obj.label_type == 0]), len([obj for obj in cc if obj.label_type == -1])
print np.mean([obj.vb for obj in cc if obj.label_type == 1 and obj.complete])
print np.mean([obj.vb for obj in cc if obj.label_type == 0 and obj.complete])

import seaborn as sns
# Plotting the ratio of daughter to mother Whi5 fluorescence at birth
filt_cc = [obj for obj in cc if obj.complete and not(obj.error) and not(obj.daughter is None) and
           not(obj.next_gen is None)]
print len(filt_cc)
inds = [obj.index for obj in cc]
mother_start = np.asarray([cc[inds.index(obj.next_gen)].int_fl[0] for obj in filt_cc])
daughter_start = np.asarray([cc[inds.index(obj.daughter)].int_fl[0] for obj in filt_cc])
fig = plt.figure(figsize=[5, 5])
sns.distplot(daughter_start/mother_start)
plt.xlabel('Fl(daughter)/Fl(mother)')
plt.title('Number of cells {0}'.format(len(filt_cc)))
fig.savefig(directory+'/'+'fluorescence_ratios.png', bbox_inches='tight', dpi=fig.dpi)

mother_startv = np.asarray([cc[inds.index(obj.next_gen)].ellipse_volume[0] for obj in filt_cc])
daughter_startv = np.asarray([cc[inds.index(obj.daughter)].ellipse_volume[0] for obj in filt_cc])
fig = plt.figure(figsize=[5, 5])
sns.distplot(daughter_startv/mother_startv)
plt.xlabel('V(daughter)/V(mother)')
plt.title('Number of cells {0}'.format(len(filt_cc)))
fig.savefig(directory+'/'+'volume_ratios.png', bbox_inches='tight', dpi=fig.dpi)

fig = plt.figure(figsize=[5, 5])
sns.distplot(daughter_start*mother_startv/(daughter_startv*mother_start))
plt.xlabel('Conc(daughter)/Conc(mother)')
plt.title('Number of cells {0}'.format(len(filt_cc)))
fig.savefig(directory+'/'+'concentration_ratios.png', bbox_inches='tight', dpi=fig.dpi)