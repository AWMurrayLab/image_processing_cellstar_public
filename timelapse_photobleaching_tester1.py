import numpy as np
import custom_image_toolkit as C
import pandas as pd
import os
import cPickle as pickle
import skimage
from skimage import io
import scipy
import scipy.io as sio
import matplotlib.pyplot as plt

# expt_ids = ['/181917_yFB79_60X_Raff_125uMGal', '/181114_yFB78_Raff_125Gal', '/190607_yFB78_60X_Raff_125uMGal']
# '/181207_yFB79_60X_Raff_125uMGal'  # appears to have a corrupted .mat file for Scene 1 Frame 45
# expt_ids = ['/181114_yFB78_Raff_125Gal']
# expt_ids = ['/181917_yFB79_60X_Raff_125uMGal', '/190607_yFB78_60X_Raff_125uMGal', '/181114_yFB78_Raff_125Gal']
# expt_ids = ['/190607_yFB110_60X_Raff_125uMGal']
expt_ids = ['/190629_yFB110_60X_Raff_125uMGal']
fluors = ['mVenNB', 'mCherry']

for i0 in range(len(expt_ids)):
    pickle_in = open("./expt_ids"+expt_ids[i0]+'.pickle',"rb")
    ep = pickle.load(pickle_in)  # this gives us the experimental parameters, so that we can load everything in an
    temp_dir = ep['base_path'] + ep['expt_path']
    temp = np.load(temp_dir+'/photobleaching.npy')
    temp = temp / np.tile(np.absolute(temp[0,:,:,:]), [temp.shape[0],1,1,1])
    time = np.arange(temp.shape[0])*ep['tstep']
    for i1 in range(2):
        fig=plt.figure(figsize=[5,5])
        plt.plot(time, np.nanmean(temp, axis=1)[:, i1, 0],
                         label='Mean Fluorescence '+fluors[i1])
        plt.fill_between(time, np.nanmean(temp, axis=1)[:, i1, 0] - np.nanstd(temp, axis=1)[:, i1, 0],
                         np.nanmean(temp, axis=1)[:, i1, 0] + np.nanstd(temp, axis=1)[:, i1, 0],alpha=0.5)
        plt.plot(time, np.nanmean(temp, axis=1)[:, i1, 1], label='Median Fluorescence '+fluors[i1])
        plt.fill_between(time, np.nanmean(temp, axis=1)[:, i1, 1] - np.nanstd(temp, axis=1)[:, i1, 1],
                         np.nanmean(temp, axis=1)[:, i1, 1] + np.nanstd(temp, axis=1)[:, i1, 1],alpha=0.5)
        plt.ylabel('Cellular Fluorescence (normalized)')
        plt.xlabel('Time (mins)')
        plt.legend()
        plt.title(expt_ids[i0][1:]+', Timestep = {0} mins'.format(str(ep['tstep'])))
        fig.savefig(temp_dir+'/plots/photobleaching_'+fluors[i1]+'.png', bbox_inches='tight')
plt.clf()

for i0 in range(len(expt_ids)):
    pickle_in = open("./expt_ids"+expt_ids[i0]+'.pickle',"rb")
    ep = pickle.load(pickle_in)  # this gives us the experimental parameters, so that we can load everything in an
    temp_dir = ep['base_path'] + ep['expt_path']
    temp = np.load(temp_dir+'/photobleaching.npy')
    temp = temp / np.tile(np.absolute(temp[0,:,:,:]), [temp.shape[0],1,1,1])
    time = np.arange(temp.shape[0])*ep['tstep']
    for i1 in range(2):
        fig=plt.figure(figsize=[5,5])
        for i2 in range(temp.shape[1]):
            plt.plot(time, temp[:, i2, i1, 1], label='Median Fluorescence scene {0}'.format(i2))
        plt.ylabel('Cellular Fluorescence (normalized)')
        plt.xlabel('Time (mins)')
        # plt.legend(loc='right')
        plt.title(expt_ids[i0][1:]+', Timestep = {0} mins'.format(str(ep['tstep'])))
        fig.savefig(temp_dir+'/plots/photobleaching_fields_'+fluors[i1]+'.png', bbox_inches='tight')
plt.clf()
fluors = ['515nm', '594nm']
for i0 in range(len(expt_ids)):
    pickle_in = open("./expt_ids"+expt_ids[i0]+'.pickle',"rb")
    ep = pickle.load(pickle_in)  # this gives us the experimental parameters, so that we can load everything in an
    temp_dir = ep['base_path'] + ep['expt_path']
    temp = np.load(temp_dir+'/bkgd_brightness.npy')
    # print temp[0, :, :]
    temp = temp / np.tile(np.absolute(temp[0,:,:,:]), [temp.shape[0],1,1,1])


    time = np.arange(temp.shape[0])*ep['tstep']
    for i1 in range(2):
        fig=plt.figure(figsize=[5,5])
        plt.plot(time, np.nanmean(temp,axis=1)[:, i1, 0],
                         label='Mean Fluorescence bkgd '+fluors[i1])
        plt.fill_between(time, np.nanmean(temp, axis=1)[:, i1, 0] - np.nanstd(temp, axis=1)[:, i1, 0],
                         np.nanmean(temp, axis=1)[:, i1, 0] + np.nanstd(temp, axis=1)[:, i1, 0], alpha=0.5)
        plt.plot(time, np.nanmean(temp,axis=1)[:, i1, 1], label='Median Fluorescence bkgd '+fluors[i1])
        plt.fill_between(time, np.nanmean(temp, axis=1)[:, i1, 1] - np.nanstd(temp, axis=1)[:, i1, 1],
                         np.nanmean(temp, axis=1)[:, i1, 1] + np.nanstd(temp, axis=1)[:, i1, 1], alpha=0.5)
        plt.ylabel('Background Fluorescence (normalized)')
        plt.xlabel('Time (mins)')
        plt.legend()
        plt.title(expt_ids[i0][1:]+', Timestep = {0} mins'.format(str(ep['tstep'])))
        fig.savefig(temp_dir+'/plots/photobleaching_bkgd_'+fluors[i1]+'.png', bbox_inches='tight')
