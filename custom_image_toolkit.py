import numpy as np
import scipy
import scipy.optimize
import scipy.io as sio
import weakref
import bisect
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from scipy import stats
import pandas as pd
import cPickle as pickle
import skimage
from skimage import measure
import math
from skimage import io
import os
from matplotlib import cm
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from shutil import copyfile
from skimage import feature
from skimage import filters
import seaborn as sns


# helper
fluor_chars_names = ['Skew', 'CV', 'Gini']


class Cell(object):
    cellCount = 0  # total number of cells

    def __init__(self, birth_parameters):  # birth parameters = [tb, celltype, parent, vb, parent_current]
        self.exists = True
        # these are present for all cells
        self.index = birth_parameters['index']  # note this starts from 1 as is defined by cellstar
        self.frames = [birth_parameters['current_frame']]
        self.position = [birth_parameters['centroid']]  # X then Y
        self.index_image = [birth_parameters['index_image']]
        self.ellipse_fit = [birth_parameters['ellipse_params']]  # major, minor, orientation, centroid
        self.ellipse_volume = [birth_parameters['ellipse_volume']]
        self.area = [birth_parameters['2d_area']]
        self.nuclear_whi5 = [0]  # we will change this later depending on fluorescence data
        self.type = -1  # if cells are undetermined we label them a -1
        self.fluor_chars = [0]  # fluorescence characteristics at each time point
        self.segment_coords = [birth_parameters['segment_coords']]  # in format (y, x)  (row, col)
        self.zproj_fluor_vals = [0]
        self.zproj_fluor_vals_c2 = [0]
        self.nuclear_coords = [None]
        self.nuclear_fluor_int = [np.nan]  # make these nan values so that if we don't have actual data we just exclude
        self.nuclear_fluor_av = [np.nan]  # make these nan values so that if we don't have actual data we just exclude
        self.cytoplasmic_fluor_av = [np.nan]
        self.cytoplasmic_fluor_int = [np.nan]
        self.nuclear_fluor_int_c2 = [np.nan]
        self.nuclear_fluor_av_c2 = [np.nan]
        self.cytoplasmic_fluor_int_c2 = [np.nan]
        self.cytoplasmic_fluor_av_c2 = [np.nan]

        # if
        # self.
        # self.celltype = birth_parameters[1]  # 0 is mother, 1 is daughter
        # self.parent = birth_parameters[2]
        # self.parent_current = birth_parameters[4]  # only used for daughter cells
        Cell.cellCount += 1

    def next_frame(self, frame_parameters):
        self.frames.append(frame_parameters['current_frame'])  # adding the current frame to the list of frames for this
        # cell
        self.position.append(frame_parameters['centroid'])
        self.index_image.append(frame_parameters['index_image'])
        self.ellipse_fit.append(frame_parameters['ellipse_params'])  # major, minor, orientation, centroid
        self.ellipse_volume.append(frame_parameters['ellipse_volume'])
        self.area.append(frame_parameters['2d_area'])
        self.nuclear_whi5.append(0)
        self.fluor_chars.append(0)  # appending a characteristics list at each timepoint
        self.zproj_fluor_vals.append(0)
        self.zproj_fluor_vals_c2.append(0)
        self.nuclear_coords.append(None)
        self.nuclear_fluor_int.append(np.nan)
        self.nuclear_fluor_av.append(np.nan)
        self.cytoplasmic_fluor_av.append(np.nan)
        self.cytoplasmic_fluor_int.append(np.nan)
        self.nuclear_fluor_int_c2.append(np.nan)
        self.nuclear_fluor_av_c2.append(np.nan)
        self.cytoplasmic_fluor_int_c2.append(np.nan)
        self.cytoplasmic_fluor_av_c2.append(np.nan)
        self.segment_coords.append(frame_parameters['segment_coords'])

    def add_fluor_placeholders(self):
        # self.int_fl = [] # deprecated
        self.coords = []
        self.pixel_thresh_coords = []
        self.pixel_thresh_fluor_vals_c2 = []
        self.pixel_thresh_fluor_vals = []
        for temp_ind in range(len(self.frames)):
            self.coords.append([])
            self.pixel_thresh_coords.append([0])
            self.pixel_thresh_fluor_vals_c2.append(0)
            self.pixel_thresh_fluor_vals.append(0)


class CellCycle(object):
    cellCount = 0  # total number of cell cycles tracked

    def __init__(self, temp_cell, temp_parameters):  # generating a placeholder cell which will be fully updated later.
        self.index = CellCycle.cellCount  # this will be the unique index of this cell cycle
        # self.index = temp_parameters['index']  # this will be the unique index of this cell cycle
        self.cell_index = temp_parameters['cell_index']  # this is the unique identifier of the cell this came from
        self.cc_num = temp_parameters['cc_num']  # this is the number of the cell cycle in the cell this came from
        self.complete = temp_parameters['complete']  # whether this corresponds to a full cell cycle or not
        self.frames = [temp_cell.frames[temp_ind] for temp_ind in temp_parameters['range']]
        self.ellipse_volume = [temp_cell.ellipse_volume[temp_ind] for temp_ind in temp_parameters['range']]
        # self.int_fl = [temp_cell.int_fl[temp_ind] for temp_ind in temp_parameters['range']] # deprecated
        self.ellipse_fit = [temp_cell.ellipse_fit[temp_ind] for temp_ind in
                            temp_parameters['range']]  # major, minor, orientation, centroid
        self.data_origin = temp_parameters['data_origin']
        self.label_type = temp_cell.type
        self.nuclear_whi5 = [temp_cell.nuclear_whi5[temp_ind] for temp_ind in temp_parameters['range']]
        self.zproj_fl = [temp_cell.zproj_fluor_vals[temp_ind] for temp_ind in temp_parameters['range']]
        self.zproj_fl_c2 = [temp_cell.zproj_fluor_vals_c2[temp_ind] for temp_ind in temp_parameters['range']]
        self.pixel_thresh_coords = [temp_cell.pixel_thresh_coords[temp_ind] for temp_ind in temp_parameters['range']]
        self.pixel_thresh_fluor_vals_c2 = [temp_cell.pixel_thresh_fluor_vals_c2[temp_ind] for temp_ind in
                                           temp_parameters['range']]
        self.pixel_thresh_fluor_vals = [temp_cell.pixel_thresh_fluor_vals[temp_ind] for temp_ind in
                                        temp_parameters['range']]
        self.nuclear_coords = [temp_cell.nuclear_coords[temp_ind] for temp_ind in temp_parameters['range']]
        self.nuclear_fluor_int = [temp_cell.nuclear_fluor_int[temp_ind] for temp_ind in temp_parameters['range']]
        self.nuclear_fluor_av = [temp_cell.nuclear_fluor_av[temp_ind] for temp_ind in temp_parameters['range']]
        self.cytoplasmic_fluor_av = [temp_cell.cytoplasmic_fluor_av[temp_ind] for temp_ind in temp_parameters['range']]
        self.cytoplasmic_fluor_int = [temp_cell.cytoplasmic_fluor_int[temp_ind] for temp_ind in temp_parameters['range']]
        self.nuclear_fluor_int_c2 = [temp_cell.nuclear_fluor_int_c2[temp_ind] for temp_ind in
                                      temp_parameters['range']]
        self.nuclear_fluor_av_c2 = [temp_cell.nuclear_fluor_av_c2[temp_ind] for temp_ind in
                                     temp_parameters['range']]
        self.cytoplasmic_fluor_int_c2 = [temp_cell.cytoplasmic_fluor_int_c2[temp_ind] for temp_ind in
                                    temp_parameters['range']]
        self.cytoplasmic_fluor_av_c2 = [temp_cell.cytoplasmic_fluor_av_c2[temp_ind] for temp_ind in
                                    temp_parameters['range']]

        # print [len(obj) for obj in self.pixel_thresh_coords], self.data_origin, self.cell_index, self.frames
        self.pixel_thresh_vol = []
        for temp_ind in temp_parameters['range']:
            if len(temp_cell.pixel_thresh_coords[temp_ind]) == 3:
                self.pixel_thresh_vol.append(len(temp_cell.pixel_thresh_coords[temp_ind][0]))
            else:
                self.pixel_thresh_vol.append(0)
                print 'Pixel thresholds unknown for Cell cycle {0}, cell {1}, frame {2}'.format(self.index,
                                                                                                self.cell_index,
                                                                                                temp_cell.frames[
                                                                                                    temp_ind])
        # self.pixel_thresh_vol = [len(temp[0]) for temp in self.pixel_thresh_coords]  # needs to be scaled by volume
        # of pixel in 3D.

        self.segment_coords = [temp_cell.segment_coords[temp_ind] for temp_ind in temp_parameters['range']]
        if temp_parameters['complete']:  # if this is a full cell cycle.
            self.tb = self.frames[0]
            self.td = self.frames[-1]
            self.vb = self.ellipse_volume[0]
            self.vd = self.ellipse_volume[-1]
        self.start = temp_parameters['start']  # gives the frame at which Start begins.
        self.error = temp_parameters['error']

        # lineage variables
        self.daughter = None  # gives the index of the daughter cell if it exists
        self.celltype = -1  # unsure of celltype to begin with
        self.next_gen = None  # gives the index of the next generation of this cell cycle if it exists
        self.prev_gen = None  # gives the index of the previous generation of this cell cycle if it exists
        self.bud = None  # gives the index of the bud for this cell if it exists
        self.parent = None  # gives the unique index of the cell cycle for the "parent" cell in which it produced the
        # current cell (if there is one)
        self.family = None
        CellCycle.cellCount += 1

    def update_cycle(self, temp_cell, temp_parameters):  # this is what we run to get full cell cycle info
        self.celltype = temp_parameters['celltype']  # this is zero if the cell cycle is a mother, 1 if a daughter, 2
        # if a bud, -1 if unsure
        self.bud = temp_parameters['bud']  # unique index of bud cell
        self.next_gen = temp_parameters['next_gen']  # unique index of next generation
        self.parent = temp_parameters['parent']  # this is None if there was no parent assigned.


def import_segments(path):
    temp = sio.loadmat(path)
    return temp


def single_frame(temp_path, temp_cells, temp_frame, tracking_csv, color_seq):
    # This function takes an existing set of cells, temp_cells, and will append cells from a given frame to them.
    # temp_path is the path of the .mat file with segments that is output automatically for each field of view from
    # cell star. temp_frame is an index for the frame number in question (starts from 1). tracking_csv is the path of
    # the tracking .csv file output by cellstar automatically. color_seq is the range of values in the color map.

    temp_im = sio.loadmat(temp_path)['segments']
    temp1 = np.amax(temp_im)
    # print temp1  # number of segmented cells in this image
    fig = plt.figure(figsize=[5.12, 5.12], frameon=False)  # note figsize is selected to scale with that of the image
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    temp_im1 = np.zeros(temp_im.shape)
    # plt.imshow(temp_im)
    temp_coords = []
    temp_outlines = np.zeros(temp_im.shape)
    for i0 in range(1, temp1 + 1):  # this just indexes the cells in the image, not the unique cell identifier
        # Note that this index appears to (very frustratingly) be different for the tracking_csv index and the value in
        # the image
        # print i0
        temp2 = temp_im == i0  # go through sequentially
        # calculating ellipse behavior
        temp_val = fit_ellipse(temp2.astype(int))
        temp_coords.append(temp_val[4])
        # print temp_coords[-1]
        # plt.plot(, color=(1.0, 1.0, 1.0, 1.0))
        temp_cell_index = find_min_distance(tracking_csv[(tracking_csv['Frame_number'] == temp_frame + 1)], temp_val[3])
        # assign tracks based on proximity of centroid
        temp3 = tracking_csv[(tracking_csv['Frame_number'] == temp_frame + 1) & (tracking_csv[' Unique_cell_number'] ==
                                                                                 temp_cell_index)]
        # Now we have the data for the cell with that proximity
        # temp3 = tracking_csv[(tracking_csv['Frame_number']==temp_frame+1) & (tracking_csv[' Cell_number']==i0)]
        # note that temp_frame has 1 added here to account for the first frame being the background
        # print c
        temp4 = {}
        temp4['index'] = temp3[' Unique_cell_number'].iloc[0]  # unique cell identifier
        temp4['index_image'] = temp3[' Cell_number'].iloc[0]  # cell identifier in that image
        temp4['file_centroid'] = np.array([temp3[' Position_X'].iloc[0], temp3[' Position_Y'].iloc[0]])  # cell centroid
        # stored as x, y
        temp4['segment_coords'] = temp_val[5]
        temp4['centroid'] = np.array(temp_val[3])
        if np.linalg.norm(temp4['file_centroid'] - temp4['centroid']) > 20:
            print temp_frame, temp4['index'], temp4['index_image'], i0
            print 'stored centroid', temp4['file_centroid']
            print 'measured centroid', temp4['centroid']
            raise ValueError('Cell classification is different somehow')
        temp4['current_frame'] = temp_frame  # adding the current frame index
        temp4['ellipse_params'] = temp_val[:4]
        temp4['ellipse_volume'] = math.pi * temp_val[0] * temp_val[1] ** 2
        temp4['2d_area'] = np.sum(temp2)  # gives the total number of pixels in the 2d plane
        if (len(temp_cells) == 0) or (not (temp3[' Unique_cell_number'].iloc[0] in [obj.index for obj in temp_cells])):
            # if this cell has not been added yet
            temp_cells.append(Cell(temp4))  # adding a new cell
        else:
            temp_ind = [obj.index for obj in temp_cells].index(temp3[' Unique_cell_number'].iloc[0])
            if ~(temp_cells[temp_ind].index == temp3[' Unique_cell_number'].iloc[0]):
                raise ValueError('Selected cell is incorrect')
            temp_cells[temp_ind].next_frame(frame_parameters=temp4)
        temp_im1[temp4['segment_coords'][:, 0], temp4['segment_coords'][:, 1]] = temp4['index']
        # making the plot of cell values
    ax.imshow(temp_im1)
    for coords in temp_coords:
        temp_x, temp_y = zip(*coords)
        # print len(temp_y)
        # print list(temp_y), list(temp_x)
        # exit()
        plt.plot(list(temp_x), list(temp_y), color=(1.0, 1.0, 1.0, 1.0))
        temp_inds1 = [i for i, e in enumerate(temp_y) if 0 > int(round(e)) or int(round(e)) > 511]
        temp_inds2 = [i for i, e in enumerate(temp_x) if 0 > int(round(e)) or int(round(e)) > 511]
        excl_inds = list(set(temp_inds1 + temp_inds2))
        # print excl_inds
        corr_inds = [i0 for i0 in range(len(temp_x)) if not (i0 in excl_inds)]
        # print corr_inds
        temp_x1 = [int(round(temp_x[i0])) for i0 in corr_inds]
        temp_y1 = [int(round(temp_y[i0])) for i0 in corr_inds]
        # print temp_outlines.shape
        temp_outlines[temp_y1, temp_x1] = 1
        # exit()
        # temp_coords[*zip(*coords)]
        # plt.plot(obj.position[-1][0], obj.position[-1][1], '.', color=color_seq[(obj.index-1)%11])\
    plt.xlim(-0.5, 511.5)
    plt.ylim(511.5, -0.5)
    return temp_cells, ax, fig, temp_outlines


def add_fluorescence_traces_v2(temp_path, temp_cells, frame_list, current_frame, z_scaling, z_offset, temp_bf_mask,
                               save_coords=False, exists_c2=None, temp_outlines=None, temp_im_c2=None):
    # This is the updated version of add_fluorescence_traces_v1, where we only consider the fluorescence values in the
    # full z-stack. We no longer store the full fluorescence distribution either, instead extracting relevant
    # characteristics from it such as the skewness, genie coefficient and standard deviation for the cell classification
    # z_scaling is deprecated but gives the ratio of the z stepsize in um to the x,y pixel size in um.
    # temp_outlines gives the outlines in format [y,x]
    # temp_im_c2 is still available at this stage to allow us to include data for the secondary channel in the nucleus
    temp_im1 = io.imread(temp_path)
    temp_bf_mask = np.tile(temp_bf_mask, [temp_im1.shape[0], 1, 1])
    bkgd = np.median(temp_im1[np.nonzero(temp_bf_mask==0)])  # gives the median background fluorescence value
    temp_im = temp_im1 - bkgd
    # subtracting the average background at each xy point. Structure is (z,y,x)
    # finding the blobs
    temp_im2 = skimage.filters.gaussian(temp_im1, sigma=1)
    temp = skimage.feature.blob_log(temp_im2, min_sigma=1.0, max_sigma=2.5, num_sigma=4, threshold=0.0015, overlap=0.01,
                                    log_scale=False)  # this has the format of a numpy array Nx4, where the columns are
    # z, y, x, r where r is the radius of the circle centered at that point
    # now we plot the fitted blobs along with the outlines of the cells at this timepoint
    fig = plt.figure(figsize=[10,10])
    ax = plt.subplot(1, 1, 1)
    plt.imshow(np.amax(temp_im2, axis=0))
    # plotting the cell outlines in white
    temp_coords = np.nonzero(temp_outlines[current_frame-1, :, :])
    plt.plot(temp_coords[1], temp_coords[0], '.', markersize=1, color=(1.0, 1.0, 1.0, 1.0))  # plotting in format x,y
    del temp_im2, temp_coords
    # plotting the circles of each blob
    for blob in temp:
        #     print blob
        z, y, x, r = blob
        c = plt.Circle((x, y), np.sqrt(3) * r, color='r', linewidth=2, fill=False)
        ax.add_patch(c)
    fig.subplots_adjust(bottom=0)
    fig.subplots_adjust(top=1)
    fig.subplots_adjust(right=1)
    fig.subplots_adjust(left=0)

    temp_yx = zip(*[temp[:, 1], temp[:, 2]])  # list of coordinates in type tuple
    # now we go through the set of cells in the current timepoint and determine which ones have this x,y point inside,
    # in addition to determining the fluorescence within each cell and storing that.
    for temp_ind in frame_list:
        i0 = temp_cells[temp_ind].frames.index(current_frame)
        # gives the index for the current frame in that cell's history
        # if current_frame != temp_cells[temp_ind].frames[i0]:
        #     print(i0, current_frame, temp_cells[temp_ind].frames[i0], temp_cells[temp_ind].frames)
        #     raise ValueError('Wrong frame selected')

        temp_coords = temp_cells[temp_ind].segment_coords[i0]  # the segmented coordinates. NOTE OVERWRITING TEMP_COORDS
        # FROM EARLIER. NO LONGER OUTLINES, NOW FULL 2D CELL MASK COORDS.
        # saving relevant fluorescence statistics
        temp_chars = list([])
        # print temp_coords
        temp_chars.append(scipy.stats.skew(temp_im[:, temp_coords[:, 0], temp_coords[:, 1]].flatten()))
        # skew
        temp_chars.append(np.std(temp_im[:, temp_coords[:, 0], temp_coords[:, 1]].flatten())/np.absolute(np.mean(temp_im[:, temp_coords[:, 0], temp_coords[:, 1]].flatten())))
        # CV
        temp_chars.append(gini(temp_im[:, temp_coords[:, 0], temp_coords[:, 1]]))  # Gini coefficient
        # print temp_chars
        temp_cells[temp_ind].fluor_chars[i0] = temp_chars
        # print temp_cells[temp_ind].fluor_chars
        # exit()
        # generating a list of fluorescence values to get good statistics. Note image formatting with z, y, x.
        temp_cells[temp_ind].zproj_fluor_vals[i0] = np.sum(temp_im[:, temp_coords[:, 0], temp_coords[:, 1]])
        # calculating the full z projection as well.
        if not (exists_c2 is None):  # if we have the data then we have to also add fluorescence integrated from
            # pixels segmented with respect to the second fluorescence channel
            temp_cells[temp_ind].pixel_thresh_fluor_vals[i0] = \
                np.sum(temp_im[temp_cells[temp_ind].pixel_thresh_coords[i0]])

        # Saving whether there was a blob inside the current cell and if so what its coords were
        temp_blobs = list(
            set(zip(*[temp_coords[:, 0], temp_coords[:, 1]])) & set(temp_yx))  # the list of coordinates that are in
        # both temp_coords and temp_yx. This gives the list of all blobs whose centers lie within the current cell
        # being considered.

        # code = coords of the one blob OR None if no blobs or more than one blob
        if len(temp_blobs) == 1:  # if there is exactly one blob within the current cell
            temp_deets = temp[temp_yx.index(temp_blobs[0]), :]
            # Saving the fluorescence within the nucleus. This part is computationally intensive but would require a lot of
            # work to generate something more effective

            # generating a cuboid grid around the center of the blob
            Z, Y, X = np.meshgrid(
                np.arange(np.maximum(temp_deets[0] - np.sqrt(3)*temp_deets[3], 0),
                          np.minimum(temp_deets[0] + np.sqrt(3)*temp_deets[3], temp_im.shape[0])).astype(int),
                np.arange(np.maximum(temp_deets[1] - np.sqrt(3) * temp_deets[3], 0),
                np.minimum(temp_deets[1] + np.sqrt(3) * temp_deets[3], temp_im.shape[1])).astype(int),
                np.arange(np.maximum(temp_deets[2] - np.sqrt(3) * temp_deets[3], 0),
                np.minimum(temp_deets[2] + np.sqrt(3) * temp_deets[3], temp_im.shape[2])).astype(int))

            # gives the list of coordinates within the sphere centered at the middle of the current blob
            coord_scaled = zip(Z.ravel(), Y.ravel(), X.ravel())
            temp_coords1 = [a for a in coord_scaled if np.linalg.norm(a - temp_deets[:3])<np.sqrt(3)*temp_deets[3]]
            coord_scaled1 = zip(*temp_coords1)  # this is what we will use to calculate the relevant values
            temp_cells[temp_ind].nuclear_coords[i0] = temp_coords1  # saving the coordinates of the nucleus at this
            # time point for use with the secondary fluorescence channel
            temp_cells[temp_ind].nuclear_fluor_int[i0] = np.sum(temp_im[coord_scaled1])
            temp_cells[temp_ind].nuclear_fluor_av[i0] = np.mean(temp_im[coord_scaled1])
            if not(exists_c2 is None):  # trying to use this data to infer values about the cytoplasm of the cell
                temp_cells[temp_ind].nuclear_fluor_int_c2[i0] = np.sum(temp_im_c2[coord_scaled1])  # nuclear fluor in c2
                temp_cells[temp_ind].nuclear_fluor_av_c2[i0] = np.mean(temp_im_c2[coord_scaled1])  # nuclear fluor in c2
                temp_coords2 = zip(*list(set(zip(*temp_cells[temp_ind].pixel_thresh_coords[i0])) - set(temp_coords1)))
                # cytoplasmic coordinates for this cell in format ([z], [y], [x]) that can be fed directly into np array
                # Calculating the cytoplasmic fluorescence values for this cell, both fluorescence channels
                temp_cells[temp_ind].cytoplasmic_fluor_int[i0] = np.sum(temp_im[temp_coords2])
                temp_cells[temp_ind].cytoplasmic_fluor_av[i0] = np.mean(temp_im[temp_coords2])
                temp_cells[temp_ind].cytoplasmic_fluor_int_c2[i0] = np.sum(temp_im_c2[temp_coords2])
                temp_cells[temp_ind].cytoplasmic_fluor_av_c2[i0] = np.mean(temp_im_c2[temp_coords2])
                # The majority of healthy cells have their nuclei almost completely contained within pixel_thresh_coords
                # at least within the first timepoint.
                # print temp_coords2
                # print 'Full list', temp_coords1
                # print zip(*temp_cells[temp_ind].pixel_thresh_coords[i0])
                # print temp_coords
    print('Number of cells = {0}'.format(temp_ind))
    return temp_cells, fig


def add_fluorescence_traces_c2_v1(temp_path, temp_cells, frame_list, current_frame, temp_bf_mask):
    # v1 differs from the base by also analyzing different ways of segmenting cell volume based on fluorescence
    # This function adds the z-projected fluorescence above and below the 2d cell segmentation for a second fluorescence
    # channel c2
    temp_im = io.imread(temp_path)
    temp_bf_mask = np.tile(temp_bf_mask, [temp_im.shape[0], 1, 1])
    # print np.amax(temp_bf_mask), np.amin(temp_bf_mask), current_frame
    temp_vals = temp_im[np.nonzero(temp_bf_mask==0)]
    median_bkgd, std_bkgd = np.median(temp_vals), np.median(temp_vals)-np.percentile(temp_vals, 16)
    temp_bkgd = median_bkgd + 4.0 * std_bkgd
    temp_mask = temp_im > temp_bkgd
    temp_bin1 = np.zeros(temp_im.shape)
    # generating the high pixel density mask
    temp_im1 = temp_im - median_bkgd  # subtracting average background measured over the whole image
    figs = []
    for temp_ind in frame_list:  # iterating through each cell to add the appropriate value
        temp_binary_image = np.zeros(temp_im.shape)
        i0 = temp_cells[temp_ind].frames.index(current_frame)
        # i0 = [i for i, e in enumerate(temp_cells[temp_ind].frames) if e == current_frame][0]
        if current_frame != temp_cells[temp_ind].frames[i0]:
            print i0, current_frame, temp_cells[temp_ind].frames[i0], temp_cells[temp_ind].frames
            raise ValueError('Wrong frame selected')
        # gives the index for the current frame in that cell's history
        temp_coords = temp_cells[temp_ind].segment_coords[i0]
        temp_cells[temp_ind].zproj_fluor_vals_c2[i0] = np.sum(temp_im1[:, temp_coords[:, 0], temp_coords[:, 1]])
        # calculating the full z projection

        # Generating the fluorescence and volume value integrated only over pixels where the fluorescence is high. Note
        # we will also do this for Whi5, since the nucleus appears to be included in the constitutive fluor channel, so
        # this should give a slightly more accurate measure for Whi5 fluorescence.
        # while the vacuoles are excluded for both
        temp_binary_image[:, temp_coords[:, 0], temp_coords[:, 1]] = 1
        temp_bin1+=temp_binary_image
        temp_binary_image *= temp_mask
        temp_cells[temp_ind].pixel_thresh_coords[i0] = np.nonzero(temp_binary_image)
        temp_cells[temp_ind].pixel_thresh_fluor_vals_c2[i0] = \
            np.sum(temp_im1[temp_cells[temp_ind].pixel_thresh_coords[i0]])
    temp_bin1=temp_bin1>0
    temp_mask*=temp_bin1  # gives the cell segments in this frame
    # plt.imshow(temp_mask[7, :, :])
    # plt.show()
    return temp_cells, figs, temp_mask, temp_im1


def find_min_distance(temp_array, temp_centroid):
    # returns the unique cell number with stored centroid closest to the calculated centroid
    # print temp_array.columns[3:]
    temp1 = temp_array.as_matrix(columns=temp_array.columns[3:])
    temp2 = temp1[:, :2]
    # print temp2
    temp3 = temp1[:, 2]
    # print temp2.shape, temp3.shape
    temp2[:, 0] = temp2[:, 0] - temp_centroid[0]
    temp2[:, 1] = temp2[:, 1] - temp_centroid[1]
    temp_dists = np.linalg.norm(temp2, axis=1)
    # print temp_dists
    temp_index = temp3[np.argmin(temp_dists)]
    # print temp_centroid, temp1[np.argmin(temp_dists), :]
    # print temp_index
    return temp_index


def rotation(temp_angle):
    # defines the SO(2) matrix for rotation of angle theta in anticlockwise direction.
    temp = np.array([[np.cos(temp_angle), -np.sin(temp_angle)], [np.sin(temp_angle), np.cos(temp_angle)]])
    return temp


def rotation_3d(temp_angle):
    # defines the SO(2) matrix for rotation of angle theta in anticlockwise direction.
    temp = np.array(
        [[np.cos(temp_angle), -np.sin(temp_angle), 0], [np.sin(temp_angle), np.cos(temp_angle), 0], [0, 0, 1]])
    return temp


def fit_ellipse(temp_im):
    temp = skimage.measure.regionprops(temp_im)
    if len(temp) != 1:
        raise ValueError('Too many regions')
    # http://scikit-image.org/docs/dev/api/skimage.measure.html
    temp1 = temp[0].major_axis_length / 2.0
    temp2 = temp[0].minor_axis_length / 2.0
    if temp1<temp2:
        raise ValueError('Ellipse volume is incorrect')
    temp3 = -temp[0].orientation  # note using negative angle here since y counts down, not up in our reckoning
    temp4 = temp[0].centroid
    temp_mat = rotation(temp3)
    # temp_coords = [reversed(np.dot(temp_mat, (temp1*np.cos(theta), temp2*np.sin(theta))) + temp4) for theta in np.linspace(0.0, 2*math.pi, 20)]
    temp_coords = [(np.dot(temp_mat, (temp1 * np.cos(theta), temp2 * np.sin(theta))) + temp4)[::-1] for theta in
                   np.linspace(0.0, 2 * math.pi, 20)]
    # print temp1, temp2, temp3
    return [temp1, temp2, temp3, temp4[::-1], temp_coords,
            temp[0].coords]  # major length, minor length, orientation, centroid,
    #  (centroid should be x, y), coordinates of ellipse in format (y, x)


def generate_ellipse_coords(temp_vals):
    # returns coordinates in [(X, Y)] format
    temp1 = temp_vals[0]
    temp2 = temp_vals[1]
    temp3 = temp_vals[2]  # note using negative angle here since y counts down, not up in our reckoning
    temp4 = temp_vals[3]
    temp_mat = rotation(temp_vals[2])
    temp_coords = [(np.dot(temp_mat, (temp1 * np.cos(theta), temp2 * np.sin(theta))) + temp4) for theta in
                   np.linspace(0.0, 2 * math.pi, 20)]
    return temp_coords


def save_object(obj, filename):
    # Code taken from:
    # https://stackoverflow.com/questions/4529815/saving-an-object-data-persistence?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

        # saving multiple objects
        # tech_companies = [
        #     Company('Apple', 114.18), Company('Google', 908.60), Company('Microsoft', 69.18)
        # ]
        # save_object(tech_companies, 'tech_companies.pkl')
        # restoring multiple objects
        # with open('tech_companies.pkl', 'rb') as input:
        #     tech_companies = pickle.load(input)


def assign_labels_1(temp_cells, temp_inds, temp_coords):
    # we expect that each coordinate in temp_coords should lie within a cell outline in temp_cells
    # temp_coords should be given in x,y format
    temp_labels = []
    temp_vals = []
    temp_centres = []
    temp_scales = []
    temp_mat = []
    assignments_made = np.zeros(len(temp_coords))
    for temp_ind1 in temp_inds:
        if temp_cells[temp_ind1].frames[0] != 1:
            print temp_cells[temp_ind1].frames[0], temp_ind1
            raise ValueError('Wrong cells')  # we want the first frame
        temp_vals.append(temp_cells[temp_ind1].ellipse_fit[0])  # getting the ellipse properties
        temp_centres.append(np.array([temp_vals[-1][3][0], temp_vals[-1][3][1]]))  # x,y format
        temp_scales.append(np.array([temp_vals[-1][0], temp_vals[-1][1]]))
        temp_mat.append(rotation(-temp_vals[-1][2]))  # inverse rotation about Z axis to determine viable points
        temp_coords1 = [np.linalg.norm((np.dot(temp_mat[-1], (a - temp_centres[-1])) / temp_scales[-1])) for a in
                        temp_coords]
        # rotating and finding scaled distance of each labeled point from the centroid of the cell in question
        temp_labels.append([i for i, e in enumerate(temp_coords1) if e < 1])
        if len(temp_labels[-1]) > 1:
            print 'multiple points assigned to a single cell', temp_labels[-1], temp_cells[temp_ind1].index
        for i, e in enumerate(temp_coords1):
            if e < 1:
                assignments_made[i] += 1  # tracking whether each point has been assigned to a cell.
        if len(temp_labels[-1]) > 0:
            temp_cells[temp_ind1].type = 1  # cells with a label in the first timeframe are assigned a 1
        else:
            temp_cells[temp_ind1].type = 0  # cells without a label in the first timeframe are assigned a 0
    temp1 = np.where(assignments_made == 0)[0]
    if len(temp1) > 0:
        print 'The following points were not assigned:'
        for i0 in temp1:
            print temp_coords[i0]
    return temp_cells


def assign_labels_2(temp_cells, temp_inds, temp_coords):
    # This version takes the minimum distance from centroids. This should be better when segmentation isn't great.
    # we expect that each coordinate in temp_coords should lie within a cell outline in temp_cells
    # temp_coords should be given in x,y format
    # print temp_coords,
    temp_vals = []
    temp_centres = []
    temp_labels = []
    assignments_made = np.zeros(len(temp_coords))
    for temp_ind1 in temp_inds:
        if temp_cells[temp_ind1].frames[0] != 1:
            print temp_cells[temp_ind1].frames[0], temp_ind1
            raise ValueError('Wrong cells')  # we want the first frame
        temp_cells[temp_ind1].type = 0  # first frame cells are all type zero unless labeled
        temp_vals.append(temp_cells[temp_ind1].ellipse_fit[0])  # getting the ellipse properties
        temp_centres.append(np.array([temp_vals[-1][3][0], temp_vals[-1][3][1]]))  # x,y format
    for temp_ind2 in range(len(temp_coords)):  # in this case we must do the assignment wrt temp_coords
        temp_coords1 = [np.linalg.norm(a - temp_coords[temp_ind2]) for a in temp_centres]
        if np.amin(temp_coords1) < 20:  # if this distance is less than the average grown cell diameter
            temp_cells[temp_inds[np.argmin(temp_coords1)]].type = 1
            temp_labels.append(np.argmin(temp_coords1))
            assignments_made[temp_ind2] = 1
        else:
            print "Too far", temp_coords[temp_ind2]
    print 'Assigned points'
    for i0 in temp_labels:
        print temp_centres[i0]
    temp1 = np.where(assignments_made == 0)[0]
    if len(temp1) > 0:
        print 'The following points were not assigned:'
        for i0 in temp1:
            print temp_coords[i0]
    return temp_cells


def assign_labels_3(temp_cells, temp_inds, temp_coords, temp_frame_num):
    # we expect that each coordinate in temp_coords should lie within a cell outline in temp_cells
    # temp_coords should be given as a list of tuples in (x,y) format
    # temp_inds should be a list of coordinates in temp_cells for cells which are present at this timestep, denoted by
    # temp_frame_num

    temp_labels = []
    temp_vals = []

    temp_timepoints = []  # temp_timepoints corresponds to the current timepoint in each cell's reference frame
    assignments_made = []  # this will detail which cells this point fits within
    for i0 in temp_inds:
        temp_timepoints.append(temp_cells[i0].frames.index(temp_frame_num))
        # gives the index at which this frame sits for the cell in question
        if temp_cells[i0].frames[temp_timepoints[-1]] != temp_frame_num:  # confirming that this is the correct index
            print temp_frame_num, temp_timepoints[-1], temp_cells[i0].frames
            raise ValueError('Wrong cells')

    for temp_ind1 in range(len(temp_coords)):

        temp_comparison = []

        # assessing the scaled distance between the coordinate and the ellipse in question
        for i0 in range(len(temp_inds)):
            temp_params = temp_cells[temp_inds[i0]].ellipse_fit[temp_timepoints[i0]]
            if ellipse_contains_point(temp_params, temp_coords[temp_ind1]):
                temp_comparison.append(i0)  # this point has been classified to fall within the current ellipse.

        if len(temp_comparison) > 1:  # if this point fits within more than one ellipse we take the one with the closest
            # centroid
            ind = compare_distances(temp_cells, temp_comparison, temp_coords[temp_ind1], temp_timepoints)
            assignments_made.append(ind)
            temp_cells[temp_inds[ind]].nuclear_whi5[temp_timepoints[ind]] = 1  # this cell has Whi5 nuclear localized.
        else:
            ind = temp_comparison[0]
            temp_cells[temp_inds[ind]].nuclear_whi5[temp_timepoints[ind]] = 1  # this cell has Whi5 nuclear localized.
        assignments_made.append(ind)
    return temp_cells, assignments_made


def automated_whi5_assignment_1(temp_cells, temp_cutoff):  # works based on skewness
    for obj in temp_cells:
        for i0 in range(len(obj.frames)):
            obj.nuclear_whi5[i0] = scipy.stats.skew(obj.zproj_fluor_vals[i0]) > temp_cutoff
    return temp_cells


def automated_whi5_assignment_2(temp_cells, temp_clf):  # works based on machine learning classification
    for obj in temp_cells:
        temp1 = np.array([obj.fluor_chars[temp_ind] for temp_ind in range(len(obj.frames))])
        y_pred = list(temp_clf.predict(temp1))
        if len(y_pred) != len(obj.nuclear_whi5):
            raise ValueError('Error 001. Whi5 classification coded wrongly.')
        else:
            obj.nuclear_whi5 = y_pred
    return temp_cells


def ellipse_contains_point(temp_params, temp_coord):
    # temp_params is the "ellipse_fit" property of cells, with parameters being the outputs from ellipse_fit [:4]
    # i.e. major length, minor length, orientation, centroid (centroid should be x, y). temp_coord is the position
    # being classified.
    # output is a boolean, True if temp_coord sits within the assigned ellipse, false otherwise.

    if len(temp_coord) == 2:
        # 2D ellipse mapping
        temp_centre = np.array([temp_params[3][0], temp_params[3][1]])  # x,y format
        temp_scales = np.array([temp_params[0], temp_params[1]])
        temp_mat = rotation(-temp_params[2])  # inverse rotation about Z axis to determine viable points
        temp_distance = np.linalg.norm((np.dot(temp_mat, (np.asarray(temp_coord) - temp_centre)) / temp_scales))

    return temp_distance < 1  # if this fits within the assigned ellipse, then we return


def compare_distances(temp_cells, temp_inds, temp_coord, temp_timepoints):
    # temp_inds is the index in temp_cells of all cells for which temp_coord fits the ellipse
    temp_centres = [temp_cells[temp_inds[i0]].position[temp_timepoints[i0]] for i0 in range(len(temp_inds))]
    temp_dists = [np.linalg.norm(np.asarray(val - temp_coord)) for val in temp_centres]
    return temp_dists.index(min(temp_dists))  # find the index of the minimum distance


def assign_md_pair(temp_cells, mother_ind, daughter_ind, temp_frame_num):
    # first just check that this division event hasn't already been tracked
    # Note that the input indices here should not be the unique cell identifiers,
    # but rather the index in temp_cells.
    if not (temp_frame_num in temp_cells[mother_ind].daughter_assignment_frames):
        temp_cells[mother_ind].daughters.append(temp_cells[daughter_ind].index)
        temp_cells[daughter_ind].parent = temp_cells[mother_ind].index
        temp_cells[mother_ind].daughter_assignment_frames.append(temp_frame_num)
        temp_cells[daughter_ind].parent_assignment_frame = temp_cells[mother_ind].index
        # enforcing that the nuclear localization state of both cells at this timepoint is correct.
        temp_cells[mother_ind].nuclear_whi5[temp_frame_num - temp_cells[mother_ind].frames[0]] = 1
        temp_cells[daughter_ind].nuclear_whi5[temp_frame_num - temp_cells[daughter_ind].frames[0]] = 1
    return temp_cells


def create_cycles(temp_cells, temp_ind, temp_cycles, ind, temp_data_origin):
    # take a cell and generate cycles for each of its individual cycles, which will be appended to temp_cycles
    temp_cells_indexes = [obj.index for obj in temp_cells]
    temp_cells[temp_ind].temp_cycle_inds = []
    temp_obj = temp_cells[temp_ind]

    # gives the cell that temp_cycles[temp_ind] came from.
    temp = np.diff(np.array(temp_obj.nuclear_whi5))  # 1 gives entering G1, -1 gives leaving G1.
    temp_first_timepoints = np.where(temp == 1)[0]
    # print temp_first_timepoints
    temp_first_timepoints = temp_first_timepoints + 1
    temp_cells[temp_ind].first_timepoints = temp_first_timepoints
    # now this gives the index in the nuclear_whi5 parameter of the appropriate event
    range_points = np.insert(temp_first_timepoints, 0, 0)
    range_points = np.append(range_points, len(temp))  # range_points can now be iterated through to generate cycles of
    # the appropriate length, with the condition that if a cell stopped being tracked just as a new cycle started, this
    # will still work.
    num_cycles = len(temp_first_timepoints) - 1
    # this gives the number of complete cell cycles tracked, ignoring for now
    # whether daughters were assigned at each timepoint.

    for temp_gen in range(len(range_points) - 1):
        range_timepoints = range(range_points[temp_gen], range_points[temp_gen + 1] + 1)
        # gives the indices of the appropriate time points for this cell cycle (relative to the starting tracking time
        # of the cell in question)
        temp_params = {'range': range_timepoints, 'complete': None, 'parent': None, 'error': False,
                       'start': None, 'cc_num': temp_gen, 'cell_index': temp_cells_indexes[temp_ind],
                       'data_origin': temp_data_origin}
        temp_cells[temp_ind].temp_cycle_inds.append(ind)  # keeping track of where this cycle is in the initial cell
        # list

        # defining whether this is a complete cycle or not
        if range_timepoints[0] in temp_first_timepoints and range_timepoints[-1] in temp_first_timepoints:
            temp_params['complete'] = True  # if the start and end points fit the description of temp_first_timepoints
            # this corresponds to a full cell cycle
        else:
            temp_params['complete'] = False  # in this case we missed either the start or the end of the cell cycle.
        # generating the next cell cycles
        #
        # if temp_gen == 0 and not(temp_obj.parent is None):
        #     temp_params['parent'] = temp_obj.parent
        # Updating the cell cycle with the appropriate info
        temp_cycles.append(CellCycle(temp_obj, temp_params))
        ind += 1
    return temp_cells, temp_cycles, ind


def stitch_cycles(temp_cells, temp_cycles, temp_ind):
    # this takes a fully populated list of cell cycles and assigns cell types and lineages based on the lineages in
    # the corresponding list of cells. Assumes that the order of the list is the same as the unique cell cycle index

    # populates the handles "next_gen", "prev_gen", "daughter", "mother", "bud" and "celltype" for all cycles in a given
    # cell. Updates their mothers as well, so that run over the whole population, this should give the full number of
    # annotations.
    temp_obj = temp_cells[temp_ind]
    temp_inds = temp_obj.temp_cycle_inds  # list of locations of cell cycles in temp_cycles associated with cell.
    # NOT unique identifier
    temp_c_inds = [obj.index for obj in temp_cells]
    temp_cc_inds = [obj.index for obj in temp_cycles]
    # tv1 = [obj.index for obj in temp_cycles]
    # tv2 = [temp_cycles[tv1.index(i0)] for i0 in temp_obj.temp_cycle_inds]  # gives a list with the linked cycles from a
    # given cell. Must be length at least 1.
    for i0 in range(len(temp_inds) - 1):
        # assigning next generation. Unique cc identifier.
        temp_cycles[temp_inds[i0]].next_gen = temp_cc_inds[temp_inds[i0 + 1]]
        # assigning parent generation. Unique cc identifier.
        temp_cycles[temp_inds[i0 + 1]].prev_gen = temp_cc_inds[temp_inds[i0]]

    if not (temp_obj.parent is None):  # in this case, the cell in question tracked had its birth tracked
        temp_cycles[temp_inds[0]].celltype = 2  # bud
        # temp_cycles[temp_inds[0]].label_type = temp_cycles
        if len(temp_inds) > 1:
            temp_cycles[temp_inds[1]].celltype = 1  # daughter
            for i1 in range(2, len(temp_inds)):
                temp_cycles[temp_inds[i1]].celltype = 0  # mother

        temp_parent = temp_cells[temp_c_inds.index(temp_obj.parent)]  # gives the parent cell
        # print len(temp_parent.temp_cycle_inds), temp_parent.daughters.index(temp_obj.index)
        # print temp_parent.temp_cycle_inds, temp_parent.daughters
        if (temp_cycles[temp_inds[0]].frames[-1] in temp_parent.daughter_assignment_frames) \
                and (temp_parent.daughters[temp_parent.daughter_assignment_frames.index(
                    temp_cycles[temp_inds[0]].frames[-1])] == temp_obj.index):
            # only assign daughters if we know its parent is recording the correct division event.
            # print temp_cycles[temp_inds[0]].frames[-1], temp_parent.daughter_assignment_frames

            temp_index = temp_parent.daughter_assignment_frames.index(temp_cycles[temp_inds[0]].frames[-1])
            # gives the index at which the parent cell gave birth to this daughter
            if not temp_parent.daughters[temp_index] == temp_obj.index:
                print 'Parent cell'
                study_cell(temp_parent)
                print 'Daughter cell'
                study_cell(temp_obj)
                print 'All daughters of parent'
                for ind1 in temp_parent.daughters:
                    study_cell(temp_cells[temp_c_inds.index(ind1)])
                raise ValueError('Indexing of daughters is incorrect')
            # to give the correct cell cycle index for the parent cell, we take the following:
            temp_index1 = bisect.bisect_left(temp_parent.first_timepoints,
                                             temp_parent.frames.index(temp_cycles[temp_inds[0]].frames[-1]))
            # gives the index of the cell cycle for the parent cell with which this birth is associated.
            temp_parent_cc_index = temp_parent.temp_cycle_inds[temp_index1]  # gives the
            # index of the cell cycle for the parent within the list of cell cycles
            # temp_parent_cc_index = temp_parent.temp_cycle_inds[temp_parent.daughters.index(temp_obj.index)]  # gives the
            # index of the cell cycle for the parent cell

            temp_cycles[temp_inds[0]].parent = temp_cc_inds[temp_parent_cc_index]
            temp_cycles[temp_parent_cc_index].bud = temp_cc_inds[temp_inds[0]]
            # defining the bud for the parent cell cycle
            if len(temp_inds) > 1:
                temp_cycles[temp_parent_cc_index].daughter = temp_cc_inds[temp_inds[1]]
                # defining the daughter for the parent cell
                temp_cycles[temp_inds[1]].mother = temp_cc_inds[temp_parent_cc_index]
                # assigning a handle for the mother cell of the daughter
    elif len(temp_inds) > 1:  # in this case, we never assigned a parent for the current cell, so latter generations are
        # mothers and the first one is "unknown"
        for i1 in range(1, len(temp_inds)):
            temp_cycles[temp_inds[i1]].celltype = 0

    return temp_cells, temp_cycles


def integrate_bud_data(temp_cycles):
    indices = [obj.index for obj in temp_cycles]
    for obj in temp_cycles:
        if not (obj.daughter is None) and obj.complete:  # this cell cycle should have been captured for the full time
            temp = np.insert(np.diff(obj.nuclear_whi5), 0, 0)
            temp1 = np.where(temp == -1)[0]
            # determining the point of start in this cell cycle:
            if len(temp1) == 1:  # testing the number of "start" events in this cell cycle
                obj.start = temp1[0]
            else:
                obj.error = True

            # obj.int_fl_bud = []  # deprecated
            obj.vbud = []
            obj.ellipse_fit_bud = []
            obj.zproj_fl_bud = []
            obj.zproj_fl_bud_c2 = []
            obj.segment_fl_bud = []
            obj.segment_fl_bud_c2 = []
            obj.segment_coords_bud = []
            obj.segment_vol_bud = []
            obj.bud_seg = []

            # Note that by definition we should not incorporate any bud data regarding cytoplasm vs nucleus, because
            # even if a blob is detected within a bud, while it remains classified as a bud, we have not marked a
            # localization event for Whi5. We therefore assume that any "nuclear detection" from a blob detection will
            # be erroneous. As such, all fluorescence within the bud should be calculated based on the
            # segment_fl_bud (both channels).

            for temp_ind in range(len(obj.frames)):
                obj.zproj_fl_bud.append(0)
                obj.vbud.append(0)
                # obj.int_fl_bud.append(0)  # deprecated
                obj.ellipse_fit_bud.append(0)
                obj.bud_seg.append(None)
                obj.zproj_fl_bud_c2.append(0)
                obj.segment_fl_bud.append(0)
                obj.segment_fl_bud_c2.append(0)
                obj.segment_coords_bud.append(0)
                obj.segment_vol_bud.append(0)
            # print temp_cycles[obj.bud].frames, obj.frames
            bud_ind = indices.index(obj.bud)
            for ind in temp_cycles[bud_ind].frames:
                if ind in obj.frames:
                    temp_ind = obj.frames.index(ind)
                    temp_bud_ind = indices.index(obj.bud)
                    obj.vbud[temp_ind] = temp_cycles[temp_bud_ind].ellipse_volume[
                        ind - temp_cycles[temp_bud_ind].frames[0]]
                    # obj.int_fl_bud[temp_ind] = temp_cycles[temp_bud_ind].int_fl[ind-temp_cycles[temp_bud_ind].frames[0]]  # deprecated
                    obj.zproj_fl_bud[temp_ind] = temp_cycles[temp_bud_ind].zproj_fl[
                        ind - temp_cycles[temp_bud_ind].frames[0]]
                    obj.ellipse_fit_bud[temp_ind] = temp_cycles[temp_bud_ind].ellipse_fit[
                        ind - temp_cycles[temp_bud_ind].frames[0]]
                    obj.bud_seg[temp_ind] = temp_cycles[temp_bud_ind].segment_coords[
                        ind - temp_cycles[temp_bud_ind].frames[0]]
                    obj.zproj_fl_bud_c2[temp_ind] = temp_cycles[temp_bud_ind].zproj_fl_c2[
                        ind - temp_cycles[temp_bud_ind].frames[0]]
                    obj.segment_fl_bud[temp_ind] = temp_cycles[temp_bud_ind].pixel_thresh_fluor_vals[
                        ind - temp_cycles[temp_bud_ind].frames[0]]
                    obj.segment_fl_bud_c2[temp_ind] = temp_cycles[temp_bud_ind].pixel_thresh_fluor_vals_c2[
                        ind - temp_cycles[temp_bud_ind].frames[0]]
                    obj.segment_coords_bud[temp_ind] = temp_cycles[temp_bud_ind].pixel_thresh_coords[
                        ind - temp_cycles[temp_bud_ind].frames[0]]
                    obj.segment_vol_bud[temp_ind] = temp_cycles[temp_bud_ind].pixel_thresh_vol[
                        ind - temp_cycles[temp_bud_ind].frames[0]]
                else:
                    # print ind, obj.frames
                    obj.error = True
                    temp_cycles[temp_bud_ind].error = True

    return temp_cycles


def generate_tree(temp_cycles, temp_ind):
    # returns a list with the unique indices of all connected cells for a cell with a given unique index
    temp_inds = [obj.index for obj in temp_cycles]
    temp_tree = [temp_inds[temp_ind]]
    for ind in temp_tree:
        temp_obj = temp_cycles[temp_inds.index(ind)]
        # get all the indices of directly related cells.
        temp_var = [temp_obj.daughter, temp_obj.next_gen, temp_obj.prev_gen, temp_obj.bud, temp_obj.parent]
        # remove None values
        temp_var = [ind1 for ind1 in temp_var if not (ind1 is None)]
        for ind1 in temp_var:
            if not (ind1 in temp_tree):  # only append cell indices if you know it is not already listed. IMPORTANT.
                # Otherwise this will hit a loop and continue indefinitely
                temp_tree.append(ind1)
                # remove any repeats to be extra careful
                # temp_tree = list(set(temp_tree))
    return temp_tree


def correct_bud_labels(temp_cycles, temp_tree):
    temp_type = None
    temp_inds = [obj.index for obj in temp_cycles]
    first_frames = [ind for ind in temp_tree if temp_cycles[temp_inds.index(ind)].frames[0] == 1]
    # check the cells which are present in the first frame.
    if len(first_frames) == 2:
        i0 = 0
        for ind in first_frames:
            if temp_cycles[temp_inds.index(ind)].celltype == 2:  # if this cell is a bud then it may have been
                # incorrectly labeled
                parent_ind = temp_cycles[temp_inds.index(ind)].parent
                if parent_ind == first_frames[i0 - 1]:  # if this first frame has a bud assigned then we use the parent
                    # label type
                    temp1 = temp_cycles[temp_inds.index(parent_ind)].label_type
                    if temp1 != temp_cycles[temp_inds.index(ind)].label_type:
                        print 'Bud label type was reassigned'
                        temp_type = temp1
            i0 += 1
            # temp_cycles[temp_inds.index(ind)].label_type = temp_cycles[parent_ind].label_type
    else:
        print 'Tree of size {0} had {1} roots'.format(len(temp_tree), len(first_frames))
        print temp_cycles[temp_inds.index(first_frames[0])].data_origin
        for ind in temp_tree:
            temp_cycles[temp_inds.index(ind)].error = True  # If we don't understand how the erroneous labeling happened
            # We change these cell cycles to be erroneous so we don't use this data.
    return temp_type


def inherit_lineage_properties(temp_cycles):
    done_cycles = []
    temp_inds = [obj.index for obj in temp_cycles]
    size_fam = []
    i0 = 0
    # families = []
    for ind in range(len(temp_cycles)):
        if not temp_cycles[ind].index in done_cycles:
            temp_tree = generate_tree(temp_cycles, ind)
            if not set(temp_tree).isdisjoint(done_cycles):
                print 'Error: Cells present in multiple lineages'
            temp_types = [temp_cycles[temp_inds.index(ind1)].label_type for ind1 in temp_tree]
            if 1 in temp_types:
                if not (0 in temp_types):  # as long as these lineages are internally consistent
                    temp_type1 = 1
                    for ind1 in temp_tree:
                        temp_cycles[temp_inds.index(ind1)].label_type = temp_type1
                    del temp_type1
                else:
                    print 'Error: Multiple cell labels in a single lineage in family {0}'.format(i0)
                    temp1 = correct_bud_labels(temp_cycles, temp_tree)
                    if not temp1 is None:
                        temp_type1 = temp1
                        for ind1 in temp_tree:
                            temp_cycles[temp_inds.index(ind1)].label_type = temp_type1
                        del temp_type1
            elif 0 in temp_types:
                if not (1 in temp_types):  # as long as these lineages are internally consistent
                    temp_type1 = 0
                    for ind1 in temp_tree:
                        temp_cycles[temp_inds.index(ind1)].label_type = temp_type1

                    del temp_type1
                else:
                    print 'Multiple cell labels in a single lineage in family {0}'.format(i0)
            for ind1 in temp_tree:
                temp_cycles[temp_inds.index(ind1)].family = i0  # indexing which families these cells belong to
            done_cycles += temp_tree
            size_fam.append(len(temp_tree))
            i0 += 1
            # families.append(temp_tree)
    print 'Average family size =', np.mean(size_fam)
    print 'Maximum family size =', np.amax(size_fam)
    print 'Median family size =', np.median(size_fam)
    print 'Number families =', i0
    return temp_cycles


def populate_cells_all_scenes_1(temp_base_path, temp_expt_path, temp_image_filename, temp_bf_filename,
                                temp_num_scenes, temp_num_frames, temp_sel_scenes=None):
    color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
    dims = 512
    for scene in range(1, temp_num_scenes):
        if not (temp_sel_scenes is None):
            scene1 = temp_sel_scenes[scene - 1]
        else:
            scene1 = scene
        # print type(scene1)
        c = []  # this will track all the cells over this timelapse
        # making necessary directory tree
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene1)
        if not os.path.exists(directory):
            os.makedirs(directory)
        if not os.path.exists(directory + '/images'):
            os.makedirs(directory + '/images')
        outlines = np.zeros([temp_num_frames[scene1 - 1], dims, dims])
        for frame_num in range(1, temp_num_frames[scene1 - 1]):
            filename = temp_image_filename + temp_bf_filename + 's{0}_t{1}_segmentation.mat'.format(
                str(scene1), str(frame_num).zfill(2))
            path = temp_base_path + temp_expt_path + '/scene_{0}/segments/'.format(scene1) + filename
            tracking_csv_file = pd.DataFrame.from_csv(
                temp_base_path + temp_expt_path + '/scene_{0}/segments/tracking/tracking.csv'.format(scene1),
                index_col=None)
            c, im, fig, temp_outlines = single_frame(path, c, frame_num, tracking_csv_file, color_seq_val)
            outlines[frame_num - 1, :, :] = temp_outlines[:, :]
            fig.subplots_adjust(bottom=0)
            fig.subplots_adjust(top=1)
            fig.subplots_adjust(right=1)
            fig.subplots_adjust(left=0)
            # extent = mpl.transforms.Bbox(((0, 0), (5, 5)))
            fig.savefig(directory + '/images/frame_{1}.tif'.format(scene1, frame_num))
            del im, fig
        np.save(directory + '/cell_outlines_scene_{0}'.format(scene1), outlines)
        save_object(c, directory + '/cells_scene_{0}.pkl'.format(scene1))
        del c


def populate_cells_all_scenes_2(expt_id):
    # this takes the cell output from script populate_cells_1.py and adds fluorescence data to it
    # temp_fl_filename_c2 is the filename of the secondary fluorescence channel (primary corresponds to Whi5).
    pickle_in = open("./expt_ids" + expt_id + '.pickle', "rb")
    ep = pickle.load(pickle_in)  # this gives us the experimental parameters, so that we can load everything in an
    # orderly fashion
    pixel_size = {'60X': 0.267, '100X': 0.16}
    z_scale, z_offset = 0.4 / pixel_size['60X'], 0

    temp_base_path, temp_expt_path, temp_image_filename, temp_fl_filename = ep['base_path'], ep['expt_path'], ep['image_filename'], ep['fl_filename']
    temp_num_scenes, temp_num_frames, temp_bkgd_scene, temp_fl_filename_c2 = ep['num_scenes'], ep['num_frames'], ep['bkgd_scene'], ep['fl_filename_c2']
    temp_bf_filename = ep['bf_filename']

    for scene in range(1, temp_num_scenes):
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene)
        with open(directory + '/cells_scene_{0}.pkl'.format(scene), 'rb') as input:
            c = pickle.load(input)
        frame_list = [obj.frames for obj in c]
        for obj in c:
            obj.add_fluor_placeholders()
        outlines = np.load(directory + '/cell_outlines_scene_{0}.npy'.format(scene))
        print 'scene = {0}'.format(scene)
        update_list = []
        for frame_num in range(1, temp_num_frames[scene - 1]):
            # loading the brightfield mask to calculate background fluorescence values
            bf_mask = io.imread(temp_base_path+temp_expt_path+'/scene_{0}/segments'.format(scene)+temp_image_filename + temp_bf_filename + 's{0}_t{1}_fgmask.png'.format(
                str(scene), str(frame_num).zfill(2)))>0
            temp = [(frame_num in temp1) for temp1 in frame_list]
            update_list.append(
                [i for i, e in enumerate(temp) if e != 0])  # gives the list of indices that have to be addressed
            # at each frame
            temp_fl_image_c2 = None  # placeholder for the bkgd subtracted C2 fluor image
            # Adding the fluorescence data for the second channel. We do this because the second channel is the one with
            # the strong constitutive fluor.
            if not (temp_fl_filename_c2 is None):  # if we have a second fluorescence channel in this experiment
                filename_fl = temp_image_filename + temp_fl_filename_c2 + 's{0}_t{1}.TIF'.format(str(scene),
                                                                                                 str(frame_num))
                # taking the mean with respect to the z axis. Do it this way since there doesn't seem
                # to be any systematic bias in that direction.
                c, figs, temp_mask, temp_fl_image_c2 = add_fluorescence_traces_c2_v1(
                    temp_path=temp_base_path + temp_expt_path + filename_fl,
                    temp_cells=c, frame_list=update_list[-1], current_frame=frame_num,
                    temp_bf_mask = bf_mask)
                temp_mask = temp_mask.astype('uint16')  # converting it to a format that FIJI can read
                temp_mask *= 65535
                io.imsave(directory + '/images/mask3d_s{0}_t{1}.TIF'.format(str(scene), str(frame_num)), temp_mask)
                # saving the figures in figs. This is only necessary if you want to track how good the segmentation is
                # since it saves a separate image for each cell.
                del figs, temp_mask
            filename_fl = temp_image_filename + temp_fl_filename + 's{0}_t{1}.TIF'.format(str(scene), str(frame_num))

            c, fig = add_fluorescence_traces_v2(temp_path=temp_base_path + temp_expt_path + filename_fl, temp_cells=c,
                                                 frame_list=update_list[-1],
                                                 current_frame=frame_num, z_scaling=z_scale, z_offset=z_offset,
                                                 temp_bf_mask = bf_mask, save_coords=False,
                                                 exists_c2=temp_fl_filename_c2, temp_outlines=outlines,
                                                      temp_im_c2=temp_fl_image_c2)
            del temp_fl_image_c2  # Now the c2 data has been fully integrated into cells we can remove it from memory
            fig.savefig(directory+'/images/spot_detection_s{0}_t{1}.TIF'.format(scene, frame_num))
            print 'done scene {0}, frame {1}'.format(scene, frame_num)
        save_object(c, directory + '/cells_fl_scene_{0}.pkl'.format(scene))
        del c, frame_list


def populate_cells_all_scenes_3(temp_base_path, temp_expt_path, temp_label_path, temp_num_scenes):
    # This script takes the cell output from populate_cells_2 and adds data showing whether cells are labeled as being
    # stained or not. Must have generated csv files with coordinates of labeled cells for this to work.
    temp_frame_num = 1
    for scene in range(1, temp_num_scenes):
        print 'scene = {0}'.format(scene)
        with open(temp_base_path + temp_expt_path + '/scene_{0}/outputs/cells_fl_scene_{0}.pkl'.format(scene),
                  'rb') as input:
            c = pickle.load(input)
        frame_list = [obj.frames for obj in c]
        temp = [(temp_frame_num in temp1) for temp1 in frame_list]
        # gives the indices of each object in the first frame
        update_list = [i for i, e in enumerate(temp) if e != 0]
        tracking_csv_file = pd.DataFrame.from_csv(temp_base_path + temp_label_path + '/scene_{0}.csv'.format(scene),
                                                  index_col=None)
        coords = []
        for i0 in range(len(tracking_csv_file['X'])):
            coords.append(np.array([tracking_csv_file['X'][i0], tracking_csv_file['Y'][i0]]))
        # coords now stores the coordinates of cell center points for the first frame
        c = assign_labels_2(c, update_list, coords)
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene)
        save_object(c, directory + '/cells_fl_lab_scene_{0}.pkl'.format(scene))


def track_localization_manual_annotation(temp_base_path, temp_expt_path, temp_image_filename, temp_fl_filename,
                                         temp_num_frames, temp_analyzed_scene, temp_threshold, temp_drange,
                                         temp_analyzed_frames, temp_label_path):
    # Note: Ensure that track_localization_final.py has run successfully with the correct parameters prior to running
    # this script.
    # this function takes the annotated datasets for Whi5 localization and integrates it into the cell data, so that
    # cells are labeled as containing Whi5 in the correct timepoints. This is then fed into analyze_whi5_distribution.
    # temp_analyzed_scene was previously the scene number that was annotated to give ground truth for Whi5 localization.
    # Now deprecated.
    # temp_num_frames is the number of frames analyzed in temp_analyzed_scene.

    directory = temp_base_path + temp_expt_path + '/whi5_analysis'
    temp_frame_indexes = np.load(directory + '/samples.npy')  # gives the indices of all the analyzed frames.
    temp_indexing = [int(np.sum(temp_num_frames[:i0])) for i0 in range(0, len(temp_num_frames))]
    temp_scenes = [np.searchsorted(temp_indexing, temp_ind, side='right') for temp_ind in temp_frame_indexes]
    # scene numbers
    temp_frames = [temp_frame_indexes[i0] - temp_indexing[temp_scenes[i0] - 1] + 1
                   for i0 in range(len(temp_frame_indexes))]  # frame numbers
    fluor_name = temp_image_filename + temp_fl_filename
    # drange = 65535.0
    c1 = []  # this will be a list with all the cells from the training frames
    y1 = []  # this will be the classification of cell type based on training data
    for i0 in range(len(temp_frame_indexes)):  # iterating through each frame
        dir1 = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(temp_scenes[i0])
        # opening cell file
        if temp_label_path is None:  # in this case we didn't do the labeling since we don't have labels
            temp_loading_path = dir1 + '/cells_fl_scene_{0}.pkl'.format(temp_scenes[i0])
        else:
            temp_loading_path = dir1 + '/cells_fl_lab_scene_{0}.pkl'.format(temp_scenes[i0])
        with open(temp_loading_path, 'rb') as input:
            c = pickle.load(input)
        frame_list = [obj.frames for obj in c]
        # loading cell outlines to confirm assignments
        outlines = np.load(dir1 + '/cell_outlines_scene_{0}.npy'.format(temp_scenes[i0]))
        temp2 = [(temp_frames[i0] in temp1) for temp1 in frame_list]  # checking whether cells are in the current frame
        update_list = [i for i, e in enumerate(temp2) if e != 0]  # indices of all cells that are in the current frame
        # tracked centroids in format x, y
        temp = np.load(
            directory + '/fl_loc_centres/scene_{0}_frame_{1}.npy'.format(temp_scenes[i0], temp_frames[i0]))
        # gives the coordinate list of all annotated cells in the current frame
        coords = zip(temp[:, 0], temp[:, 1])  # format=(x,y)
        # print coords, temp
        # coords now stores the coordinates of cell center points for the analyzed scene, given the current frame
        # num
        c, assignments = assign_labels_3(c, update_list, coords, temp_frames[i0])
        save_object(c, dir1 + '/cells_scene_{0}_v1.pkl'.format(temp_scenes[i0]))
        c1 += [c[temp_ind1].fluor_chars[temp_frames[i0] - c[temp_ind1].frames[0]] for temp_ind1 in update_list]
        y1 += [c[temp_ind1].nuclear_whi5[temp_frames[i0] - c[temp_ind1].frames[0]] for temp_ind1 in update_list]
        temp_new = [ind for ind in range(len(c)) if
                    temp_frames[i0] in c[ind].frames and c[ind].nuclear_whi5[temp_frames[i0] - c[ind].frames[0]]]
        # if i0 == 0:
        #     print temp_new, c[temp_new[0]].nuclear_whi5, temp_scenes[i0], temp_frames[i0], c[temp_new[0]].nuclear_whi5[temp_frames[i0]-c[temp_new[0]].frames[0]]
        # exit()
        # organizing figure data
        temp_im = io.imread(
            temp_base_path + temp_expt_path + fluor_name + 's{0}_t{1}.TIF'.format(temp_scenes[i0], temp_frames[i0]))
        temp_im1 = temp_im / temp_drange
        # print temp_im.shape
        if np.sum(temp_im > temp_threshold) > 0:
            temp_im1 = np.log(np.amax(temp_im1, axis=0) / np.amax(temp_im1))
            # using a log scale in this case because this image contains the frustrating high energy pixels that
            # sometimes crop up
        else:
            temp_im1 = np.amax(temp_im1, axis=0) / np.amax(temp_im1)

        temp_im1 *= outlines[temp_frames[i0] - 1, :, :] == 0
        temp_im1 += outlines[temp_frames[i0] - 1, :, :].astype('uint16')
        # print c[ind].frames
        temp_cell_coords = [c[ind].position[temp_frames[i0] - c[ind].frames[0]] for ind in update_list
                            if c[ind].nuclear_whi5[temp_frames[i0] - c[ind].frames[0]] == 1]
        # store the centroids of the correct cells.

        # plotting figures
        fig = plt.figure(figsize=[10, 10], frameon=False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        ax.imshow(temp_im1)
        temp_x, temp_y = zip(*temp_cell_coords)
        plt.plot(temp_x, temp_y, 'x', color='r', linestyle='None')
        fig.subplots_adjust(bottom=0)
        fig.subplots_adjust(top=1)
        fig.subplots_adjust(right=1)
        fig.subplots_adjust(left=0)
        fig.savefig(directory +
                    '/images/whi5_assignments_scene_{0}_frame_{1}.tif'.format(temp_scenes[i0], temp_frames[i0]))
    np.save(directory + '/training_cells.npy', np.asarray(c1))
    np.save(directory + '/training_cells_labels.npy', np.asarray(y1))


def analyze_whi5_distribution(temp_base_path, temp_expt_path, temp_image_filename, temp_fl_filename, temp_num_frames,
                              temp_analyzed_scene, temp_num_scenes, temp_threshold, temp_drange, temp_analyzed_frames,
                              temp_label_path, temp_sel_scenes=None):
    import seaborn as sns
    # temp_label_path determines whether we use the output from populate_cells_all_scenes_3 or 2.
    temp_dir = temp_base_path + temp_expt_path + '/whi5_analysis'
    temp_frame_indexes = np.load(temp_dir + '/samples.npy')  # gives the indices of all the analyzed frames.
    temp_indexing = [int(np.sum(temp_num_frames[:i0])) for i0 in range(0, len(temp_num_frames))]
    temp_scenes = [np.searchsorted(temp_indexing, temp_ind, side='right') for temp_ind in temp_frame_indexes]
    # scene numbers
    temp_frames = [temp_frame_indexes[i0] - temp_indexing[temp_scenes[i0] - 1] + 1
                   for i0 in range(len(temp_frame_indexes))]  # frame numbers
    if not os.path.exists(temp_dir + '/plots'):
        os.makedirs(temp_dir + '/plots')
    X = np.load(temp_dir + '/training_cells.npy')
    y = np.load(temp_dir + '/training_cells_labels.npy')
    lower_bound = 0.01
    inds1 = np.nonzero(y)
    inds2 = np.nonzero(y == 0)
    for i0 in range(len(fluor_chars_names)):
        v1 = X[inds1, i0]
        v2 = X[inds2, i0]
        print v1.shape
        make_dist_plots_1 = 1
        if make_dist_plots_1:
            fig = plt.figure(figsize=[5, 5])
            sns.distplot(v1, label='Nuclear')
            sns.distplot(v2, label='Cytoplasmic')
            plt.legend()
            plt.xlabel(fluor_chars_names[i0])
            plt.title(fluor_chars_names[i0] + ' in manually annotated cell populations')
            fig.savefig(temp_dir + '/plots/' + fluor_chars_names[i0])
            del fig

    fig = plt.figure(figsize=[5, 5])
    kde = scipy.stats.gaussian_kde(v1)
    kde1 = scipy.stats.gaussian_kde(v2)
    xpoints = np.linspace(-2.0, 5.0, 100)
    vals = [kde.integrate_box_1d(low=-10, high=lim) for lim in xpoints]
    sel = [i for i, e in enumerate(vals) if e > lower_bound]
    skew_cutoff = xpoints[sel[0]]
    plt.axvline(x=xpoints[sel[0]], label='99% bound = {0}'.format(str(np.round(skew_cutoff, decimals=3))))
    print kde1.integrate_box_1d(low=skew_cutoff, high=10.0)
    plt.plot(xpoints, kde.evaluate(xpoints), label='Nuclear')
    plt.plot(xpoints, kde1.evaluate(xpoints), label='Cytoplasmic')
    plt.plot(xpoints, vals, label='Nuclear CDF')
    plt.legend()
    fig.savefig(temp_dir + '/plots/skewness_bounding')

    temp_clf = SVC(kernel="linear", C=0.025)  # linear svm for classification
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5, random_state=42)
    temp_clf.fit(X_train, y_train)
    score = temp_clf.score(X_test, y_test)
    print "Learning algorithm score on test data:", score
    y_test_pred = temp_clf.predict(X_test)
    out = summary_statistics(y_test_pred, y_test)
    print 'Accuracy, Sensitivity, Specificity'
    print np.around(out, 4)

    for scene1 in range(1, temp_num_scenes):
        if not (temp_sel_scenes is None):  # if we have only selected a certain number of scenes to look at here.
            scene = temp_sel_scenes[scene1 - 1]
        else:
            scene = scene1
        print 'Working on Scene number {0}'.format(scene)
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene)
        if temp_label_path is None:  # in this case we didn't do the labeling since we don't have labels
            temp_loading_path = directory + '/cells_fl_scene_{0}.pkl'.format(scene)
        else:
            temp_loading_path = directory + '/cells_fl_lab_scene_{0}.pkl'.format(scene)
        with open(temp_loading_path, 'rb') as input:
            c = pickle.load(input)
        c = automated_whi5_assignment_2(c, temp_clf)  # ML classifier
        # c = automated_whi5_assignment_1(c, skew_cutoff)  # skewness cutoff classifier
        save_object(c, directory + '/cells_scene_{0}_v2.pkl'.format(scene))
        # producing plots to ensure that this tracking is adequate
        outlines = np.load(directory + '/cell_outlines_scene_{0}.npy'.format(scene))
        frame_list = [obj.frames for obj in c]

        # organizing figure data
        for frame_num in range(1, temp_num_frames[scene - 1]):
            # checking which cells to look at for this timepoint
            temp2 = [(frame_num in temp1) for temp1 in frame_list]
            update_list = [i for i, e in enumerate(temp2) if e != 0]
            # loading image
            temp_im = io.imread(
                temp_base_path + temp_expt_path + temp_image_filename + temp_fl_filename + 's{0}_t{1}.TIF'.format(
                    str(scene),
                    str(frame_num)))
            temp_im1 = temp_im / temp_drange
            if np.sum(temp_im > temp_threshold) > 0:
                temp_im1 = np.log(np.amax(temp_im1, axis=0) / np.amax(temp_im1))
                # using a log scale in this case because this image contains the frustrating high energy pixels that
                # sometimes crop up
            else:
                temp_im1 = np.amax(temp_im1, axis=0) / np.amax(temp_im1)
            temp_im1 *= outlines[frame_num - 1, :, :] == 0
            temp_im1 += outlines[frame_num - 1, :, :].astype('uint16')
            temp_cell_coords = [c[ind1].position[frame_num - c[ind1].frames[0]] for ind1 in update_list if
                                (c[ind1].nuclear_whi5[
                                     frame_num - c[ind1].frames[0]])]  # store the centroids of the correct cells.
            # plotting figures
            fig = plt.figure(figsize=[10.0, 10.0], frameon=False)
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)
            ax.imshow(temp_im1, cmap='viridis')
            if len(temp_cell_coords) > 0:
                temp_x, temp_y = zip(*temp_cell_coords)
                plt.plot(temp_x, temp_y, '.', markersize=4, color='r', linestyle='None')
            else:
                print 'No cells in G1 during frame {0}'.format(frame_num)
            fig.subplots_adjust(bottom=0)
            fig.subplots_adjust(top=1)
            fig.subplots_adjust(right=1)
            fig.subplots_adjust(left=0)
            fig.savefig(directory + '/images/automated_whi5_assignments_frame_{0}.tif'.format(frame_num))
    for i0 in range(len(temp_frame_indexes)):
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(temp_scenes[i0])
        name = '/automated_whi5_assignments_frame_{0}.tif'.format(temp_frames[i0])
        copyfile(directory + '/images' + name,
                 temp_dir + '/images' + '/whi5_assignments_scene_{0}_frame_{1}_auto.tif'.format(temp_scenes[i0],
                                                                                                temp_frames[i0]))
        # so that we can compare them directly


def assign_troublesome_pair(temp_posn, temp_frame_num, temp_dir, temp_cell_num, temp_scene):
    # Unfinished, but this would prompt the user to find the position of the partner cell at that timepoint and add it to a .txt tab
    # delimited file for this scene and frame number. The format should be
    # cell  X   Y
    # num   x1  y1
    print 'Assign pair cell for scene {0}, frame {1}, cell number {2}'.format(temp_scene, temp_frame_num,
                                                                              temp_cell_num)
    temp_path = temp_dir + '/cell_coords/scene_{0}_frame_{1}.txt'.format(temp_scene, temp_frame_num)
    print 'position: ', temp_posn
    if not os.path.exists(temp_dir + '/cell_coords'):
        os.makedirs(temp_dir + '/cell_coords')
    if not os.path.exists(temp_path):
        f = open(temp_path, 'w')
        f.write('cell\tX\tY')
        f.close()
        print 'Current cell partner storage:'
    with open(temp_path, 'r') as f:
        print f.read()
    temp_complete = input('Cell partner assigned: ')
    if temp_complete == 'y':
        temp_posns_pd = pd.read_csv(
            temp_dir + '/cell_coords/scene_{0}_frame_{1}.txt'.format(temp_scene, temp_frame_num),
            sep='\t')
        temp1 = temp_posns_pd.cell == temp_cell_num
        temp_posn_partner = np.array([temp_posns_pd[temp1]['X'].iloc[0], temp_posns_pd[temp1]['Y'].iloc[0]])
    else:
        temp_posn_partner = None
    return temp_posn_partner


def determine_mothers(temp_cells, temp_frame_num, temp_inds):
    # This will take two cells, and assign mother vs daughter, preferentially using previously known lineage data but
    # otherwise using size data. If both cells are known to be mother cells, we return an error.
    temp_inds1 = None
    temp_mother_vec = np.array([(len(temp_cells[temp_ind].daughters) > 0) for temp_ind in temp_inds])
    if np.sum(temp_mother_vec) == 1:
        temp_inds2 = [i for i, x in enumerate(temp_mother_vec) if x]
        temp_inds1 = [temp_inds[temp_inds2[0]], temp_inds[temp_inds2[0] - 1]]  # gives the index of mother first, then
        # daughter
    elif np.sum(temp_mother_vec) == 0:  # if no cells are mothers we report the mother as the larger of the two cells.
        temp_val = [temp_cells[temp_inds[0]].ellipse_volume[temp_frame_num - temp_cells[temp_inds[0]].frames[0]],
                    temp_cells[temp_inds[1]].ellipse_volume[temp_frame_num - temp_cells[temp_inds[1]].frames[0]]]
        # print temp_val
        m_ind = np.argmax(temp_val)
        temp_inds1 = [temp_inds[m_ind], temp_inds[m_ind - 1]]  # gives the index of mother first, then
        # daughter
    if temp_inds1 is None:
        print 'Tried to assign a mother-mother pairing at division in frame {0}'.format(temp_frame_num)
    return temp_inds1


def test_tree(temp_conn_mat, temp_assigned, temp_ind):
    # classifies the branch and classification structure. temp_type = 0 must be added to pending_assignment, temp_type=1
    # must be assigned to its pair, and temp_type=2 requires manual assignment
    temp_sum = np.sum(temp_conn_mat[temp_ind, :], axis=0)
    if temp_sum == 0:
        temp_type = 0  # if there is no prospective pair
    elif temp_sum == 1:
        temp_ind1 = np.nonzero(temp_conn_mat[temp_ind, :])[0][0]
        if not (temp_assigned[temp_ind1]):  # if the prospective pair cell has not been assigned yet
            if np.sum(temp_conn_mat[temp_ind1, :]) == 1:
                temp_type = 1  # if this is a regular tree with two pairs that are unassigned.
            elif np.sum(temp_assigned[np.nonzero(temp_conn_mat[temp_ind1, :])] == 0) == 1:  # if there is still only one
                # unassigned pair for the pair then we treat this as normal.
                temp_type = 1
            else:
                temp_type = 2  # if the other cell has more than one pair, and more than one is unassigned, then this
                # needs to be manually assigned
        else:  # if the prospective pair has already been assigned to another cell
            temp_type = 0
    elif temp_sum > 1:
        if np.sum(temp_assigned[np.nonzero(temp_conn_mat[temp_ind, :])] == 0) == 1:  # if there are multiple potential
            # pairs but only one is unassigned
            temp_type = -1  # in this case we need to select the pair as unassigned one.
        elif np.sum(temp_assigned[np.nonzero(temp_conn_mat[temp_ind, :])] == 0) > 1:
            temp_type = 2  # if there is more than one prospective pair and more than one unassigned then this requires
            # manual assignment
        else:
            temp_type = 0  # if there are no unassigned pairs then this should be flagged for appending to
            # pending_assignment
    return temp_type


def assign_lineages(temp_base_path, temp_expt_path, temp_image_filename, temp_fl_filename, temp_num_frames,
                    temp_analyzed_scene, temp_num_scenes, temp_threshold, temp_drange,
                    temp_manual_annotation):
    # This function takes as inputs the current cells with specified Whi5 localization, and assigns lineages
    # temp_manual_annotation is a boolean value for whether this will involve you manually annotating these datasets.
    distance_cutoff = 40  # cutoff for centroid distance above which mother-daughter pairings cannot be assigned
    num_divisions = 0
    num_nuclear_events = 0
    for scene in range(1, temp_num_scenes):
        print 'Starting scene', scene
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene)
        with open(directory + '/cells_scene_{0}_v2.pkl'.format(scene), 'rb') as input:
            c = pickle.load(input)
        # cell_indices = [obj.index for obj in c]
        for obj in c:
            obj.daughters = []
            obj.parent = None
            obj.daughter_assignment_frames = []
            obj.parent_assignment_frame = None
        frame_list = [obj.frames for obj in c]
        pending_assignment = []
        outlines = np.load(directory + '/cell_outlines_scene_{0}.npy'.format(scene))
        for frame_num in range(1, temp_num_frames[
                    scene - 1]):  # we go through frame by frame, starting with frame # 2 since we need a previous
            # timepoint to compare it to.
            if frame_num % 10 == 0:
                print 'Reached Scene {0}, frame {1}'.format(scene, frame_num)

            # Now we check which cells to analyze for each frame
            assigned_inds = [[], []]  # this will be used later to generate a figure showing the mother and daughter
            # pairs for each cell type
            pending_assignment.append([])
            if frame_num > 1:
                temp2 = [(frame_num in temp1) for temp1 in frame_list]
                update_list = [i for i, e in enumerate(temp2) if e != 0]  # list of cells in the current frame
                # we now consider only cells for whom Whi5 is in the nucleus. This happens coincidentally for
                # mother and daughter cells, and we use spatial proximity to determine which daughter belongs to which
                # mother.
                update_list1 = [i for i in update_list if c[i].nuclear_whi5[frame_num - c[i].frames[0]] == 1]
                temp_new = [[], [], []]
                temp_started_G1 = [[], [], []]
                for ind1 in update_list1:
                    temp_ind = frame_num - c[ind1].frames[0]  # cell age
                    if temp_ind == 0:  # if this cell was segmented for the first time in this frame then we don't use
                        # that since erroneous classification happens quite a lot for very small cells due to poor
                        # segmentation.
                        print 'Newborn cell {0} with nuclear Whi5 found in timepoint {1}'.format(c[ind1].index,
                                                                                                 frame_num)
                        temp_new[0].append(ind1)
                        temp_new[1].append(temp_ind)
                        # print c[ind1].position
                        temp_new[2].append(c[ind1].position[temp_ind])
                        num_nuclear_events += 1
                    elif c[ind1].nuclear_whi5[temp_ind - 1] == 0:
                        # if temp_ind-2==0 or c[ind1].nuclear_whi5[temp_ind-2] == 0:  # we require that either the cell only
                        # started being
                        temp_started_G1[0].append(ind1)  # these are the cells we have to assign. Note this is the
                        # index in c, NOT the unique cell identifier.
                        temp_started_G1[1].append(temp_ind)
                        temp_started_G1[2].append(c[ind1].position[temp_ind])
                        num_nuclear_events += 1
                # NOTE:
                # temp_started_G1 now contains details about those cells present in the current frame, with current Whi5
                # localization that did not have Whi5 localized in the nucleus in the previous frame. The details stored
                # are: the index in c, age at current frame, position at current frame.
                temp_mothers = np.array([(len(c[temp_ind1].daughters) > 0) for temp_ind1 in temp_started_G1[0]])
                # tracking whether these cells have already been mothers.
                to_remove = [[], []]
                assigned = np.zeros(len(temp_started_G1[0]))  # This will track whether a cell has already been assigned
                # We use this to avoid multiple assignments
                dist_mat = np.zeros(
                    [len(temp_started_G1[0]), len(temp_started_G1[0])])  # matrix of distances between cell
                # i and j for dist_mat[i,j]
                dist_mat += (distance_cutoff + 1000.0) * np.identity(
                    len(temp_started_G1[0]))  # prevent cells partnering
                # themselves.
                for ind1 in range(len(temp_started_G1[0])):
                    for ind2 in range(len(temp_started_G1[0])):
                        if ind1 != ind2:  # prevent cells partnering themselves.
                            dist_mat[ind1, ind2] = np.linalg.norm(temp_started_G1[2][ind1] - temp_started_G1[2][ind2])
                conn_mat = dist_mat < distance_cutoff
                # We now go through each cell in temp_started_G1 and see what other cell it is closest to.
                for ind1 in range(len(temp_started_G1[0])):
                    if not (assigned[ind1]):  # if the current cell has not already been assigned
                        # evaluating the number of candidate partners this cell has
                        is_pair = test_tree(dist_mat < distance_cutoff, assigned, ind1)
                        if is_pair == 1:  # this is the regular case
                            temp_ind1 = np.nonzero(conn_mat[ind1, :])[0][0]
                            inds1 = determine_mothers(c, frame_num,
                                                      [temp_started_G1[0][ind1], temp_started_G1[0][temp_ind1]])
                            if not (inds1 is None):  # if there wasn't a problem in assigning mother vs daughter
                                assigned_inds[0].append(inds1[0])
                                assigned_inds[1].append(inds1[1])
                                c = assign_md_pair(c, mother_ind=inds1[0],
                                                   daughter_ind=inds1[1],
                                                   temp_frame_num=frame_num)
                                assigned[ind1] = 1
                                assigned[temp_ind1] = 1
                                num_divisions += 1
                        elif is_pair == 2 and temp_manual_annotation:
                            # If the pairing is ambiguous and this cell has not already been assigned then we have to
                            # manually assign this.
                            temp_loc = assign_troublesome_pair(temp_started_G1[2][ind1], frame_num,
                                                               directory, temp_started_G1[0][ind1], scene)
                            if not (temp_loc is None):  # if we were able to find a candidate
                                # first we search for this within the list of cells for which Whi5 has just entered the
                                # nucleus
                                d1 = np.ones(len(temp_started_G1[0]))
                                for ind2 in range(len(temp_started_G1[0])):
                                    d1[ind2] = np.linalg.norm(
                                        temp_started_G1[2][ind2] - temp_loc)
                                if np.amin(d1) < 20 and conn_mat[ind1, np.argmin(d1)]:
                                    # if we can assign the point to a cell that has newly started G1 and this selected
                                    # cell is sufficiently close to our cell in question
                                    inds1 = determine_mothers(c, frame_num,
                                                              [temp_started_G1[0][ind1],
                                                               temp_started_G1[0][np.argmin(d1)]])
                                    if not inds1 is None:  # if there was no problem with assignment mother vs. daughter
                                        assigned_inds[0].append(inds1[0])
                                        assigned_inds[1].append(inds1[1])
                                        c = assign_md_pair(c, mother_ind=inds1[0],
                                                           daughter_ind=inds1[1],
                                                           temp_frame_num=frame_num)
                                        # tracking that these cells have been assigned
                                        assigned[ind1] = 1
                                        assigned[np.argmin(d1)] = 1
                                        num_divisions += 1
                            else:  # if we can't find a pair in the newly G1 population we simply add the current
                                # cell to pending_assignment
                                print 'Unable to assign a pair in scene {0}, frame {1}, cell {2}'.format(
                                    scene, frame_num, temp_started_G1[0][ind1])
                                pending_assignment[-1].append(
                                    temp_started_G1[0][ind1])  # if we can't assign this well, we compare it to
                                # the previous timepoint and the following one.
                        elif is_pair == -1:  # in this case, the cell has multiple potential partners, but only one is
                            # unassigned. We therefore can assign the cell the unassigned potential partner.
                            temp_num = np.nonzero(
                                assigned[np.nonzero(conn_mat[ind1, :])] == 0)  # the index in nonzero ent
                            # of conn_mat in which the cell has not been assigned yet
                            ind2 = np.nonzero(conn_mat[ind1, :])[0][temp_num][0]  # gives the index in "assigned" of the
                            # pair
                            inds1 = determine_mothers(c, frame_num,
                                                      [temp_started_G1[0][ind1], temp_started_G1[0][ind2]])
                            if not inds1 is None:  # if there was no problem with assignment mother vs. daughter
                                assigned_inds[0].append(inds1[0])
                                assigned_inds[1].append(inds1[1])
                                c = assign_md_pair(c, mother_ind=inds1[0],
                                                   daughter_ind=inds1[1],
                                                   temp_frame_num=frame_num)
                                # tracking that these cells have been assigned
                                assigned[ind1] = 1
                                assigned[ind2] = 1
                                num_divisions += 1
                        elif is_pair == 0:
                            pending_assignment[-1].append(
                                temp_started_G1[0][ind1])  # if we can't assign this well, we compare it to
                            # the previous timepoint and the following one
                            print 'Unable to assign a pair in scene {0}, frame {1}, cell {2}'.format(
                                scene, frame_num, temp_started_G1[0][ind1])
                # we take any cell that is pending assignment in the current frame, and compare to the previous frame to
                # see if there is an appropriate assignment there.
                dist_mat = (distance_cutoff + 1000.0) * np.ones(
                    [len(pending_assignment[-1]), len(pending_assignment[-2])])
                # note that here we should never have the same cell appear in subsequent timesteps since we only append
                # cells to this vector that are tracked as having had Whi5 not in the nucleus in the nucleus in the prev
                # timestep. I.e. if it is in pending_assignment[-2] it cannot be in pending_assignment[-1].
                for i0 in range(len(pending_assignment[-1])):
                    coord1 = c[pending_assignment[-1][i0]].position[
                        frame_num - c[pending_assignment[-1][i0]].frames[0]]
                    for i1 in range(len(pending_assignment[-2])):
                        if frame_num in c[pending_assignment[-2][i1]].frames:  # if the cell in the previous timept
                            # still exists in the current one we can assign it
                            # compare the positions in the current frame and the previous frame
                            coord2 = c[pending_assignment[-2][i1]].position[frame_num - 1 -
                                                                            c[pending_assignment[-2][i1]].frames[0]]
                            dist_mat[i0, i1] = np.linalg.norm(coord1 - coord2)
                conn_mat = dist_mat < distance_cutoff
                assigned = [np.zeros(conn_mat.shape[0]), np.zeros(conn_mat.shape[1])]
                for i0 in range(len(pending_assignment[-1])):
                    if not (assigned[0][i0]):  # if the cell has not already been assigned
                        if np.sum(conn_mat[i0, :]) == 1:  # if the current cell has only one potential pair
                            i1 = np.nonzero(conn_mat[i0, :])[0][0]  # index of the paired cell in pending_assignment[-2]
                            if np.sum(conn_mat[:, i1]) == 1 and not (assigned[1][i1]):  # if the potential pair has only
                                # one potential pair and is not already assigned (similar to is_pair==1 above).
                                inds1 = determine_mothers(c, frame_num,
                                                          [pending_assignment[-1][i0], pending_assignment[-2][i1]])
                                if not inds1 is None:  # if there was no problem with assignment mother vs. daughter
                                    assigned_inds[0].append(inds1[0])
                                    assigned_inds[1].append(inds1[1])
                                    c = assign_md_pair(c, mother_ind=inds1[0],
                                                       daughter_ind=inds1[1],
                                                       temp_frame_num=frame_num - 1)
                                    # we assign these cells paired in the previous frame, since that's when
                                    # the nuclear localization event must have occurred for both. However, they will appear
                                    # in the output images as being tracked in the current frame.
                                    assigned[0][i0] = 1
                                    assigned[1][i1] = 1
                                    num_divisions += 2
                                    to_remove[0].append(pending_assignment[-1][i0])
                                    to_remove[1].append(pending_assignment[-2][i1])
                                    print 'Assigned M {0}, D {1} in frame {2}'.format(inds1[0], inds1[1], frame_num)
                                    # removing the newly tracked items from pending_assignment. Cells from the previous
                                    # timestep will still appear as not being tracked in the previous frame.
                # removing the assigned cells from the pending_assignment list
                pending_assignment[-1] = list(set(pending_assignment[-1]) - set(to_remove[0]))
                pending_assignment[-2] = list(set(pending_assignment[-2]) - set(to_remove[1]))

            ######################## Producing a figure to track the results so far
            # producing a figure to track the result of this so far.
            # Note that the "bad cells" shown in red will only be listed for the cells that have not been assigned
            # retroactively
            temp_im = io.imread(
                temp_base_path + temp_expt_path + temp_image_filename + temp_fl_filename +
                's{0}_t{1}.TIF'.format(str(scene), str(frame_num))) / temp_drange
            temp_im1 = np.amax(temp_im, axis=0) / np.amax(temp_im)  # scaling for visualization purposes.
            temp_im1 *= outlines[frame_num - 1, :, :] == 0
            temp_im1 += outlines[frame_num - 1, :, :].astype('uint16')
            m_cell_coords = [c[ind1].position[frame_num - c[ind1].frames[0]] for ind1 in assigned_inds[0]]
            # for ind1 in assigned_inds[1]:
            #     print frame_num, c[ind1].frames[0], len(c[ind1].frames), c[ind1].frames
            d_cell_coords = [c[ind1].position[frame_num - c[ind1].frames[0]] for ind1 in assigned_inds[1]]
            bad_cell_coords = [c[ind1].position[frame_num - c[ind1].frames[0]] for ind1 in pending_assignment[-1]]
            # store the centroids of the correct cells.
            # print temp_cell_coords

            # print temp_cell_coords
            # print temp_cell_coords

            # plotting figures
            fig = plt.figure(figsize=[5.12, 5.12], frameon=False)
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)
            ax.imshow(temp_im1, cmap='viridis')
            if len(m_cell_coords) > 0:
                temp_x, temp_y = zip(*m_cell_coords)
                plt.plot(temp_x, temp_y, '.', markersize=6, color='w', linestyle='None')
            if len(d_cell_coords) > 0:
                temp_x, temp_y = zip(*d_cell_coords)
                plt.plot(temp_x, temp_y, '.', markersize=4, color='k', linestyle='None')
            if len(bad_cell_coords) > 0:
                temp_x, temp_y = zip(*bad_cell_coords)
                plt.plot(temp_x, temp_y, '.', markersize=2, color='r', linestyle='None')
            fig.subplots_adjust(bottom=0)
            fig.subplots_adjust(top=1)
            fig.subplots_adjust(right=1)
            fig.subplots_adjust(left=0)
            fig.savefig(directory + '/images/automated_lineage_assignments_frame_{0}.tif'.format(frame_num))
        save_object(c, directory + '/cells_scene_{0}_v3.pkl'.format(scene))
        # print 'I got here', directory + '/images/cells_scene_{0}_v3.pkl'.format(scene)


def create_cycles_full(temp_base_path, temp_expt_path, temp_num_scenes):
    # This function takes the output from assign_lineages and creates individual cell cycles from that
    cc = []
    ind = 0  # this indexes the cell cycles we are generating.
    # this will be repeated with each new iteration

    for scene in range(1, temp_num_scenes):
        print 'Scene Number: {0}'.format(scene)
        data_index = [temp_base_path + temp_expt_path, scene]
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene)
        with open(directory + '/cells_scene_{0}_v3.pkl'.format(scene), 'rb') as input:
            c = pickle.load(input)
        for i0 in range(len(c)):
            # if not(i0 == c[i0].index-1):
            #     print 'Cell indexing error in list number {0}, index number {1}'.format(i0, c[i0].index)
            #     raise ValueError('Cell indexing distinct from cell list indexing')
            # else:
            c, cc, ind = create_cycles(c, i0, cc, ind,
                                       temp_data_origin=data_index)  # note that the unique index for each
            # cell cycle is unique across scenes also
        for i0 in range(len(c)):
            c, cc = stitch_cycles(c, cc, i0)
        print len(cc)
        # cc += cc
    cc = inherit_lineage_properties(cc)  # inheriting the lineage properties of related cells
    print len(cc), np.sum([1 for obj in cc if obj.complete and not (obj.daughter is None)])
    print np.sum([1 for obj in cc if obj.complete and not (obj.daughter is None) and obj.celltype == 1])
    print np.sum(
        [1 for obj in cc if obj.complete and not (obj.daughter is None) and obj.celltype == 1 and not obj.error])
    print np.sum([1 for obj in cc if obj.complete and not (obj.daughter is None) and obj.celltype == 0])
    print np.sum(
        [1 for obj in cc if obj.complete and not (obj.daughter is None) and obj.celltype == 0 and not obj.error])
    cc = integrate_bud_data(cc)
    print len(cc), np.sum([1 for obj in cc if obj.complete and not (obj.daughter is None)])
    print np.sum([1 for obj in cc if obj.complete and not (obj.daughter is None) and obj.celltype == 1])
    print np.sum(
        [1 for obj in cc if obj.complete and not (obj.daughter is None) and obj.celltype == 1 and not obj.error])
    print np.sum([1 for obj in cc if obj.complete and not (obj.daughter is None) and obj.celltype == 0])
    print np.sum(
        [1 for obj in cc if obj.complete and not (obj.daughter is None) and obj.celltype == 0 and not obj.error])

    save_object(cc, temp_base_path + temp_expt_path + '/cell_cycles_compiled.pkl'.format(scene))


def gini(array):  # https://github.com/oliviaguest/gini/blob/master/gini.py, accessed 02/7/2019
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq:
    # http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from:
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    # All values are treated equally, arrays must be 1d:
    array = array.flatten()
    if np.amin(array) < 0:
        # Values cannot be negative:
        array -= np.amin(array)
    # Values cannot be 0:
    array += 0.0000001
    # Values must be sorted:
    array = np.sort(array)
    # Index per array element:
    index = np.arange(1, array.shape[0] + 1)
    # Number of array elements:
    n = array.shape[0]
    # Gini coefficient:
    return (np.sum((2 * index - n - 1) * array)) / (n * np.sum(array))


def summary_statistics(test_out, test_y):
    num_fp = np.sum((test_out - test_y) > 0)
    num_tp = np.sum((test_out * test_y) > 0)
    num_tn = np.sum((test_out - 1) * (test_y - 1) > 0)
    num_fn = np.sum((test_out - test_y) < 0)
    acc = (num_tp + num_tn) * 1.0 / (num_fp + num_tp + num_tn + num_fn)
    sens = num_tp * 1.0 / (num_tp + num_fn)
    spec = num_tn * 1.0 / (num_tn + num_fp)
    print 'FP', 'TP', 'TN', 'FN'
    print num_fp, num_tp, num_tn, num_fn
    return acc, sens, spec


def study_cell(temp_cell):
    print 'Cell number: ', temp_cell.index
    print 'Index image', temp_cell.index_image
    print 'Frames', temp_cell.frames
    print 'position', temp_cell.position
    print 'nuclear whi5', temp_cell.nuclear_whi5
    print 'type', temp_cell.type
    print 'area', temp_cell.area
    print 'daughter assignment frames', temp_cell.daughter_assignment_frames
    print 'daughters', temp_cell.daughters
    print 'parent', temp_cell.parent_assignment_frame


def validate_cycles(temp_base_path, temp_expt_path, temp_image_filename, temp_bf_filename, temp_num_frames,
                    temp_num_scenes):
    with open(temp_base_path + temp_expt_path + '/cell_cycles_filtered.pkl', 'rb') as input:
        temp_cycles = pickle.load(input)
    filt_cc = [obj for obj in temp_cycles if obj.complete and not (obj.error) and not (obj.daughter is None)]
    for scene in range(1, temp_num_scenes):
        scene_cycles = [obj for obj in filt_cc if obj.data_origin[1] == scene]  # selecting the cycles in this scene
        frame_list = [obj.frames for obj in scene_cycles]
        for frame_num in range(1, temp_num_frames[scene - 1]):
            update_list = [i for i, e in enumerate(frame_list) if frame_num in e]  # gives the number of objects in the
            # current frame
            # load the current image
            filename = temp_image_filename + temp_bf_filename + 's{0}_t{1}.TIF'.format(
                str(scene), str(frame_num).zfill(2))
            temp_im = io.imread(temp_base_path + temp_expt_path + '/scene_{0}'.format(scene) + filename) / 65535.0
            # convert this grayscale to rgb
            temp_im1 = np.repeat(temp_im[:, :, np.newaxis], 3, axis=2)
            for temp_ind in update_list:  # going through to add segmented visualization
                temp_coords = zip(*scene_cycles[temp_ind].segment_coords[frame_num - scene_cycles[temp_ind].frames[0]])
                temp_mask_im = np.zeros(temp_im.shape)
                temp_mask_im[temp_coords] = 1  # creating the mask for the main cell
                if not (scene_cycles[temp_ind].bud_seg[frame_num - scene_cycles[temp_ind].frames[0]] is None):
                    # if we have coords for the bud at this stage we add that data
                    temp_coords = zip(*scene_cycles[temp_ind].bud_seg[frame_num - scene_cycles[temp_ind].frames[0]])
                    temp_mask_im[temp_coords] = 1  # creating the mask for the bud
                temp_mask_im1 = np.repeat(temp_mask_im[:, :, np.newaxis], 3, axis=2)  # now in rgb shape!
                temp_rgb_vals = cm.tab20(scene_cycles[temp_ind].index % 20)  # cyclic color map so that we track cell
                # cycles
                temp_mask_im2 = temp_mask_im1 == 0
                temp_im1 *= temp_mask_im2
                temp_im1 += temp_mask_im1 * np.asarray(temp_rgb_vals[:3])
            fig = plt.figure(figsize=[5.12, 5.12],
                             frameon=False)  # note figsize is selected to scale with that of the image
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.imshow(temp_im1)
            ax.set_axis_off()
            fig.add_axes(ax)
            fig.subplots_adjust(bottom=0)
            fig.subplots_adjust(top=1)
            fig.subplots_adjust(right=1)
            fig.subplots_adjust(left=0)
            # extent = mpl.transforms.Bbox(((0, 0), (5, 5)))
            fig.savefig(temp_base_path + temp_expt_path + '/scene_{0}'.format(scene) + '/outputs' +
                        '/images/cell_cycle_validation_frame_{1}.tif'.format(scene, frame_num))
            del fig


def label_boundaries(temp_base_path, temp_expt_path, temp_image_filename, temp_bf_filename, temp_num_scenes,
                     temp_num_frames, temp_sel_scenes=None):
    cells = []
    bdys = set([0, 511])
    for scene1 in range(1, temp_num_scenes):
        if not (temp_sel_scenes is None):
            scene = temp_sel_scenes[scene1 - 1]
        else:
            scene = scene1
        with open(temp_base_path + temp_expt_path + '/scene_{0}/outputs/cells_fl_scene_{0}.pkl'.format(scene),
                  'rb') as input:
            temp_cells = pickle.load(input)
        for obj in temp_cells:
            obj.edge_cell = False  # default setting
            obj.scene_num = scene  # track which scene a cell comes from
            for i0 in range(len(obj.frames)):
                temp_x, temp_y = zip(*obj.segment_coords[i0])
                truval = 0
                if bdys.intersection(set(temp_x)):  # if this cell intersects the image boundaries
                    truval = 1
                elif bdys.intersection(set(temp_y)):  # if this cell intersects the image boundaries
                    truval = 1
                if truval:
                    obj.edge_cell = True  # record this for the cell
        cells += temp_cells
    save_object(cells, temp_base_path + temp_expt_path + '/cells_compiled.pkl')


def filter_cycles(temp_base_path, temp_expt_path, temp_scale, temp_size_thresh=None):
    # temp_height is 3.5um for haploid cell_asics plates, 5.0 for diploid
    # temp_scale is the pixel size in the XY plane
    # determining whether that cycle shows connection with the image boundary
    bdys = set([0, 511])
    with open(temp_base_path + temp_expt_path + '/cell_cycles_compiled.pkl', 'rb') as input:
        temp_cycles = pickle.load(input)
    temp_cycles1 = []
    temp_areas = np.array([obj.segment_coords[-1].shape[0] for obj in temp_cycles])
    for obj in temp_cycles:
        obj.edge_cycle = False
        obj.pancake = False
        for i1 in range(len(obj.segment_coords)):
            temp_coords = obj.segment_coords[i1]
            if not (
                        obj.daughter is None) and obj.complete:  # this cell cycle should have been captured for the full time
                if not (obj.bud_seg[i1] is None):
                    # print 'hi'
                    # print temp_coords.shape, obj.bud_seg[i1]
                    temp_coords = np.concatenate((temp_coords, obj.bud_seg[i1]), axis=0)  # if this cell has a bud we
                    # append the segment coordinates of the bud to this area so that we flag that too.
            for i0 in range(2):
                if bdys.intersection(zip(*temp_coords)[i0]):
                    obj.edge_cycle = True  # recording that this cell cycle was segmented with a connection to the image
                    # boundary
    # determining whether that cell is too big for it to fit without being squashed

    if not (temp_size_thresh is None):
        temp_thresh = np.mean(temp_areas) + 1.5 * np.std(temp_areas)  # 1.5 sd above mean area at division.
        # temp_thresh = math.pi*(50.0*temp_height/(2*temp_scale))**2  # the area of a circle that is too big to fit
        # properly in 3D (in units of pixels squared). Note this is difficult to calibrate.
        for obj in temp_cycles:

            for temp_coords in obj.segment_coords:
                if temp_coords.shape[0] > temp_thresh:  # the length of temp_coords corresponds to the area in units of
                    # pixels squared
                    obj.pancake = True  # in this case we now know not to trust this cell cycle.
    for obj in temp_cycles:
        if not (obj.edge_cycle) and not (obj.pancake) and obj.complete and not (obj.error) and not (
            obj.daughter is None):
            # if this is a high quality cell cycle
            temp_cycles1.append(obj)
    save_object(temp_cycles1, temp_base_path + temp_expt_path + '/cell_cycles_filtered.pkl')
    print 'Number of high quality cell cycles: {0}'.format(len(temp_cycles1))
    return temp_cycles1


########################################################################################################################
# SINGLE TIMEPOINT EXPERIMENTAL ANALYSIS FROM HERE ON
########################################################################################################################


def timepoint_populate_cells_all_scenes_1(temp_base_path, temp_expt_paths, temp_expt_conds, temp_image_filenames,
                                          temp_bf_filename, temp_completed_scenes, temp_num_scenes, temp_df):
    # first we need to generate brightfield cells for each brightfield image

    # to be run from timepoint_analysis.py
    # temp_base_path is base path
    # temp_expt_paths is expt_paths
    # temp_expt_conds is expt_conds
    # temp_image_filenames is image_filenames
    # temp_bf_filename is bf_filename
    # temp_completed_scenes is completed_scenes
    color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
    dims = 512

    for i0 in range(len(temp_image_filenames)):
        temp1 =temp_df.loc[temp_df.script == 'timepoint_populate_cells_all_scenes_1', 'run_{0}'.format(i0)][0]
        if temp_completed_scenes[i0] and temp1==0:
            c = []  # this will track all the cells over this condition
            # making necessary directory tree
            directory = temp_base_path + temp_expt_paths[i0] + '/outputs'
            if not os.path.exists(directory):
                os.makedirs(directory)
            if not os.path.exists(directory + '/images'):
                os.makedirs(directory + '/images')
            if not os.path.exists(directory + '/images/bf_segments'):
                os.makedirs(directory + '/images/bf_segments')
            if not os.path.exists(directory + '/images/spot_detection'):
                os.makedirs(directory + '/images/spot_detection')
            outlines = np.zeros([temp_num_scenes[i0], dims, dims])
            for frame_num in range(1, temp_num_scenes[i0]+1):
                filename = temp_image_filenames[i0] + temp_bf_filename + '_s{0}_segmentation.mat'.format(
                    str(frame_num).zfill(3))
                path = temp_base_path + temp_expt_paths[i0] + '/bf_ims/segments' + filename
                c, im, fig, temp_outlines = timepoint_single_frame(path, c, frame_num, color_seq_val)
                outlines[frame_num - 1, :, :] = temp_outlines[:, :]
                fig.subplots_adjust(bottom=0)
                fig.subplots_adjust(top=1)
                fig.subplots_adjust(right=1)
                fig.subplots_adjust(left=0)
                # extent = mpl.transforms.Bbox(((0, 0), (5, 5)))
                fig.savefig(directory + '/images/bf_segments/scene_{0}.tif'.format(frame_num))
                del im, fig
            np.save(directory + '/cell_outlines', outlines)
            save_object(c, directory + '/cells.pkl')
            del c
            temp_df.at[temp_df.script=='timepoint_populate_cells_all_scenes_1','run_{0}'.format(i0)]=1
            print temp_df
    return temp_df


def timepoint_single_frame(temp_path, temp_cells, temp_frame, color_seq):
    # This function takes an existing set of cells, temp_cells, and will append cells from a given frame to them.
    # temp_path is the path of the .mat file with segments that is output automatically for each field of view from
    # cell star. temp_frame is an index for the frame number in question (starts from 1). tracking_csv is the path of
    # the tracking .csv file output by cellstar automatically. color_seq is the range of values in the color map.

    temp_im = sio.loadmat(temp_path)['segments']
    temp1 = np.amax(temp_im)
    # print temp1  # number of segmented cells in this image
    fig = plt.figure(figsize=[5.12, 5.12], frameon=False)  # note figsize is selected to scale with that of the image
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    temp_im1 = np.zeros(temp_im.shape)
    # plt.imshow(temp_im)
    temp_coords = []
    temp_outlines = np.zeros(temp_im.shape)
    for i0 in range(1, temp1 + 1):  # this just indexes the cells in the image, not the unique cell identifier
        temp2 = temp_im == i0  # go through sequentially
        # calculating ellipse behavior
        temp_val = fit_ellipse(temp2.astype(int))
        temp_coords.append(temp_val[4])
        temp4 = {}
        temp4['index'] = len(temp_cells)  # unique cell identifier in the order it was appended to the list of cells
        temp4['index_image'] = i0  # cell identifier in that image
        temp4['file_centroid'] = None  # placeholder
        temp4['segment_coords'] = temp_val[5]
        temp4['centroid'] = np.array(temp_val[3])
        temp4['current_frame'] = temp_frame  # adding the current frame index
        temp4['ellipse_params'] = temp_val[:4]
        temp4['ellipse_volume'] = math.pi * temp_val[0] * temp_val[1] ** 2  # note for some reason this neglected the
        # 4/3 value, so this is added later on
        temp4['2d_area'] = np.sum(temp2)  # gives the total number of pixels in the 2d plane
        temp_cells.append(Cell(temp4))  # adding a new cell
        temp_im1[temp4['segment_coords'][:, 0], temp4['segment_coords'][:, 1]] = temp4['index']
        # making the plot of cell values
    ax.imshow(temp_im1)
    for coords in temp_coords:
        temp_x, temp_y = zip(*coords)
        plt.plot(list(temp_x), list(temp_y), color=(1.0, 1.0, 1.0, 1.0))
        temp_inds1 = [i for i, e in enumerate(temp_y) if 0 > int(round(e)) or int(round(e)) > 511]
        temp_inds2 = [i for i, e in enumerate(temp_x) if 0 > int(round(e)) or int(round(e)) > 511]
        excl_inds = list(set(temp_inds1 + temp_inds2))
        # print excl_inds
        corr_inds = [i0 for i0 in range(len(temp_x)) if not (i0 in excl_inds)]
        # print corr_inds
        temp_x1 = [int(round(temp_x[i0])) for i0 in corr_inds]
        temp_y1 = [int(round(temp_y[i0])) for i0 in corr_inds]
        # print temp_outlines.shape
        temp_outlines[temp_y1, temp_x1] = 1
    plt.xlim(-0.5, 511.5)
    plt.ylim(511.5, -0.5)
    return temp_cells, ax, fig, temp_outlines


def timepoint_populate_cells_all_scenes_2(temp_base_path, temp_expt_paths, temp_expt_conds, temp_image_filenames,
                                          temp_bf_filename, temp_completed_scenes, temp_num_scenes, temp_bkgd_scenes,
                                          temp_fl_filename, temp_fl_filename_c2, temp_thresh, temp_df):
    # this takes the cell output from script populate_cells_1.py and adds fluorescence data to it for both channels
    # temp_fl_filename_c2 is the filename of the secondary fluorescence channel (primary corresponds to Whi5).
    pixel_size = {'60X': 0.267, '100X': 0.16}
    z_scale, z_offset = 0.4 / pixel_size['60X'], 0

    color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
    generate_masks = True
    bkgd_ims, bkgd_ims_c2 = [], [[], []]  #bkgd_ims_c2 has mean and then std_deviation. Note this stores it for all conds
    # loading backgrounds
    for cond in range(len(temp_image_filenames)):
        temp1 = temp_df.loc[temp_df.script == 'timepoint_populate_cells_all_scenes_2', 'run_{0}'.format(cond)][1]
        if temp_completed_scenes[cond] and temp1==0:  # if we have analyzed this scene
            filename_fl_bkgd = temp_image_filenames[cond] + temp_fl_filename_c2 + '_s{0}.TIF'.format(str(temp_bkgd_scenes[cond]))
            # channel 2
            # format is z, y, x
            temp_bkgd_im = io.imread(temp_base_path + temp_expt_paths[cond] + filename_fl_bkgd)
            bkgd_ims_c2[0].append(np.mean(temp_bkgd_im, axis=0))
            bkgd_ims_c2[1].append(np.std(temp_bkgd_im, axis=0))
            del temp_bkgd_im
            # channel 1
            filename_fl_bkgd = temp_image_filenames[cond] + temp_fl_filename + '_s{0}.TIF'.format(str(temp_bkgd_scenes[cond]))
            # format is z, y, x
            temp_bkgd_im = io.imread(temp_base_path + temp_expt_paths[cond] + filename_fl_bkgd)
            bkgd_ims.append(np.mean(temp_bkgd_im, axis=0))
            # redefining the background image to be the mean along the z axis
            del temp_bkgd_im

            # loading appropriate data for each condition
            directory = temp_base_path + temp_expt_paths[cond] + '/outputs'
            with open(directory + '/cells.pkl', 'rb') as input:
                c = pickle.load(input)
            frame_list = [obj.frames[0] for obj in c]  # gives the appropriate scene for each cell in c
            for obj in c:
                obj.add_fluor_placeholders()
            outlines = np.load(directory + '/cell_outlines.npy')
            for frame_num in range(1, temp_num_scenes[cond]+1):
                print 'scene = {0}'.format(frame_num)
                update_list=[i for i, e in enumerate(frame_list) if e==frame_num]
                # gives the list of indices that have to be addressed at each frame
                temp_fl_image_c2 = None  # placeholder for the bkgd subtracted C2 fluor image
                # Adding the fluorescence data for the second channel. We do this because the second channel is the one with
                # the strong constitutive fluor.
                filename_fl = temp_image_filenames[cond] + temp_fl_filename_c2 + '_s{0}.TIF'.format(str(frame_num))
                # taking the mean with respect to the z axis. Do it this way since there doesn't seem
                # to be any systematic bias in that direction.
                c, figs, temp_masks, temp_fl_image_c2 = timepoint_add_fluorescence_traces_c2_v1(
                    temp_path=temp_base_path + temp_expt_paths[cond] + filename_fl,
                    temp_cells=c, frame_list=update_list, current_frame=frame_num,
                    mean_bkgd=bkgd_ims_c2[0][-1],
                    std_bkgd=bkgd_ims_c2[1][-1])
                # saving the figures in figs. This is only necessary if you want to track how good the segmentation is
                # since it saves a separate image for each cell.
                del figs, temp_masks
                filename_fl = temp_image_filenames[cond] + temp_fl_filename + '_s{0}.TIF'.format(str(frame_num))

                c, mask, fig = timepoint_add_fluorescence_traces_v2(temp_path=temp_base_path + temp_expt_paths[cond] + filename_fl, temp_cells=c,
                                                     frame_list=update_list,
                                                     current_frame=frame_num, z_scaling=z_scale, z_offset=z_offset,
                                                     bkgd=bkgd_ims[-1], save_coords=False,
                                                     exists_c2=temp_fl_filename_c2, temp_outlines=outlines,
                                                          temp_im_c2=temp_fl_image_c2, temp_thresh=temp_thresh,
                                                                    temp_cond=cond)
                del temp_fl_image_c2  # Now the c2 data has been fully integrated into cells we can remove it from memory
                fig.savefig(directory+'/images/spot_detection/spot_detection_s{0}.TIF'.format(frame_num))
                del mask
                print 'done condition '+temp_expt_conds[cond][1:]+', scene {0}'.format(frame_num)
                # if frame_num>1:
                #     print c[update_list[0]].fluor_chars, update_list[0]
            save_object(c, directory + '/cells_fl.pkl')
            del c, frame_list
            temp_df.at[temp_df.script == 'timepoint_populate_cells_all_scenes_2', 'run_{0}'.format(cond)] = 1
            print temp_df


def timepoint_populate_cells_all_scenes_2_V1(temp_base_path, temp_expt_paths, temp_expt_conds, temp_image_filenames,
                                          temp_bf_filename, temp_completed_scenes, temp_num_scenes, temp_bkgd_scenes,
                                          temp_fl_filename, temp_fl_filename_c2, temp_thresh, temp_df,temp_unfinished_scenes,
                                             temp_expt_path1):
    # this takes the cell output from script populate_cells_1.py and adds fluorescence data to it for both channels
    # this differs from timepoint_populate_cells_all_scenes_2 in that it does not assume a second fluorescence channel c2
    pixel_size = {'60X': 0.267, '100X': 0.16}
    z_scale, z_offset = 0.4 / pixel_size['60X'], 0

    color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
    generate_masks = True
    bkgd_ims = []
    if not (temp_fl_filename_c2 is None):
        bkgd_ims_c2 = [[], []] #bkgd_ims_c2 has mean and then std_deviation. Note this stores it for all conds
    # loading backgrounds.
    for cond in temp_unfinished_scenes:
        temp1 = temp_df.loc[temp_df.script == 'timepoint_populate_cells_all_scenes_2', 'run_{0}'.format(cond)][1]
        if temp_completed_scenes[cond] and temp1==0:  # if we have analyzed this scene
            if not(temp_fl_filename_c2 is None):
                filename_fl_bkgd = temp_image_filenames[cond] + temp_fl_filename_c2 + '_s{0}.TIF'.format(str(temp_bkgd_scenes[cond]))
                # channel 2
                # format is z, y, x
                temp_bkgd_im = io.imread(temp_base_path + temp_expt_paths[cond] + filename_fl_bkgd)
                bkgd_ims_c2[0].append(np.mean(temp_bkgd_im, axis=0))
                bkgd_ims_c2[1].append(np.std(temp_bkgd_im, axis=0))
                del temp_bkgd_im
            # channel 1
            filename_fl_bkgd = temp_image_filenames[cond] + temp_fl_filename + '_s{0}.TIF'.format(str(temp_bkgd_scenes[cond]))
            # format is z, y, x
            temp_bkgd_im = io.imread(temp_base_path + temp_expt_paths[cond] + filename_fl_bkgd)
            bkgd_ims.append(np.mean(temp_bkgd_im, axis=0))
            # redefining the background image to be the mean along the z axis
            del temp_bkgd_im

            # loading appropriate data for each condition
            directory = temp_base_path + temp_expt_paths[cond] + '/outputs'
            with open(directory + '/cells.pkl', 'rb') as input:
                c = pickle.load(input)
            frame_list = [obj.frames[0] for obj in c]  # gives the appropriate scene for each cell in c
            for obj in c:
                obj.add_fluor_placeholders()
            outlines = np.load(directory + '/cell_outlines.npy')
            for frame_num in range(1, temp_num_scenes[cond]+1):
                print 'scene = {0}'.format(frame_num)
                update_list=[i for i, e in enumerate(frame_list) if e==frame_num]
                # gives the list of indices that have to be addressed at each frame
                temp_fl_image_c2 = None  # placeholder for the bkgd subtracted C2 fluor image
                # Adding the fluorescence data for the second channel. We do this because the second channel is the one with
                # the strong constitutive fluor.
                if not (temp_fl_filename_c2 is None):
                    filename_fl = temp_image_filenames[cond] + temp_fl_filename_c2 + '_s{0}.TIF'.format(str(frame_num))
                    # taking the mean with respect to the z axis. Do it this way since there doesn't seem
                    # to be any systematic bias in that direction.
                    c, figs, temp_masks, temp_fl_image_c2 = timepoint_add_fluorescence_traces_c2_v1(
                        temp_path=temp_base_path + temp_expt_paths[cond] + filename_fl,
                        temp_cells=c, frame_list=update_list, current_frame=frame_num,
                        mean_bkgd=bkgd_ims_c2[0][-1],
                        std_bkgd=bkgd_ims_c2[1][-1])
                    # saving the figures in figs. This is only necessary if you want to track how good the segmentation is
                    # since it saves a separate image for each cell.
                    del figs, temp_masks
                    # Note that if there is no second fluorescence channel then pixel_thresh_coords will be empty, as will
                    # all c2 values.


                filename_fl = temp_image_filenames[cond] + temp_fl_filename + '_s{0}.TIF'.format(str(frame_num))

                c, mask, fig = timepoint_add_fluorescence_traces_v2(temp_path=temp_base_path + temp_expt_paths[cond] + filename_fl, temp_cells=c,
                                                     frame_list=update_list,
                                                     current_frame=frame_num, z_scaling=z_scale, z_offset=z_offset,
                                                     bkgd=bkgd_ims[-1], save_coords=False,
                                                     exists_c2=temp_fl_filename_c2, temp_outlines=outlines,
                                                          temp_im_c2=temp_fl_image_c2, temp_thresh=temp_thresh,
                                                                    temp_cond=cond)
                del temp_fl_image_c2  # Now the c2 data has been fully integrated into cells we can remove it from memory
                fig.savefig(directory+'/images/spot_detection/spot_detection_s{0}.TIF'.format(frame_num))
                del mask
                print 'done condition '+temp_expt_conds[cond][1:]+', scene {0}'.format(frame_num)
            save_object(c, directory + '/cells_fl.pkl')
            del c, frame_list
            temp_df.at[temp_df.script == 'timepoint_populate_cells_all_scenes_2', 'run_{0}'.format(cond)] = 1
            print temp_df
            temp_df.to_csv(temp_base_path + temp_expt_path1 + '/progress_report', sep='\t')


def timepoint_add_fluorescence_traces_c2_v1(temp_path, temp_cells, frame_list, current_frame, mean_bkgd, std_bkgd):
    # v1 differs from the base by also analyzing different ways of segmenting cell volume based on fluorescence
    # This function adds the z-projected fluorescence above and below the 2d cell segmentation for a second fluorescence
    # channel c2
    temp_im = io.imread(temp_path)
    temp_mask = temp_im > np.tile(mean_bkgd + 3 * std_bkgd, [temp_im.shape[0], 1, 1])
    # generating the high pixel density mask
    temp_im1 = temp_im - mean_bkgd  # subtracting average background at each pixel value.
    figs = []
    temp_masks = []
    for temp_ind in frame_list:
        temp_binary_image = np.zeros(temp_im.shape)
        i0 = 0  # only one frame stored for each cell in a timepoint expt
        if current_frame != temp_cells[temp_ind].frames[i0]:
            print temp_ind, current_frame, temp_cells[temp_ind].frames[i0], temp_cells[temp_ind].frames
            raise ValueError('Wrong frame selected')
        temp_coords = temp_cells[temp_ind].segment_coords[i0]
        temp_cells[temp_ind].zproj_fluor_vals_c2[i0] = np.sum(temp_im1[:, temp_coords[:, 0], temp_coords[:, 1]])
        # calculating the full z projection

        # Generating the fluorescence and volume value integrated only over pixels where the fluorescence is high. Note
        # we will also do this for Whi5, since the nucleus appears to be included in the constitutive fluor channel, so
        # this should give a slightly more accurate measure for Whi5 fluorescence.
        # while the vacuoles are excluded for both
        temp_binary_image[:, temp_coords[:, 0], temp_coords[:, 1]] = 1
        temp_binary_image *= temp_mask
        temp_cells[temp_ind].pixel_thresh_coords[i0] = np.nonzero(temp_binary_image)
        temp_cells[temp_ind].pixel_thresh_fluor_vals_c2[i0] = \
            np.sum(temp_im1[temp_cells[temp_ind].pixel_thresh_coords[i0]])
        temp_masks.append(temp_binary_image)
    return temp_cells, figs, temp_masks, temp_im1


def timepoint_add_fluorescence_traces_v2(temp_path, temp_cells, frame_list, current_frame, z_scaling, z_offset, bkgd,
                               save_coords=False, exists_c2=None, temp_outlines=None, temp_im_c2=None,
                                         temp_thresh=None, temp_cond=None):
    # This is the timepoint version of add_fluorescence_traces_v2
    # We extract relevant characteristics from the fluorescence distribution such as the skewness, genie coefficient and
    # standard deviation for the cell classification.
    # z_scaling is deprecated but gives the ratio of the z stepsize in um to the x,y pixel size in um.
    # temp_outlines gives the outlines in format [y,x]
    # temp_im_c2 is still available at this stage to allow us to include data for the secondary channel in the nucleus
    temp_im1 = io.imread(temp_path)
    temp_im = temp_im1 - np.tile(bkgd, [temp_im1.shape[0], 1, 1])
    # subtracting the average background at each xy point. Structure is (z,y,x)
    temp_mask = np.zeros(temp_im.shape)  # this will be a binary mask to track the location of each cell
    # finding the blobs
    temp_im2 = skimage.filters.gaussian(temp_im1, sigma=1)
    temp = skimage.feature.blob_log(temp_im2, min_sigma=1.0, max_sigma=2.5, num_sigma=4, threshold=temp_thresh[temp_cond], overlap=0.01,
                                    log_scale=False)  # this has the format of a numpy array Nx4, where the columns are
    # z, y, x, r where r is the radius of the circle centered at that point
    # now we plot the fitted blobs along with the outlines of the cells at this timepoint
    fig = plt.figure(figsize=[10,10])
    ax = plt.subplot(1, 1, 1)
    plt.imshow(np.amax(temp_im2, axis=0))
    # plotting the cell outlines in white
    temp_coords = np.nonzero(temp_outlines[current_frame-1, :, :])
    plt.plot(temp_coords[1], temp_coords[0], '.', markersize=1, color=(1.0, 1.0, 1.0, 1.0))  # plotting in format x,y
    del temp_im2, temp_coords
    # plotting the circles of each blob
    for blob in temp:
        #     print blob
        z, y, x, r = blob
        c = plt.Circle((x, y), np.sqrt(3) * r, color='r', linewidth=2, fill=False)
        ax.add_patch(c)
    fig.subplots_adjust(bottom=0)
    fig.subplots_adjust(top=1)
    fig.subplots_adjust(right=1)
    fig.subplots_adjust(left=0)

    temp_yx = zip(*[temp[:, 1], temp[:, 2]])  # list of coordinates in type tuple
    # now we go through the set of cells in the current timepoint and determine which ones have this x,y point inside,
    # in addition to determining the fluorescence within each cell and storing that.
    for temp_ind in frame_list:
        i0 = temp_cells[temp_ind].frames.index(current_frame)  # should be zero
        # gives the index for the current frame in that cell's history
        # if current_frame != temp_cells[temp_ind].frames[i0]:
        #     print(i0, current_frame, temp_cells[temp_ind].frames[i0], temp_cells[temp_ind].frames)
        #     raise ValueError('Wrong frame selected')

        temp_coords = temp_cells[temp_ind].segment_coords[i0]  # the segmented coordinates. NOTE OVERWRITING TEMP_COORDS
        # FROM EARLIER. NO LONGER OUTLINES, NOW FULL 2D CELL MASK COORDS.
        # saving relevant fluorescence statistics
        temp_chars = list([])
        # print temp_coords
        temp_chars.append(scipy.stats.skew(temp_im[:, temp_coords[:, 0], temp_coords[:, 1]].flatten()))
        # skew
        temp_chars.append(np.std(temp_im[:, temp_coords[:, 0], temp_coords[:, 1]].flatten())/np.absolute(np.mean(temp_im[:, temp_coords[:, 0], temp_coords[:, 1]].flatten())))
        # CV
        temp_chars.append(gini(temp_im[:, temp_coords[:, 0], temp_coords[:, 1]]))  # Gini coefficient
        # print temp_chars
        temp_cells[temp_ind].fluor_chars[i0] = temp_chars
        # print temp_cells[temp_ind].fluor_chars
        # exit()
        # generating a list of fluorescence values to get good statistics. Note image formatting with z, y, x.
        temp_cells[temp_ind].zproj_fluor_vals[i0] = np.sum(temp_im[:, temp_coords[:, 0], temp_coords[:, 1]])
        # calculating the full z projection as well.
        for temp2 in temp_coords:
            temp_mask[:, temp2[1], temp2[0]] = 1  # storing these coordinates in a binary array
        if not (exists_c2 is None):  # if we have the data then we have to also add fluorescence integrated from
            # pixels segmented with respect to the second fluorescence channel
            temp_cells[temp_ind].pixel_thresh_fluor_vals[i0] = \
                np.sum(temp_im[temp_cells[temp_ind].pixel_thresh_coords[i0]])

        # Saving whether there was a blob inside the current cell and if so what its coords were
        temp_blobs = list(
            set(zip(*[temp_coords[:, 0], temp_coords[:, 1]])) & set(temp_yx))  # the list of coordinates that are in
        # both temp_coords and temp_yx. This gives the list of all blobs whose centers lie within the current cell
        # being considered.

        # code = coords of the one blob OR None if no blobs or more than one blob
        if len(temp_blobs) == 1:  # if there is exactly one blob within the current cell
            temp_deets = temp[temp_yx.index(temp_blobs[0]), :]
            # Saving the fluorescence within the nucleus. This part is computationally intensive but would require a lot of
            # work to generate something more effective

            # generating a cuboid grid around the center of the blob
            Z, Y, X = np.meshgrid(
                np.arange(np.maximum(temp_deets[0] - np.sqrt(3)*temp_deets[3], 0),
                          np.minimum(temp_deets[0] + np.sqrt(3)*temp_deets[3], temp_im.shape[0])).astype(int),
                np.arange(np.maximum(temp_deets[1] - np.sqrt(3) * temp_deets[3], 0),
                np.minimum(temp_deets[1] + np.sqrt(3) * temp_deets[3], temp_im.shape[1])).astype(int),
                np.arange(np.maximum(temp_deets[2] - np.sqrt(3) * temp_deets[3], 0),
                np.minimum(temp_deets[2] + np.sqrt(3) * temp_deets[3], temp_im.shape[2])).astype(int))

            # gives the list of coordinates within the sphere centered at the middle of the current blob
            coord_scaled = zip(Z.ravel(), Y.ravel(), X.ravel())
            temp_coords1 = [a for a in coord_scaled if np.linalg.norm(a - temp_deets[:3])<np.sqrt(3)*temp_deets[3]]
            coord_scaled1 = zip(*temp_coords1)  # this is what we will use to calculate the relevant values
            temp_cells[temp_ind].nuclear_coords[i0] = temp_coords1  # saving the coordinates of the nucleus at this
            # time point for use with the secondary fluorescence channel
            temp_cells[temp_ind].nuclear_fluor_int[i0] = np.sum(temp_im[coord_scaled1])
            temp_cells[temp_ind].nuclear_fluor_av[i0] = np.mean(temp_im[coord_scaled1])
            if not(exists_c2 is None):  # trying to use this data to infer values about the cytoplasm of the cell
                temp_cells[temp_ind].nuclear_fluor_int_c2[i0] = np.sum(temp_im_c2[coord_scaled1])  # nuclear fluor in c2
                temp_cells[temp_ind].nuclear_fluor_av_c2[i0] = np.mean(temp_im_c2[coord_scaled1])  # nuclear fluor in c2
                temp_coords2 = zip(*list(set(zip(*temp_cells[temp_ind].pixel_thresh_coords[i0])) - set(temp_coords1)))
                # cytoplasmic coordinates for this cell in format ([z], [y], [x]) that can be fed directly into np array
                # Calculating the cytoplasmic fluorescence values for this cell, both fluorescence channels
                temp_cells[temp_ind].cytoplasmic_fluor_int[i0] = np.sum(temp_im[temp_coords2])
                temp_cells[temp_ind].cytoplasmic_fluor_av[i0] = np.mean(temp_im[temp_coords2])
                temp_cells[temp_ind].cytoplasmic_fluor_int_c2[i0] = np.sum(temp_im_c2[temp_coords2])
                temp_cells[temp_ind].cytoplasmic_fluor_av_c2[i0] = np.mean(temp_im_c2[temp_coords2])
                # The majority of healthy cells have their nuclei almost completely contained within pixel_thresh_coords
                # at least within the first timepoint.
                # print temp_coords2
                # print 'Full list', temp_coords1
                # print zip(*temp_cells[temp_ind].pixel_thresh_coords[i0])
                # print temp_coords
    # print('Number of cells = {0}'.format(len(frame_list)))
    return temp_cells, temp_mask, fig


def timepoint_track_localization_manual_annotation(temp_base_path, temp_expt_paths, temp_expt_conds,
                                                   temp_image_filenames, temp_completed_scenes, temp_fl_filename,
                                                   temp_threshold, temp_drange, temp_df=None):
    # Note: Ensure that timepoint_track_localization_final.py has run successfully with the correct parameters prior to running
    # this script.
    # this function takes the annotated datasets for Whi5 localization and integrates it into the cell data, so that
    # cells are labeled as containing Whi5 in the correct timepoints. This is then fed into analyze_whi5_distribution.
    # temp_analyzed_scene was previously the scene number that was annotated to give ground truth for Whi5 localization.
    # Now deprecated.
    # temp_num_frames is the number of frames analyzed in temp_analyzed_scene.

    for cond in range(len(temp_expt_conds)):
        # print temp_df.loc[temp_df.script == 'timepoint_track_localization_manual_annotation', 'run_{0}'.format(cond)]
        temp1 = temp_df.loc[temp_df.script == 'timepoint_track_localization_manual_annotation', 'run_{0}'.format(cond)][2]
        if temp_completed_scenes[cond] and temp1==0:
            directory = temp_base_path + temp_expt_paths[cond] + '/outputs/whi5_analysis'
            temp_frame_indexes = np.load(directory + '/samples.npy')  # gives the indices of all the analyzed frames.
            # print temp_frame_indexes, len(temp_frame_indexes)
            # scene numbers
            fluor_name = temp_image_filenames[cond] + temp_fl_filename
            # drange = 65535.0
            c1 = []  # this will be a list with all the cells from the training frames
            y1 = []  # this will be the classification of cell type based on training data

            dir1 = temp_base_path + temp_expt_paths[cond] + '/outputs'
            # opening cell file
            temp_loading_path = dir1 + '/cells_fl.pkl'
            with open(temp_loading_path, 'rb') as input:
                c = pickle.load(input)
            frame_list = [obj.frames[0] for obj in c]
            outlines = np.load(dir1 + '/cell_outlines.npy')
            # print c[34].fluor_chars
            for i0 in range(len(temp_frame_indexes)):  # iterating through each frame that we have considered
                # loading cell outlines to confirm assignments
                update_list = [i for i, e in enumerate(frame_list) if e == temp_frame_indexes[i0]]
                # indices of all cells that are in the current frame
                # tracked centroids in format x, y
                temp = np.load(
                    directory + '/fl_loc_centres/scene_{0}.npy'.format(temp_frame_indexes[i0]))
                # gives the coordinate list of all annotated cells in the current frame
                coords = zip(temp[:, 0], temp[:, 1])  # format=(x,y)

                # print coords, update_list
                # print coords, temp
                # coords now stores the coordinates of cell center points for the analyzed scene, given the current frame
                # num
                # print [c[temp_ind1].pixel_thresh_fluor_vals[temp_frame_indexes[i0] - c[temp_ind1].frames[0]] for temp_ind1 in update_list]
                c, assignments = timepoint_assign_labels_3(c, update_list, coords, temp_frame_indexes[i0])
                # print [c[temp_ind1].pixel_thresh_fluor_vals[temp_frame_indexes[i0] - c[temp_ind1].frames[0]] for temp_ind1 in
                #        update_list]
                # exit()
                save_object(c, dir1 + '/cells_v1.pkl')
                c1 += [c[temp_ind1].fluor_chars[temp_frame_indexes[i0] - c[temp_ind1].frames[0]] for temp_ind1 in update_list]
                # print c1
                y1 += [c[temp_ind1].nuclear_whi5[temp_frame_indexes[i0] - c[temp_ind1].frames[0]] for temp_ind1 in update_list]
                # print y1
                temp_new = [ind for ind in range(len(c)) if
                            temp_frame_indexes[i0] in c[ind].frames and c[ind].nuclear_whi5[temp_frame_indexes[i0] - c[ind].frames[0]]]

                # organizing figure data
                temp_im = io.imread(
                    temp_base_path + temp_expt_paths[cond] + fluor_name + '_s{0}.TIF'.format(temp_frame_indexes[i0]))
                temp_im1 = temp_im / temp_drange
                # print temp_im.shape
                if np.sum(temp_im > temp_threshold) > 0:
                    temp_im1 = np.log(np.amax(temp_im1, axis=0) / np.amax(temp_im1))
                    # using a log scale in this case because this image contains the frustrating high energy pixels that
                    # sometimes crop up
                else:
                    temp_im1 = np.amax(temp_im1, axis=0) / np.amax(temp_im1)

                temp_im1 *= outlines[temp_frame_indexes[i0] - 1, :, :] == 0
                temp_im1 += outlines[temp_frame_indexes[i0] - 1, :, :].astype('uint16')
                # print c[ind].frames
                temp_cell_coords = [c[ind].position[temp_frame_indexes[i0] - c[ind].frames[0]] for ind in update_list
                                    if c[ind].nuclear_whi5[temp_frame_indexes[i0] - c[ind].frames[0]] == 1]
                # store the centroids of the correct cells.

                # plotting figures
                fig = plt.figure(figsize=[10, 10], frameon=False)
                ax = plt.Axes(fig, [0., 0., 1., 1.])
                ax.set_axis_off()
                fig.add_axes(ax)
                ax.imshow(temp_im1)
                temp_x, temp_y = zip(*temp_cell_coords)
                plt.plot(temp_x, temp_y, 'x', color='r', linestyle='None')
                fig.subplots_adjust(bottom=0)
                fig.subplots_adjust(top=1)
                fig.subplots_adjust(right=1)
                fig.subplots_adjust(left=0)
                fig.savefig(directory +
                            '/images/whi5_assignments_scene_{0}.tif'.format(temp_frame_indexes[i0]))
            np.save(directory + '/training_cells.npy', np.asarray(c1))
            np.save(directory + '/training_cells_labels.npy', np.asarray(y1))
            temp_df.at[temp_df.script == 'timepoint_track_localization_manual_annotation', 'run_{0}'.format(cond)] = 1
            print temp_df


def timepoint_assign_labels_3(temp_cells, temp_inds, temp_coords, temp_frame_num):
    # we expect that each coordinate in temp_coords should lie within a cell outline in temp_cells
    # temp_coords should be given as a list of tuples in (x,y) format
    # temp_inds should be a list of coordinates in temp_cells for cells which are present at this timestep, denoted by
    # temp_frame_num

    temp_labels = []
    temp_vals = []

    temp_timepoints = []  # temp_timepoints corresponds to the current timepoint in each cell's reference frame
    assignments_made = []  # this will detail which cells this point fits within
    for i0 in temp_inds:
        temp_timepoints.append(temp_cells[i0].frames.index(temp_frame_num))
        # gives the index at which this frame sits for the cell in question
        if temp_cells[i0].frames[temp_timepoints[-1]] != temp_frame_num:  # confirming that this is the correct index
            print temp_frame_num, temp_timepoints[-1], temp_cells[i0].frames
            raise ValueError('Wrong cells')

    # print temp_inds, temp_coords
    for temp_ind1 in range(len(temp_coords)):

        temp_comparison = []

        # assessing the scaled distance between the coordinate and the ellipse in question
        for i0 in range(len(temp_inds)):
            temp_params = temp_cells[temp_inds[i0]].ellipse_fit[temp_timepoints[i0]]
            if ellipse_contains_point(temp_params, temp_coords[temp_ind1]):
                temp_comparison.append(i0)  # this point has been classified to fall within the current ellipse.

        if len(temp_comparison) > 1:  # if this point fits within more than one ellipse we take the one with the closest
            # centroid
            ind = compare_distances(temp_cells, temp_comparison, temp_coords[temp_ind1], temp_timepoints)
            assignments_made.append(ind)
            temp_cells[temp_inds[ind]].nuclear_whi5[temp_timepoints[ind]] = 1  # this cell has Whi5 nuclear localized.
        else:
            ind = temp_comparison[0]
            temp_cells[temp_inds[ind]].nuclear_whi5[temp_timepoints[ind]] = 1  # this cell has Whi5 nuclear localized.
        assignments_made.append(ind)
    return temp_cells, assignments_made


#(temp_base_path, temp_expt_paths, temp_expt_conds,
                                                   #temp_image_filenames, temp_completed_scenes, temp_fl_filename,
                                                   #temp_threshold, temp_drange)

def timepoint_analyze_whi5_distribution(temp_base_path, temp_expt_paths, temp_expt_conds, temp_image_filenames, temp_completed_scenes,
                                        temp_fl_filename, temp_num_scenes, temp_threshold, temp_drange, temp_df):
    import seaborn as sns
    # temp_label_path determines whether we use the output from populate_cells_all_scenes_3 or 2.
    for cond in range(len(temp_expt_conds)):
        temp1 = temp_df.loc[temp_df.script == 'timepoint_analyze_whi5_distribution', 'run_{0}'.format(cond)][3]
        if temp_completed_scenes[cond] and temp1==0:
            temp_dir = temp_base_path + temp_expt_paths[cond] + '/outputs/whi5_analysis'
            temp_scenes = np.load(temp_dir + '/samples.npy')  # gives the indices of all the analyzed frames.
            # scene numbers
            if not os.path.exists(temp_dir + '/plots'):
                os.makedirs(temp_dir + '/plots')
            if not os.path.exists(temp_base_path + temp_expt_paths[cond] + '/outputs/images/auto_whi5_detection'):
                os.makedirs(temp_base_path + temp_expt_paths[cond] + '/outputs/images/auto_whi5_detection')
            X = np.load(temp_dir + '/training_cells.npy')
            y = np.load(temp_dir + '/training_cells_labels.npy')
            print cond, X
            # print X.shape, y.shape
            # exit()
            lower_bound = 0.01
            inds1 = np.nonzero(y)
            inds2 = np.nonzero(y == 0)
            for i0 in range(len(fluor_chars_names)):
                v1 = X[inds1, i0]
                v2 = X[inds2, i0]
                print v1.shape
                make_dist_plots_1 = 1
                if make_dist_plots_1:
                    fig = plt.figure(figsize=[5, 5])
                    sns.distplot(v1, label='Nuclear')
                    sns.distplot(v2, label='Cytoplasmic')
                    plt.legend()
                    plt.xlabel(fluor_chars_names[i0])
                    plt.title(fluor_chars_names[i0] + ' in manually annotated cell populations')
                    fig.savefig(temp_dir + '/plots/' + fluor_chars_names[i0])
                    del fig

            fig = plt.figure(figsize=[5, 5])
            kde = scipy.stats.gaussian_kde(v1)
            kde1 = scipy.stats.gaussian_kde(v2)
            xpoints = np.linspace(-2.0, 5.0, 100)
            vals = [kde.integrate_box_1d(low=-10, high=lim) for lim in xpoints]
            sel = [i for i, e in enumerate(vals) if e > lower_bound]
            skew_cutoff = xpoints[sel[0]]
            plt.axvline(x=xpoints[sel[0]], label='99% bound = {0}'.format(str(np.round(skew_cutoff, decimals=3))))
            print kde1.integrate_box_1d(low=skew_cutoff, high=10.0)
            plt.plot(xpoints, kde.evaluate(xpoints), label='Nuclear')
            plt.plot(xpoints, kde1.evaluate(xpoints), label='Cytoplasmic')
            plt.plot(xpoints, vals, label='Nuclear CDF')
            plt.legend()
            fig.savefig(temp_dir + '/plots/skewness_bounding')

            ##### THIS PART NEEDS TO CHANGE TO IMPLEMENT MACHINE LEARNING
            temp_clf = SVC(kernel="linear", C=0.025)  # linear svm for classification
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5, random_state=42)
            temp_clf.fit(X_train, y_train)
            score = temp_clf.score(X_test, y_test)
            print "Learning algorithm score on test data:", score
            y_test_pred = temp_clf.predict(X_test)
            out = summary_statistics(y_test_pred, y_test)
            print 'Accuracy, Sensitivity, Specificity'
            print np.around(out, 4)

            directory = temp_base_path + temp_expt_paths[cond] + '/outputs'
            temp_loading_path = directory + '/cells_fl.pkl'
            with open(temp_loading_path, 'rb') as input:
                c = pickle.load(input)
            c = automated_whi5_assignment_2(c, temp_clf)  # ML classifier
            # c = automated_whi5_assignment_1(c, skew_cutoff)  # skewness cutoff classifier
            save_object(c, directory + '/cells_v2.pkl')
            # producing plots to ensure that this tracking is adequate
            outlines = np.load(directory + '/cell_outlines.npy')
            frame_list = [obj.frames for obj in c]
            # organizing figure data
            for scene in range(1, temp_num_scenes[cond]+1):
                # checking which cells to look at for this timepoint
                temp2 = [(scene in temp1) for temp1 in frame_list]
                update_list = [i for i, e in enumerate(temp2) if e != 0]
                # loading image
                temp_im = io.imread(
                    temp_base_path + temp_expt_paths[cond] + temp_image_filenames[cond] + temp_fl_filename + '_s{0}.TIF'.format(
                        str(scene)))
                temp_im1 = temp_im / temp_drange
                if np.sum(temp_im > temp_threshold) > 0:
                    temp_im1 = np.log(np.amax(temp_im1, axis=0) / np.amax(temp_im1))
                    # using a log scale in this case because this image contains the frustrating high energy pixels that
                    # sometimes crop up
                else:
                    temp_im1 = np.amax(temp_im1, axis=0) / np.amax(temp_im1)
                temp_im1 *= outlines[scene - 1, :, :] == 0
                temp_im1 += outlines[scene - 1, :, :].astype('uint16')
                temp_cell_coords = [c[ind1].position[scene - c[ind1].frames[0]] for ind1 in update_list if
                                    (c[ind1].nuclear_whi5[
                                         scene - c[ind1].frames[0]])]  # store the centroids of the correct cells.
                # plotting figures
                fig = plt.figure(figsize=[10.0, 10.0], frameon=False)
                ax = plt.Axes(fig, [0., 0., 1., 1.])
                ax.set_axis_off()
                fig.add_axes(ax)
                ax.imshow(temp_im1, cmap='viridis')
                if len(temp_cell_coords) > 0:
                    temp_x, temp_y = zip(*temp_cell_coords)
                    plt.plot(temp_x, temp_y, '.', markersize=4, color='r', linestyle='None')
                else:
                    print 'No cells in G1 during scene {0}'.format(scene)
                fig.subplots_adjust(bottom=0)
                fig.subplots_adjust(top=1)
                fig.subplots_adjust(right=1)
                fig.subplots_adjust(left=0)
                fig.savefig(directory + '/images/auto_whi5_detection/frame_{0}.tif'.format(scene))
            for i0 in range(len(temp_scenes)):
                directory = temp_base_path + temp_expt_paths[cond] + '/outputs'
                name = '/frame_{0}.tif'.format(temp_scenes[i0])
                copyfile(directory + '/images/auto_whi5_detection' + name,
                         temp_dir + '/images' + '/whi5_assignments_scene_{0}_auto.tif'.format(temp_scenes[i0]))
                # so that we can compare them directly
            temp_df.at[temp_df.script == 'timepoint_analyze_whi5_distribution', 'run_{0}'.format(cond)] = 1
            print temp_df



def timepoint_pandas_export(temp_base_path, temp_expt_paths, temp_expt_conds, temp_expt_id, temp_completed_scenes,
                            temp_df,temp_fl_c2):
    # this function takes the whi5-assigned data and converts it into a convenient pandas array
    # loading our cells
    temp_ind = 0
    dfs=[]
    volume_element = 0.267**2*0.7  # number of microns cubed per pixel
    for cond in range(len(temp_expt_conds)):
        temp1 = temp_df.loc[temp_df.script == 'timepoint_pandas_export', 'run_{0}'.format(cond)][4]
        if temp_completed_scenes[cond] and temp1==0:
            directory = temp_base_path+temp_expt_paths[cond]+'/outputs'
            temp_loading_path = directory + '/cells_v2.pkl'  # this version of cells has Whi5 annotations made
            outline_set = set([0, 511])
            with open(temp_loading_path, 'rb') as input:
                c = pickle.load(input)
            # Adding any additional variables we need
            c1 = [obj for obj in c if not(set(zip(*obj.segment_coords[0])[0]).intersection(outline_set)) and
                  not(set(zip(*obj.segment_coords[0])[1]).intersection(outline_set))]  # the list of cells that do not
            # intersect the outline of the image
            for obj in c1:
                obj.ellipse_volume = [1.33*obj.ellipse_volume[0]*0.267**3]
                if not(temp_fl_c2 is None):
                    obj.pixel_thresh_vol = [len(obj.pixel_thresh_coords[0][0])*volume_element]
                if obj.nuclear_coords[0] is None:
                    obj.nuclear_vol = [np.nan]
                else:
                    obj.nuclear_vol = [len(obj.nuclear_coords[0])*volume_element]
            print 'Total number of cells: ', len(c), 'Total number of non-edge cells: ', len(c1)
            # print dir(c1[0])
            # generating the column vectors
            temp = dir(c1[0])[19:]  # this way we avoid all the inbuilt functions and we can just get the relevant variables
            # now we will filter for variables that are numbers
            # print temp
            cols = []
            for att in temp:
                if type(getattr(c1[0], att)) is list:
                    # print att, getattr(c1[0], att)
                    if (type(getattr(c1[0], att)[0]) in set([np.int64, np.float, np.float64, float, int])):  # if this is
                        # the right kind of variable for me (i.e. a scalar), let's get that bad boy into a pandas array
                        cols.append(att)  # we save all the relevant variables in this way
            for att in cols:
                data_vals = [[getattr(obj, att)[0] for att in cols] for obj in c1]
            cols.append('Strain')  # so that we can identify these cells
            # print cols
            df = pd.DataFrame(columns=cols)
            for i in range(len(data_vals)):
                df.loc[i] = data_vals[i]+[temp_expt_conds[cond][1:]]
            # temp_ind += i
            dfs.append(df)
            # print df, 'hello'
            # print df
    df_new = pd.concat(dfs)
    for cond in range(len(temp_expt_conds)):  # here we update both conditions simultaneously since it's not computation
        # intensive and otherwise we have to do a lot of rejigging.
        temp_df.at[temp_df.script == 'timepoint_pandas_export', 'run_{0}'.format(cond)] = 1
        print temp_df
    df_new.to_pickle('./expt_ids' + temp_expt_id + '.pkl')


def timepoint_plot_distribution(temp_base_path, temp_expt_paths, temp_expt_conds, temp_image_filenames, temp_completed_scenes,
                                        temp_fl_filename, temp_num_scenes, temp_threshold, temp_drange):
    for cond in range(len(temp_expt_conds)):
        if temp_completed_scenes[cond]:
            temp_path = temp_base_path+temp_expt_paths[cond]+'/outputs/'
            path = temp_base_path+temp_expt_paths[cond]+'/outputs/cells_v2.pkl'
            with open(path, 'rb') as input:
                c = pickle.load(input)
            y = [obj.pixel_thresh_fluor_vals[0] for obj in c if obj.nuclear_whi5[0]]  # if whi5 is in the nucleus
            z = [obj.pixel_thresh_fluor_vals_c2[0] for obj in c if obj.nuclear_whi5[0]]  # if whi5 is in the nucleus
            x = [len(obj.pixel_thresh_coords[0][0]) for obj in c if obj.nuclear_whi5[0]]  # if whi5 is in the nucleus
            fig=plt.figure(figsize=[10,5])
            plt.title(temp_expt_conds[cond])
            plt.subplot(1,2,1)
            plt.plot(x, np.asarray(y)/np.asarray(x),'.', label='Whi5-mVenNB',alpha=0.5)
            plt.xlabel('Volume')
            plt.ylabel('Integrated fluorescence')
            vals=scipy.stats.linregress(np.asarray(x),np.asarray(y)/np.asarray(x))
            print 'Whi5 corr', vals[0]
            plt.subplot(1, 2, 2)
            plt.plot(x, np.asarray(z)/np.asarray(x), '.', label='Constitutive fluor', alpha=0.5)
            vals = scipy.stats.linregress(x, np.asarray(z)/np.asarray(x))
            print 'Const corr', vals[0]
            plt.xlabel('Volume')
            plt.ylabel('Integrated fluorescence')
            plt.legend()
            fig.savefig(temp_path+'fig_cond_{0}.png'.format(cond))


def processing(df):

    df['$c2_{b,seg}$'] = df['$F2_{b,seg}$']/df['$V_{b,seg}$']
    df['$c2_{b,ell}$'] = df['$F2_{b,zproj}$']/df['$V_{b,ell}$']
    df['$c2_{s,seg}$'] = df['$F2_{s,seg}$']/df['$V_{s,seg}$']
    df['$c2_{s,ell}$'] = df['$F2_{s,zproj}$']/df['$V_{s,ell}$']

    df['$c1_{b,seg}$'] = df['$F1_{b,seg}$']/df['$V_{b,seg}$']
    df['$c1_{b,ell}$'] = df['$F1_{b,zproj}$']/df['$V_{b,ell}$']
    df['$c1_{s,seg}$'] = df['$F1_{s,seg}$']/df['$V_{s,seg}$']
    df['$c1_{s,ell}$'] = df['$F1_{s,zproj}$']/df['$V_{s,ell}$']

    df['$F2_{k,b,seg}$'] = (df['$F2_{s,seg}$']-df['$F2_{b,seg}$'])/df['$t_{G1}$']
    df['$F1_{k,b,seg}$'] = (df['$F1_{s,seg}$']-df['$F1_{b,seg}$'])/df['$t_{G1}$']
    df['$F2_{k,b,zproj}$'] = (df['$F2_{s,zproj}$']-df['$F2_{b,zproj}$'])/df['$t_{G1}$']
    df['$F1_{k,b,zproj}$'] = (df['$F1_{s,zproj}$']-df['$F1_{b,zproj}$'])/df['$t_{G1}$']

    df['$V_{div,seg}$'] = df['$V_{d,seg}$']+df['$V_{bud,seg}$']
    df['$\Delta V_{budded,seg}$'] = df['$V_{div,seg}$']-df['$V_{s,seg}$']
    df['$F1_{div,seg}$'] = df['$F1_{d,seg}$']+df['$F1_{bud,seg}$']
    df['$F2_{div,seg}$'] = df['$F2_{d,seg}$']+df['$F2_{bud,seg}$']
    df['$c1_{div,seg}$'] = df['$F1_{div,seg}$']/df['$V_{div,seg}$']
    df['$c2_{div,seg}$'] = df['$F2_{div,seg}$']/df['$V_{div,seg}$']
    df['$c1_{b,zproj}$'] = df['$F1_{b,zproj}$']/df['$V_{b,seg}$']
    df['$c2_{b,zproj}$'] = df['$F2_{b,zproj}$']/df['$V_{b,seg}$']
    df['$F1_{div,zproj}$'] = df['$F1_{d,zproj}$']+df['$F1_{bud,zproj}$']
    df['$F2_{div,zproj}$'] = df['$F2_{d,zproj}$']+df['$F2_{bud,zproj}$']
    df['$c1_{div,zproj}$'] = df['$F1_{div,zproj}$']/df['$V_{div,seg}$']
    df['$c2_{div,zproj}$'] = df['$F2_{div,zproj}$']/df['$V_{div,seg}$']
    df['$V_{div,ell}$'] = df['$V_{d,ell}$']+df['$V_{bud,ell}$']
    df['$\Delta V_{budded,ell}$'] = df['$V_{div,ell}$'] - df['$V_{s,ell}$']
    df['$\Delta V_{G1,ell}$'] = df['$V_{s,ell}$']-df['$V_{b,ell}$']
    df['$\Delta V_{full,ell}$'] = df['$V_{div,ell}$'] - df['$V_{b,ell}$']
    df['$c1_{div,ell}$'] = df['$F1_{div,zproj}$']/df['$V_{div,ell}$']
    df['$c2_{div,ell}$'] = df['$F2_{div,zproj}$']/df['$V_{div,ell}$']
    df['$c1_{div,ellseg}$'] = df['$F1_{div,seg}$']/df['$V_{div,ell}$']
    df['$c2_{div,ellseg}$'] = df['$F2_{div,seg}$']/df['$V_{div,ell}$']
    df['$t_{div}$']=df['$t_{G1}$']+df['$t_{budded}$']

    df['$\mathrm{d}c_1/\mathrm{d}t$'] = (df['$F1_{div,seg}$']-df['$F1_{b,seg}$'])/df['$t_{div}$']
    df['$\mathrm{d}c_{1,bud}/\mathrm{d}t$'] = (df['$F1_{div,seg}$']-df['$F1_{s,seg}$'])/df['$t_{budded}$']
    df['$\mathrm{d}c_2/\mathrm{d}t$'] = (df['$F2_{div,seg}$']-df['$F2_{b,seg}$'])/df['$t_{div}$']
    df['$\mathrm{d}c_{2,bud}/\mathrm{d}t$'] = (df['$F2_{div,seg}$']-df['$F2_{s,seg}$'])/df['$t_{budded}$']
    df['$\mathrm{d}V/\mathrm{d}t$'] = (df['$V_{div,ell}$'] - df['$V_{b,ell}$']) / df['$t_{div}$']
    df['$\mathrm{d}V/\mathrm{d}t/V_b$'] = (df['$V_{div,ell}$'] - df['$V_{b,ell}$']) / (df['$t_{div}$']*df['$V_{b,ell}$'])
    return df


def normalize(df):
    for temp in df.expt.unique():
        x=df.expt==temp
        temp_vol_norm = np.mean(df[x]['$V_{b,seg}$'])
        temp_F1_norm = np.mean(df[x]['$F1_{b,seg}$'])
        temp_c1_norm = np.mean(df[x]['$c1_{b,seg}$'])
        temp_F2_norm = np.mean(df[x]['$F2_{b,seg}$'])
        temp_c2_norm = np.mean(df[x]['$c2_{b,seg}$'])
        df.loc[x,'$c2_{b,seg,norm}$'] =df[x]['$c2_{b,seg}$']/temp_c2_norm
        df.loc[x,'$c2_{s,seg,norm}$'] = df[x]['$c2_{s,seg}$']/temp_c2_norm
        df.loc[x,'$c2_{div,seg,norm}$'] = df[x]['$c2_{div,seg}$']/temp_c2_norm
        df.loc[x,'$c1_{b,seg,norm}$'] = df[x]['$c1_{b,seg}$']/temp_c1_norm
        df.loc[x,'$c1_{s,seg,norm}$'] = df[x]['$c1_{s,seg}$']/temp_c1_norm
        df.loc[x,'$c1_{div,seg}$'] =df[x]['$c1_{div,seg}$']/temp_c1_norm

        df.loc[x,'$F2_{b,seg,norm}$'] = df[x]['$F2_{b,seg}$'] / temp_F2_norm
        df.loc[x,'$F1_{b,seg,norm}$'] = df[x]['$F1_{b,seg}$'] / temp_F1_norm
        df.loc[x,'$F2_{s,seg,norm}$'] = df[x]['$F2_{s,seg}$'] / temp_F2_norm
        df.loc[x,'$F1_{s,seg,norm}$'] = df[x]['$F1_{s,seg}$'] / temp_F1_norm
        df.loc[x,'$F1_{div,seg,norm}$'] = df[x]['$F1_{div,seg}$']/temp_F1_norm
        df.loc[x,'$F2_{div,seg,norm}$'] = df[x]['$F2_{div,seg}$']/temp_F2_norm

        df.loc[x,'$V_{div,seg,norm}$'] = df[x]['$V_{div,seg}$']/temp_vol_norm
        df.loc[x,'$\Delta V_{budded,norm,seg}$'] =df[x]['$\Delta V_{budded,seg}$']/temp_vol_norm
        df.loc[x,'$V_{b,seg,norm}$'] = df[x]['$V_{b,seg}$'] / temp_vol_norm
        df.loc[x,'$V_{s,seg,norm}$'] = df[x]['$V_{s,seg}$'] / temp_vol_norm
    return df


def plotting_heatmap_timepoint(xv,yv,temp_df,binrange,xlab=None,ylab=None,gridsize=15,temp_title=None):
    sns.set(font_scale=1.5)
    sns.set_style("white")
    # binrange = (15, 40)
    # xv, yv = '$V_{b,ell}$', '$t_{G1}$'
    xv1 = np.linspace(binrange[0], binrange[1], 100)
    bins = scipy.stats.binned_statistic(temp_df[xv], temp_df[yv], statistic='mean', range=binrange)
    bins1 = scipy.stats.binned_statistic(temp_df[xv], temp_df[yv], statistic='std', range=binrange)
    bins2 = scipy.stats.binned_statistic(temp_df[xv], temp_df[yv], statistic='count', range=binrange)
    fig = plt.figure()
    plt.errorbar(x=(bins.bin_edges[1:] + bins.bin_edges[:-1]) / 2, y=bins.statistic,
                 yerr=bins1.statistic / np.sqrt(bins2.statistic), fmt='.', color='r', elinewidth=3, capthick=3,
                 capsize=5, label='Binned means')
    if not (temp_title is None):
        plt.title(temp_title)
    plt.legend()
    plt.hexbin(x=temp_df[xv], y=temp_df[yv], gridsize=gridsize, cmap="BuGn")
    bottom,top = plt.ylim()
    left,right = plt.xlim()
    # sns_plot = sns.regplot(x=xv, y=yv, scatter=False, data=temp_df, ci=95, truncate=True)
    plt.ylim(bottom,top)
    plt.xlim(left,right)
    if xlab is None:
        plt.xlabel(xv)
    else:
        plt.xlabel(xlab)
    if ylab is None:
        plt.ylabel(yv)
    else:
        plt.ylabel(ylab)
    return fig


def plotting_heatmap(xv,yv,temp_df,binrange,xlab=None,ylab=None,gridsize=15,temp_title=None,temp_cmap=None,temp_xlim=None,temp_ylim=None):
    # temp_xlim should be in format (left, right)
    # temp_ylim should be in format (bottom, top)
    # xv and yv are the column labels for data to be plotted. binrange is a tuple with the lower and upper x bounds to
    # bin data within. Decreasing gridsize causes the grid to become larger.
    sns.set(font_scale=1.5)
    sns.set_style("white")
    # binrange = (15, 40)
    # xv, yv = '$V_{b,ell}$', '$t_{G1}$'
    xv1 = np.linspace(binrange[0], binrange[1], 100)
    bins = scipy.stats.binned_statistic(temp_df[xv], temp_df[yv], statistic='mean', range=binrange)
    bins1 = scipy.stats.binned_statistic(temp_df[xv], temp_df[yv], statistic='std', range=binrange)
    bins2 = scipy.stats.binned_statistic(temp_df[xv], temp_df[yv], statistic='count', range=binrange)
    fig = plt.figure()
    plt.errorbar(x=(bins.bin_edges[1:] + bins.bin_edges[:-1]) / 2, y=bins.statistic,
                 yerr=bins1.statistic / np.sqrt(bins2.statistic), fmt='.', color='k', elinewidth=3, capthick=3,
                 capsize=5, label='Binned means')
    if not (temp_title is None):
        plt.title(temp_title)
    plt.legend()
    if temp_cmap is None:
        temp_cmap1='BuGn'
    else:
        temp_cmap1=temp_cmap
    plt.hexbin(x=temp_df[xv], y=temp_df[yv], gridsize=gridsize, cmap=temp_cmap1)
    if not(temp_xlim is None):
        plt.xlim(temp_xlim[0],temp_xlim[1])
    if not(temp_ylim is None):
        plt.ylim(temp_ylim[0],temp_ylim[1])

    bottom,top = plt.ylim()
    left,right = plt.xlim()
    sns_plot = sns.regplot(x=xv, y=yv, scatter=False, data=temp_df, ci=95, truncate=True)
    plt.ylim(bottom,top)
    plt.xlim(left,right)
    if xlab is None:
        plt.xlabel(xv)
    else:
        plt.xlabel(xlab)
    if ylab is None:
        plt.ylabel(yv)
    else:
        plt.ylabel(ylab)
    return fig


# def compare_PCC()
