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

# helper

class Cell(object):
    cellCount = 0  # total number of cells

    def __init__(self, birth_parameters):  # birth parameters = [tb, celltype, parent, vb, parent_current]
        self.exists = True
        # these are present for all cells
        self.index = birth_parameters['index']
        self.frames = [birth_parameters['current_frame']]
        self.position = [birth_parameters['centroid']]  # X then Y
        self.index_image = [birth_parameters['index_image']]
        self.ellipse_fit = [birth_parameters['ellipse_params']]  # major, minor, orientation, centroid
        self.ellipse_volume = [birth_parameters['ellipse_volume']]
        self.area = [birth_parameters['2d_area']]
        self.nuclear_whi5 = [0]  # we will change this later depending on fluorescence data
        self.type = -1  # if cells are undetermined we label them a -1
        self.fluor_vals = [0]  # list of fluorescence values at each timepoint.
        self.segment_coords = [birth_parameters['segment_coords']]
        self.zproj_fluor_vals = [0]
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
        self.fluor_vals.append(0)
        self.zproj_fluor_vals.append(0)
        self.segment_coords.append(frame_parameters['segment_coords'])

    def add_fluor_placeholders(self):
        self.int_fl = []
        self.coords = []
        for temp_ind in range(len(self.frames)):
            self.int_fl.append([])
            self.coords.append([])


class CellCycle(object):
    cellCount = 0  # total number of cell cycles tracked

    def __init__(self, temp_cell, temp_parameters):  # generating a placeholder cell which will be fully updated later.
        self.index = temp_parameters['index']  # this will be the unique index of this cell cycle
        self.cell_index = temp_parameters['cell_index']  # this is the index of the cell this came from
        self.cc_num = temp_parameters['cc_num']  # this is the number of the cell cycle in the cell this came from
        self.complete = temp_parameters['complete']  # whether this corresponds to a full cell cycle or not
        self.frames = [temp_cell.frames[temp_ind] for temp_ind in temp_parameters['range']]
        self.ellipse_volume = [temp_cell.ellipse_volume[temp_ind] for temp_ind in temp_parameters['range']]
        self.int_fl = [temp_cell.int_fl[temp_ind] for temp_ind in temp_parameters['range']]
        self.ellipse_fit = [temp_cell.ellipse_fit[temp_ind] for temp_ind in temp_parameters['range']]  # major, minor, orientation, centroid
        self.data_origin = temp_parameters['data_origin']
        self.label_type = temp_cell.type
        self.nuclear_whi5 = [temp_cell.nuclear_whi5[temp_ind] for temp_ind in temp_parameters['range']]
        self.zproj_fl = [temp_cell.zproj_fluor_vals[temp_ind] for temp_ind in temp_parameters['range']]
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
        self.parent = None  # gives the index of the cell cycle for the "parent" cell in which it produced the current
        # cell (if there is one)
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
    fig = plt.figure(figsize=[5.12,5.12],frameon=False)  # note figsize is selected to scale with that of the image
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    temp_im1 = np.zeros(temp_im.shape)
    # plt.imshow(temp_im)
    temp_coords = []
    temp_outlines = np.zeros(temp_im.shape)
    for i0 in range(1, temp1+1):  # this just indexes the cells in the image, not the unique cell identifier
        # Note that this index appears to (very frustratingly) be different for the tracking_csv index and the value in
        # the image
        # print i0
        temp2 = temp_im == i0  # go through sequentially
        # calculating ellipse behavior
        temp_val = fit_ellipse(temp2.astype(int))\

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
        if np.linalg.norm(temp4['file_centroid']-temp4['centroid'])>20:
            print temp_frame, temp4['index'], temp4['index_image'], i0
            print 'stored centroid', temp4['file_centroid']
            print 'measured centroid', temp4['centroid']
            raise ValueError('Cell classification is different somehow')
        temp4['current_frame'] = temp_frame  # adding the current frame index
        temp4['ellipse_params'] = temp_val[:4]
        temp4['ellipse_volume'] = math.pi*temp_val[0]*temp_val[1]**2
        temp4['2d_area'] = np.sum(temp2)  # gives the total number of pixels in the 2d plane
        if (len(temp_cells) == 0) or (not(temp3[' Unique_cell_number'].iloc[0] in [obj.index for obj in temp_cells])):
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
        temp_inds1 = [i for i, e in enumerate(temp_y) if 0>int(round(e)) or int(round(e))>511]
        temp_inds2 = [i for i, e in enumerate(temp_x) if 0 > int(round(e)) or int(round(e)) > 511]
        excl_inds = list(set(temp_inds1+temp_inds2))
        # print excl_inds
        corr_inds = [i0 for i0 in range(len(temp_x)) if not(i0 in excl_inds)]
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


def add_fluorescence_traces_v0(temp_path, temp_cells, frame_list, current_frame, z_scaling, z_offset, bkgd):
    # z_scaling gives the ratio of the z stepsize in um to the x,y pixel size in um.
    temp_im = io.imread(temp_path)-bkgd  # subtracting the average background at each xy point
    print temp_im.shape
    for temp_ind in frame_list:
        i0 = [i for i, e in enumerate(temp_cells[temp_ind].frames) if e == current_frame][0]
        if current_frame!=temp_cells[temp_ind].frames[i0]:
            print i0, current_frame, temp_cells[temp_ind].frames[i0], temp_cells[temp_ind].frames
            raise ValueError('Wrong frame selected')
        # gives the index for the current frame in that cell's history
        temp_vals = temp_cells[temp_ind].ellipse_fit[i0]
        temp_centre = np.array([temp_vals[3][1], temp_vals[3][0], temp_im.shape[2]/2+1+z_offset])  # reversed from centroid to
        # give x, y, z here. z_offset accounts for the shift in focal plane of brightfield image relative to the
        # ideal focal plane.
        temp_scale = np.array([temp_vals[0], temp_vals[1], temp_vals[1]/z_scaling])
        Z, Y, X = np.meshgrid(range(0,temp_im.shape[0]), range(temp_im.shape[1]), range(temp_im.shape[2]))
        xlist, ylist, zlist = X.ravel(), Y.ravel(), Z.ravel()
        coord_scaled = zip(xlist, ylist, zlist)
        temp_mat = rotation_3d(-temp_vals[2])  # inverse rotation about Z axis to determine viable points
        print 'I got here'
        temp_coords = [np.linalg.norm((np.dot(temp_mat, (a - temp_centre))/temp_scale)) for a in coord_scaled]
        print 'I got further'
        # confirming that these points sit within the volume defined by the prolate ellipsoid
        temp1 = [i for i, e in enumerate(temp_coords) if e <= 1]  # indices of the appropriate points
        saved_coords = [coord_scaled[i] for i in temp1]  # saved in x, y, z format
        temp_cells[temp_ind].coords[i0] = saved_coords
        temp_int_fl = np.sum([temp_im[(temp2[2], temp2[1],temp2[0])] for temp2 in saved_coords])  # integrating the
        # fluorescence. Note inserting this in z, y, x format to be consistent with image format
        temp_cells[temp_ind].int_fl[i0] = temp_int_fl
        print 'finished cell {0}'.format(temp_ind)
    return temp_cells


def add_fluorescence_traces_v1(temp_path, temp_cells, frame_list, current_frame, z_scaling, z_offset, bkgd, save_coords=True):
    # This is the updated version of add_fluorescence_traces, where instead of comparing all points in the image
    # like a sucker, we now compare points in a bounding box around the cell roughly 1/100 the size of the full image.
    # bkgd is a np array of the same dimensions as the image in question which stores the average background at each
    # pixel
    # save_coords is a boolean value which determines whether we store the pixel coordinates of each cell.
    # z_scaling gives the ratio of the z stepsize in um to the x,y pixel size in um.
    temp_im = io.imread(temp_path)-bkgd  # subtracting the average background at each xy point. Structure is (z,y,x)
    print temp_im.shape
    temp_mask = np.zeros(temp_im.shape)  # this will be a binary mask to track the location of each cell

    for temp_ind in frame_list:
        i0 = [i for i, e in enumerate(temp_cells[temp_ind].frames) if e == current_frame][0]
        if current_frame!=temp_cells[temp_ind].frames[i0]:
            print i0, current_frame, temp_cells[temp_ind].frames[i0], temp_cells[temp_ind].frames
            raise ValueError('Wrong frame selected')
        # gives the index for the current frame in that cell's history
        temp_vals = temp_cells[temp_ind].ellipse_fit[i0]
        temp_centre = np.array([temp_vals[3][0], temp_vals[3][1], temp_im.shape[0]/2+1+z_offset])  # reversed from
        # centroid to give x, y, z here. z_offset accounts for the shift in focal plane of brightfield image relative to
        # the ideal focal plane.
        tc, temp_maj= np.floor(temp_centre), np.ceil(1.2*temp_vals[0])
        temp_scale = np.array([temp_vals[0], temp_vals[1], temp_vals[1]/z_scaling])

        # This part is different from add_fluorescence_traces_v0. We define a set of coordinates in a bounding box
        # around the cell of interest. This allows us to greatly reduce the computation time.
        X, Y, Z = np.meshgrid(np.arange(np.maximum(tc[0]-temp_maj, 0), np.minimum(tc[0]+temp_maj, temp_im.shape[2])).astype(int),
                              np.arange(np.maximum(tc[1]-temp_maj, 0), np.minimum(tc[1]+temp_maj, temp_im.shape[1])).astype(int),
                              np.arange(0, temp_im.shape[0]).astype(int))

        xlist, ylist, zlist = X.ravel(), Y.ravel(), Z.ravel()
        coord_scaled = zip(xlist, ylist, zlist)
        # print temp_centre, coord_scaled
        # raise ValueError('fucked up')
        temp_mat = rotation_3d(-temp_vals[2])  # inverse rotation about Z axis to determine viable points
        # print 'I got here'
        temp_coords = [np.linalg.norm((np.dot(temp_mat, (a - temp_centre))/temp_scale)) for a in coord_scaled]

        # print 'I got further'
        # confirming that these points sit within the volume defined by the prolate ellipsoid
        temp1 = [i for i, e in enumerate(temp_coords) if e <= 1]  # indices of the appropriate points
        saved_coords = [coord_scaled[i] for i in temp1]  # saved in x, y, z format
        xtemp, ytemp, ztemp = zip(*saved_coords)  # generating lists of points

        # saving coordinates if necessary
        if save_coords:
            temp_cells[temp_ind].coords[i0] = saved_coords
        # saving fluorescence if necessary
        temp_cells[temp_ind].fluor_vals[i0] = temp_im[ztemp, ytemp, xtemp]  # generating a list of fluorescence values
        # to get good statistics. Note image formatting with z, y, x.

        temp_int_fl = np.sum(temp_im[ztemp, ytemp, xtemp])  # integrating the
        # fluorescence. Note inserting this in z, y, x format to be consistent with image format
        temp_cells[temp_ind].int_fl[i0] = temp_int_fl
        print 'finished cell {0}'.format(temp_ind), temp_int_fl, len(saved_coords)
        for temp2 in saved_coords:
            temp_mask[(temp2[2], temp2[1], temp2[0])] = 1  # storing these coordinates in an array
    return temp_cells, temp_mask


def add_fluorescence_traces_v2(temp_path, temp_cells, frame_list, current_frame, z_scaling, z_offset, bkgd, save_coords=True):
    # This is the updated version of add_fluorescence_traces_v1, where we add both the pixels within the cell volume,
    # and also do a z-sum for the fluorescence above the 2d cell segmentation
    # save_coords is a boolean value which determines whether we store the pixel coordinates of each cell.
    # z_scaling gives the ratio of the z stepsize in um to the x,y pixel size in um.
    temp_im = io.imread(temp_path)-bkgd  # subtracting the average background at each xy point. Structure is (z,y,x)
    # temp_im1 = np.sum(temp_im, axis=0)
    # print temp_im.shape
    temp_mask = np.zeros(temp_im.shape)  # this will be a binary mask to track the location of each cell

    for temp_ind in frame_list:
        i0 = [i for i, e in enumerate(temp_cells[temp_ind].frames) if e == current_frame][0]
        if current_frame != temp_cells[temp_ind].frames[i0]:
            print i0, current_frame, temp_cells[temp_ind].frames[i0], temp_cells[temp_ind].frames
            raise ValueError('Wrong frame selected')
        # gives the index for the current frame in that cell's history
        temp_vals = temp_cells[temp_ind].ellipse_fit[i0]
        temp_centre = np.array([temp_vals[3][0], temp_vals[3][1], temp_im.shape[0]/2+1+z_offset])  # reversed from
        # centroid to give x, y, z here. z_offset accounts for the shift in focal plane of brightfield image relative to
        # the ideal focal plane.
        tc, temp_maj = np.floor(temp_centre), np.ceil(1.2*temp_vals[0])
        temp_scale = np.array([temp_vals[0], temp_vals[1], temp_vals[1]/z_scaling])

        # This part is different from add_fluorescence_traces_v0. We define a set of coordinates in a bounding box
        # around the cell of interest. This allows us to greatly reduce the computation time.
        X, Y, Z = np.meshgrid(np.arange(np.maximum(tc[0]-temp_maj, 0), np.minimum(tc[0]+temp_maj, temp_im.shape[2])).astype(int),
                              np.arange(np.maximum(tc[1]-temp_maj, 0), np.minimum(tc[1]+temp_maj, temp_im.shape[1])).astype(int),
                              np.arange(0, temp_im.shape[0]).astype(int))

        xlist, ylist, zlist = X.ravel(), Y.ravel(), Z.ravel()
        coord_scaled = zip(xlist, ylist, zlist)
        # print temp_centre, coord_scaled
        # raise ValueError('fucked up')
        temp_mat = rotation_3d(-temp_vals[2])  # inverse rotation about Z axis to determine viable points
        # print 'I got here'
        temp_coords = [np.linalg.norm((np.dot(temp_mat, (a - temp_centre))/temp_scale)) for a in coord_scaled]

        # print 'I got further'
        # confirming that these points sit within the volume defined by the prolate ellipsoid
        temp1 = [i for i, e in enumerate(temp_coords) if e <= 1]  # indices of the appropriate points
        saved_coords = [coord_scaled[i] for i in temp1]  # saved in x, y, z format
        xtemp, ytemp, ztemp = zip(*saved_coords)  # generating lists of points

        # saving coordinates if necessary
        if save_coords:
            temp_cells[temp_ind].coords[i0] = saved_coords
        # saving fluorescence if necessary
        temp_cells[temp_ind].fluor_vals[i0] = temp_im[ztemp, ytemp, xtemp]  # generating a list of fluorescence values
        # to get good statistics. Note image formatting with z, y, x.

        temp_int_fl = np.sum(temp_im[ztemp, ytemp, xtemp])  # integrating the
        # fluorescence. Note inserting this in z, y, x format to be consistent with image format
        temp_cells[temp_ind].int_fl[i0] = temp_int_fl
        temp_coords = temp_cells[temp_ind].segment_coords[i0]
        temp_cells[temp_ind].zproj_fluor_vals[i0] = np.sum(temp_im[:, temp_coords[:, 0], temp_coords[:, 1]])  # calculating
        # the full z projection as well.
        # print 'finished cell {0}'.format(temp_ind), temp_int_fl, len(saved_coords)
        for temp2 in saved_coords:
            temp_mask[(temp2[2], temp2[1], temp2[0])] = 1  # storing these coordinates in an array
    print 'Number of cells = {0}'.format(temp_ind)
    return temp_cells, temp_mask


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
    temp = np.array([[np.cos(temp_angle), -np.sin(temp_angle)],[np.sin(temp_angle), np.cos(temp_angle)]])
    return temp


def rotation_3d(temp_angle):
    # defines the SO(2) matrix for rotation of angle theta in anticlockwise direction.
    temp = np.array([[np.cos(temp_angle), -np.sin(temp_angle), 0],[np.sin(temp_angle), np.cos(temp_angle), 0], [0,0,1]])
    return temp


def fit_ellipse(temp_im):
    temp = skimage.measure.regionprops(temp_im)
    if len(temp) != 1:
        raise ValueError('Too many regions')
    # http://scikit-image.org/docs/dev/api/skimage.measure.html
    temp1 = temp[0].major_axis_length/2.0
    temp2 = temp[0].minor_axis_length/2.0
    temp3 = -temp[0].orientation  # note using negative angle here since y counts down, not up in our reckoning
    temp4 = temp[0].centroid
    temp_mat = rotation(temp3)
    # temp_coords = [reversed(np.dot(temp_mat, (temp1*np.cos(theta), temp2*np.sin(theta))) + temp4) for theta in np.linspace(0.0, 2*math.pi, 20)]
    temp_coords = [(np.dot(temp_mat, (temp1 * np.cos(theta), temp2 * np.sin(theta))) + temp4)[::-1] for theta in
                   np.linspace(0.0, 2 * math.pi, 20)]
    # print temp1, temp2, temp3
    return [temp1, temp2, temp3, temp4[::-1], temp_coords, temp[0].coords]  # major length, minor length, orientation, centroid,
    #  (centroid should be x, y), coordinates of ellipse


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
    temp_mat=[]
    assignments_made = np.zeros(len(temp_coords))
    for temp_ind1 in temp_inds:
        if temp_cells[temp_ind1].frames[0] != 1:
            print temp_cells[temp_ind1].frames[0], temp_ind1
            raise ValueError('Wrong cells')  # we want the first frame
        temp_vals.append(temp_cells[temp_ind1].ellipse_fit[0])  # getting the ellipse properties
        temp_centres.append(np.array([temp_vals[-1][3][0], temp_vals[-1][3][1]]))  # x,y format
        temp_scales.append(np.array([temp_vals[-1][0], temp_vals[-1][1]]))
        temp_mat.append(rotation(-temp_vals[-1][2]))  # inverse rotation about Z axis to determine viable points
        temp_coords1 = [np.linalg.norm((np.dot(temp_mat[-1], (a - temp_centres[-1])) / temp_scales[-1])) for a in temp_coords]
        # rotating and finding scaled distance of each labeled point from the centroid of the cell in question
        temp_labels.append([i for i, e in enumerate(temp_coords1) if e < 1])
        if len(temp_labels[-1]) > 1:
            print 'multiple points assigned to a single cell', temp_labels[-1], temp_cells[temp_ind1].index
        for i, e in enumerate(temp_coords1):
            if e<1:
                assignments_made[i] += 1  # tracking whether each point has been assigned to a cell.
        if len(temp_labels[-1]) > 0:
            temp_cells[temp_ind1].type = 1  # cells with a label in the first timeframe are assigned a 1
        else:
            temp_cells[temp_ind1].type = 0  # cells without a label in the first timeframe are assigned a 0
    temp1 = np.where(assignments_made==0)[0]
    if len(temp1)>0:
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
        temp_cells[temp_ind1].type=0  # first frame cells are all type zero unless labeled
        temp_vals.append(temp_cells[temp_ind1].ellipse_fit[0])  # getting the ellipse properties
        temp_centres.append(np.array([temp_vals[-1][3][0], temp_vals[-1][3][1]]))  # x,y format
    for temp_ind2 in range(len(temp_coords)):  # in this case we must do the assignment wrt temp_coords
        temp_coords1 = [np.linalg.norm(a - temp_coords[temp_ind2]) for a in temp_centres]
        if np.amin(temp_coords1) < 20:  # if this distance is less than the average grown cell diameter
            temp_cells[temp_inds[np.argmin(temp_coords1)]].type = 1
            temp_labels.append(np.argmin(temp_coords1))
            assignments_made[temp_ind2]=1
        else:
            print "Too far", temp_coords[temp_ind2]
    print 'Assigned points'
    for i0 in temp_labels:
        print temp_centres[i0]
    temp1 = np.where(assignments_made==0)[0]
    if len(temp1)>0:
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


def automated_whi5_assignment_1(temp_cells, temp_cutoff):
    for obj in temp_cells:
        for i0 in range(len(obj.frames)):
            obj.nuclear_whi5[i0] = scipy.stats.skew(obj.fluor_vals[i0])>temp_cutoff
    return temp_cells


def ellipse_contains_point(temp_params, temp_coord):
    # temp_params is the "ellipse_fit" property of cells, with parameters being the outputs from ellipse_fit [:4]
    # i.e. major length, minor length, orientation, centroid (centroid should be x, y). temp_coord is the position
    # being classified.
    # output is a boolean, True if temp_coord sits within the assigned ellipse, false otherwise.

    if len(temp_coord)==2:
        # 2D ellipse mapping
        temp_centre = np.array([temp_params[3][0], temp_params[3][1]])  # x,y format
        temp_scales = np.array([temp_params[0], temp_params[1]])
        temp_mat = rotation(-temp_params[2])  # inverse rotation about Z axis to determine viable points
        temp_distance = np.linalg.norm((np.dot(temp_mat, (np.asarray(temp_coord) - temp_centre)) / temp_scales))

    return temp_distance<1  # if this fits within the assigned ellipse, then we return


def compare_distances(temp_cells, temp_inds, temp_coord, temp_timepoints):
    # temp_inds is the index in temp_cells of all cells for which temp_coord fits the ellipse
    temp_centres = [temp_cells[temp_inds[i0]].position[temp_timepoints[i0]] for i0 in range(len(temp_inds))]
    temp_dists = [np.linalg.norm(np.asarray(val-temp_coord)) for val in temp_centres]
    return temp_dists.index(min(temp_dists))  # find the index of the minimum distance


def assign_md_pair(temp_cells, mother_ind, daughter_ind, temp_frame_num):
    # first just check that this division event hasn't already been tracked
    if not(temp_frame_num in temp_cells[mother_ind].daughter_assignment_frames):
        temp_cells[mother_ind].daughters.append(temp_cells[daughter_ind].index)
        temp_cells[daughter_ind].parent = temp_cells[mother_ind].index
        temp_cells[mother_ind].daughter_assignment_frames.append(temp_frame_num)
        temp_cells[daughter_ind].parent_assignment_frame = temp_cells[mother_ind].index
        # enforcing that the nuclear localization state of both cells at this timepoint is correct.
        temp_cells[mother_ind].nuclear_whi5[temp_frame_num-temp_cells[mother_ind].frames[0]] = 1
        temp_cells[daughter_ind].nuclear_whi5[temp_frame_num - temp_cells[daughter_ind].frames[0]] = 1
    return temp_cells


def create_cycles(temp_cells, temp_ind, temp_cycles, ind, temp_data_origin):
    # take a cell and generate cycles for each of its individual cycles, which will be appended to temp_cycles
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

    for temp_gen in range(len(range_points)-1):
        range_timepoints = range(range_points[temp_gen], range_points[temp_gen+1]+1)
        # gives the indices of the appropriate time points for this cell cycle (relative to the starting tracking time
        # of the cell in question)
        temp_params = {'range': range_timepoints, 'complete': None, 'parent': None, 'error': False,
                       'start': None, 'index': ind, 'cc_num': temp_gen, 'cell_index': temp_ind,
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
    temp_inds = temp_obj.temp_cycle_inds
    temp_c_inds = [obj.index for obj in temp_cells]
    temp_cc_inds = [obj.index for obj in temp_cycles]
    tv1 = [obj.index for obj in temp_cycles]
    tv2 = [temp_cycles[tv1.index(i0)] for i0 in temp_obj.temp_cycle_inds]  # gives a list with the linked cycles from a
    # given cell. Must be length at least 1.
    for i0 in range(len(temp_inds)-1):
        # assigning next generation
        temp_cycles[temp_inds[i0]].next_gen = temp_inds[i0+1]
        # assigning parent generation
        temp_cycles[temp_inds[i0+1]].prev_gen = temp_inds[i0]

    if not(temp_obj.parent is None):  # in this case, the cell in question tracked had its birth tracked
        temp_cycles[temp_inds[0]].celltype = 2  # bud
        # temp_cycles[temp_inds[0]].label_type = temp_cycles
        if len(temp_inds)>1:
            temp_cycles[temp_inds[1]].celltype = 1  # daughter
            for i1 in range(2, len(temp_inds)):
                temp_cycles[temp_inds[i1]].celltype = 0  # mother

        temp_parent = temp_cells[temp_c_inds.index(temp_obj.parent)]  # gives the parent cell
        # print len(temp_parent.temp_cycle_inds), temp_parent.daughters.index(temp_obj.index)
        # print temp_parent.temp_cycle_inds, temp_parent.daughters
        if temp_cycles[temp_inds[0]].frames[-1] in temp_parent.daughter_assignment_frames:  # only assign daughters
            # print temp_cycles[temp_inds[0]].frames[-1], temp_parent.daughter_assignment_frames
            # if this is correct
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
            # index of the cell cycle for the parent cellst
            # temp_parent_cc_index = temp_parent.temp_cycle_inds[temp_parent.daughters.index(temp_obj.index)]  # gives the
            # index of the cell cycle for the parent cell

            temp_cycles[temp_inds[0]].parent = temp_parent_cc_index
            temp_cycles[temp_parent_cc_index].bud = temp_inds[0]  # defining the bud for the parent cell
            if len(temp_inds) > 1:
                temp_cycles[temp_parent_cc_index].daughter = temp_inds[1]  # defining the daughter for the parent cell
                temp_cycles[temp_inds[1]].mother = temp_parent_cc_index  # assigning a handle for the mother cell of the
                # daughter
    elif len(temp_inds) > 1:  # in this case, we never assigned a parent for the current cell, so latter generations are
        # mothers and the first one is "unknown"
        for i1 in range(1, len(temp_inds)):
            temp_cycles[temp_inds[i1]].celltype = 0

    return temp_cells, temp_cycles


def integrate_bud_data(temp_cycles):
    indices = [obj.index for obj in temp_cycles]
    for obj in temp_cycles:
        if not(obj.daughter is None) and obj.complete:  # this cell cycle should have been captured for the full time
            temp = np.insert(np.diff(obj.nuclear_whi5), 0, 0)
            temp1 = np.where(temp==-1)[0]
            # determining the point of start in this cell cycle:
            if len(temp1) == 1:  # testing the number of "start" events in this cell cycle
                obj.start = temp1[0]
            else:
                obj.error = True

            obj.int_fl_bud = []
            obj.vbud = []
            obj.ellipse_fit_bud = []
            obj.zproj_fl_bud = []
            obj.bud_seg = []
            for temp_ind in range(len(obj.frames)):
                obj.zproj_fl_bud.append(0)
                obj.vbud.append(0)
                obj.int_fl_bud.append(0)
                obj.ellipse_fit_bud.append(0)
                obj.bud_seg.append(None)
            # print temp_cycles[obj.bud].frames, obj.frames
            bud_ind = indices.index(obj.bud)
            for ind in temp_cycles[bud_ind].frames:
                if ind in obj.frames:
                    temp_ind = obj.frames.index(ind)
                    obj.vbud[temp_ind] = temp_cycles[obj.bud].ellipse_volume[ind-temp_cycles[obj.bud].frames[0]]
                    obj.int_fl_bud[temp_ind] = temp_cycles[obj.bud].int_fl[ind-temp_cycles[obj.bud].frames[0]]
                    obj.zproj_fl_bud[temp_ind] = temp_cycles[obj.bud].zproj_fl[ind - temp_cycles[obj.bud].frames[0]]
                    obj.ellipse_fit_bud[temp_ind] = temp_cycles[obj.bud].ellipse_fit[ind - temp_cycles[obj.bud].frames[0]]
                    obj.bud_seg[temp_ind] = temp_cycles[obj.bud].segment_coords[ind-temp_cycles[obj.bud].frames[0]]
                else:
                    # print ind, obj.frames
                    obj.error = True
                    temp_cycles[obj.bud].error = True

    return temp_cycles


def generate_tree(temp_cycles, temp_ind):
    # returns a list with the unique indices of all connected cells for a cell with a given unique index
    temp_inds = [obj.index for obj in temp_cycles]
    temp_tree = [temp_inds.index(temp_ind)]
    for ind in temp_tree:
        temp_obj = temp_cycles[temp_inds.index(ind)]
        # get all the indices of directly related cells.
        temp_var = [temp_obj.daughter, temp_obj.next_gen, temp_obj.prev_gen, temp_obj.bud, temp_obj.parent]
        # remove None values
        temp_var = [ind1 for ind1 in temp_var if not(ind1 is None)]
        for ind1 in temp_var:
            if not(ind1 in temp_tree):  # only append cell indices if you know it is not already listed. IMPORTANT.
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
        i0=0
        for ind in first_frames:
            if temp_cycles[temp_inds.index(ind)].celltype == 2: # if this cell is a bud then it may have been
                # incorrectly labeled
                parent_ind = temp_cycles[temp_inds.index(ind)].parent
                if parent_ind == first_frames[i0-1]:  # if this first frame has a bud assigned then we use the parent
                    # label type
                    temp1 = temp_cycles[temp_inds.index(parent_ind)].label_type
                    if temp1 != temp_cycles[temp_inds.index(ind)].label_type:
                        print 'Bud label type was reassigned'
                        temp_type = temp1
            i0+=1
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
                if not(0 in temp_types):  # as long as these lineages are internally consistent
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
            i0+=1
            # families.append(temp_tree)
    print 'Average family size =', np.mean(size_fam)
    print 'Maximum family size =', np.amax(size_fam)
    print 'Median family size =', np.median(size_fam)
    print 'Number families =', i0
    return temp_cycles


def populate_cells_all_scenes_1(temp_base_path, temp_expt_path, temp_image_filename, temp_bf_filename,
                                temp_num_scenes, temp_num_frames):
    color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
    dims = 512
    for scene in range(1, temp_num_scenes):
        c = []  # this will track all the cells over this timelapse
        # making necessary directory tree
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene)
        if not os.path.exists(directory):
            os.makedirs(directory)
        if not os.path.exists(directory + '/images'):
            os.makedirs(directory + '/images')
        outlines = np.zeros([temp_num_frames[scene - 1], dims, dims])
        for frame_num in range(1, temp_num_frames[scene - 1]):
            filename = temp_image_filename+temp_bf_filename+'s{0}_t{1}_segmentation.mat'.format(
                str(scene), str(frame_num).zfill(2))
            path = temp_base_path + temp_expt_path + '/scene_{0}/segments/'.format(scene) + filename
            tracking_csv_file = pd.DataFrame.from_csv(
                temp_base_path + temp_expt_path + '/scene_{0}/segments/tracking/tracking.csv'.format(scene), index_col=None)
            c, im, fig, temp_outlines = single_frame(path, c, frame_num, tracking_csv_file, color_seq_val)
            outlines[frame_num - 1, :, :] = temp_outlines[:, :]

            fig.subplots_adjust(bottom=0)
            fig.subplots_adjust(top=1)
            fig.subplots_adjust(right=1)
            fig.subplots_adjust(left=0)
            # extent = mpl.transforms.Bbox(((0, 0), (5, 5)))
            fig.savefig(directory + '/images/frame_{1}.tif'.format(scene, frame_num))
            del im, fig
        np.save(directory + '/cell_outlines_scene_{0}'.format(scene), outlines)
        save_object(c, directory + '/cells_scene_{0}.pkl'.format(scene))
        del c


def populate_cells_all_scenes_2(temp_base_path, temp_expt_path, temp_image_filename, temp_fl_filename,
                                temp_num_scenes, temp_num_frames, temp_bkgd_scene):
    # this takes the cell output from script populate_cells_1.py and adds fluorescence data to it
    pixel_size = {'60X': 0.267, '100X': 0.16}
    z_scale, z_offset = 0.4 / pixel_size['60X'], 0

    color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
    generate_masks = True
    for scene in range(1, temp_num_scenes):
        print 'scene = {0}'.format(scene)
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene)
        with open(directory + '/cells_scene_{0}.pkl'.format(scene), 'rb') as input:
            c = pickle.load(input)
        frame_list = [obj.frames for obj in c]
        update_list = []
        for obj in c:
            obj.add_fluor_placeholders()
        for frame_num in range(1, temp_num_frames[scene - 1]):
            temp = [(frame_num in temp1) for temp1 in frame_list]
            update_list.append(
                [i for i, e in enumerate(temp) if e != 0])  # gives the list of indices that have to be addressed
            # at each frame
            # print update_list
            filename_fl = temp_image_filename+temp_fl_filename+ 's{0}_t{1}.TIF'.format(str(scene), str(frame_num))
            filename_fl_bkgd = temp_image_filename+temp_fl_filename + 's{0}_t{1}.TIF'.format(str(temp_bkgd_scene), str(frame_num))
            # format is z, y, x
            bkgd_im = io.imread(temp_base_path + temp_expt_path + filename_fl_bkgd)
            bkgd_im1 = np.zeros(bkgd_im.shape)
            temp1 = np.mean(bkgd_im, axis=0)
            # taking the mean with respect to the z axis. Do it this way since there doesn't seem
            # to be any systematic bias in that direction.
            for i0 in range(bkgd_im.shape[0]):
                bkgd_im1[i0, :, :] = temp1[:, :]
            del temp1, bkgd_im
            c, mask = add_fluorescence_traces_v2(temp_path=temp_base_path + temp_expt_path + filename_fl, temp_cells=c,
                                                   frame_list=update_list[-1],
                                                   current_frame=frame_num, z_scaling=z_scale, z_offset=z_offset,
                                                   bkgd=bkgd_im1, save_coords=False)
            io.imsave(directory + '/images/mask3d_s{0}_t{1}.TIF'.format(str(scene), str(frame_num)), mask)
            print 'done scene {0}, frame {1}'.format(scene, frame_num)
        save_object(c, directory + '/cells_fl_scene_{0}.pkl'.format(scene))


def populate_cells_all_scenes_2_v2(temp_base_path, temp_expt_path, temp_image_filename, temp_fl_filename,
                                temp_num_scenes, temp_num_frames, temp_bkgd_scene, temp_bkgd_details=None):
    # This is the same as populate_cells_all_scenes_2 with the difference that it allows one to not have a separate
    # background frame, but to instead define a set of coordinates in a single scene to calculate the background from.
    # To do this, temp_bkgd_details should have the format [scene_num, [xmin, xmax], [ymin, ymax]]. This is only
    # necessary if temp_bkgd_scene is None
    # this takes the cell output from script populate_cells_1.py and adds fluorescence data to it
    pixel_size = {'60X': 0.267, '100X': 0.16}
    z_scale, z_offset = 0.4 / pixel_size['60X'], 0

    color_seq_val = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
    generate_masks = True
    for scene in range(1, temp_num_scenes):
        print 'scene = {0}'.format(scene)
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene)
        with open(directory + '/cells_scene_{0}.pkl'.format(scene), 'rb') as input:
            c = pickle.load(input)
        frame_list = [obj.frames for obj in c]
        update_list = []
        for obj in c:
            obj.add_fluor_placeholders()
        for frame_num in range(1, temp_num_frames[scene - 1]):
            temp = [(frame_num in temp1) for temp1 in frame_list]
            update_list.append(
                [i for i, e in enumerate(temp) if e != 0])  # gives the list of indices that have to be addressed
            # at each frame
            # print update_list
            filename_fl = temp_image_filename+temp_fl_filename+ 's{0}_t{1}.TIF'.format(str(scene), str(frame_num))

            if not(temp_bkgd_scene is None):  # if there is a dedicated background scene
                filename_fl_bkgd = temp_image_filename+temp_fl_filename + 's{0}_t{1}.TIF'.format(str(temp_bkgd_scene), str(frame_num))
                # format is z, y, x
                bkgd_im = io.imread(temp_base_path + temp_expt_path + filename_fl_bkgd)
                bkgd_im1 = np.zeros(bkgd_im.shape)
                temp1 = np.mean(bkgd_im, axis=0)
                # taking the mean with respect to the z axis. Do it this way since there doesn't seem
                # to be any systematic bias in that direction.
                for i0 in range(bkgd_im.shape[0]):
                    bkgd_im1[i0, :, :] = temp1[:, :]
                del temp1, bkgd_im
            else:  # in this case we must provide a dedicated set of coordinates for the background
                filename_fl_bkgd = temp_image_filename + temp_fl_filename + 's{0}_t{1}.TIF'.format(str(temp_bkgd_details[0]),
                                                                                                   str(frame_num))
                bkgd_im = io.imread(temp_base_path + temp_expt_path + filename_fl_bkgd)
                print bkgd_im.shape
                temp_xmin, temp_xmax, temp_ymin, temp_ymax = temp_bkgd_details[1][0], temp_bkgd_details[1][1],\
                    temp_bkgd_details[2][0], temp_bkgd_details[2][1]
                # print temp_bkgd_details
                # print bkgd_im[:, temp_ymin:temp_ymax, temp_xmin:temp_xmax]
                temp1 = np.mean(bkgd_im[:, temp_ymin:temp_ymax, temp_xmin:temp_xmax])
                bkgd_im1 = temp1*np.ones(bkgd_im.shape)
                del temp1, bkgd_im, temp_xmin, temp_xmax, temp_ymin, temp_ymax
                # print bkgd_im1
                # exit()
            c, mask = add_fluorescence_traces_v2(temp_path=temp_base_path + temp_expt_path + filename_fl, temp_cells=c,
                                                   frame_list=update_list[-1],
                                                   current_frame=frame_num, z_scaling=z_scale, z_offset=z_offset,
                                                   bkgd=bkgd_im1, save_coords=False)
            io.imsave(directory + '/images/mask3d_s{0}_t{1}.TIF'.format(str(scene), str(frame_num)), mask)
            print 'done scene {0}, frame {1}'.format(scene, frame_num)
        save_object(c, directory + '/cells_fl_scene_{0}.pkl'.format(scene))


def populate_cells_all_scenes_3(temp_base_path, temp_expt_path, temp_label_path, temp_num_scenes):
    # This script takes the cell output from populate_cells_2 and adds data showing whether cells are labeled as being
    # stained or not. Must have generated csv files with coordinates of labeled cells for this to work.
    temp_frame_num = 1
    for scene in range(1, temp_num_scenes):
        print 'scene = {0}'.format(scene)
        with open(temp_base_path + temp_expt_path + '/scene_{0}/outputs/cells_fl_scene_{0}.pkl'.format(scene), 'rb') as input:
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
    # temp_analyzed_scene is the scene number that was annotated to give ground truth for Whi5 localization.
    # temp_num_frames is the number of frames analyzed in temp_analyzed_scene.
    directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(temp_analyzed_scene)
    fluor_name = temp_image_filename + temp_fl_filename
    # drange = 65535.0
    # opening cell file
    if temp_label_path is None:  # in this case we didn't do the labeling since we don't have labels
        temp_loading_path = directory + '/cells_fl_scene_{0}.pkl'.format(temp_analyzed_scene)
    else:
        temp_loading_path = directory + '/cells_fl_lab_scene_{0}.pkl'.format(temp_analyzed_scene)
    with open(temp_loading_path, 'rb') as input:
        c = pickle.load(input)
    frame_list = [obj.frames for obj in c]
    # loading cell outlines to confirm assignments
    outlines = np.load(directory + '/cell_outlines_scene_{0}.npy'.format(temp_analyzed_scene))
    for frame_num in range(1, temp_analyzed_frames):
        temp2 = [(frame_num in temp1) for temp1 in frame_list]
        update_list = [i for i, e in enumerate(temp2) if
                       e != 0]  # gives the indices of each object in the current frame

        temp = np.load(
            directory + '/fl_loc_centres/scene_{0}_frame_{1}.npy'.format(temp_analyzed_scene, frame_num))
        # tracked centroids in format x, y
        coords = zip(temp[:, 0], temp[:, 1])
        # print coords, temp
        # coords now stores the coordinates of cell center points for the analyzed scene, given the current frame
        # num
        c, assignments = assign_labels_3(c, update_list, coords, frame_num)
        save_object(c, directory + '/cells_scene_{0}_v1.pkl'.format(temp_analyzed_scene))

        # organizing figure data
        temp_im = io.imread(
            temp_base_path + temp_expt_path + fluor_name + 's{0}_t{1}.TIF'.format(str(temp_analyzed_scene),
                                                                                  str(frame_num)))
        temp_im1 = temp_im/temp_drange
        if np.sum(temp_im > temp_threshold) > 0:
            temp_im1 = np.log(np.amax(temp_im1, axis=0) / np.amax(temp_im1))
            # using a log scale in this case because this image contains the frustrating high energy pixels that sometimes
            # crop up
        else:
            temp_im1 = np.amax(temp_im1, axis=0) / np.amax(temp_im1)
        # print temp_im.dtype
        # exit()
        # temp_im1 = np.amax(temp_im, axis=0) / np.amax(temp_im)  # scaling for visualization purposes.
        # print np.amax(temp_im1)
        # exit()
        temp_im1 *= outlines[frame_num - 1, :, :] == 0
        temp_im1 += outlines[frame_num - 1, :, :].astype('uint16')
        # print c[ind].frames
        temp_cell_coords = [c[ind].position[frame_num - c[ind].frames[0]] for ind in update_list if c[ind].nuclear_whi5[
            frame_num - c[ind].frames[0]] == 1]  # store the centroids of the correct cells.
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
        fig.savefig(directory + '/images/whi5_assignments_frame_{0}.tif'.format(frame_num))


def analyze_whi5_distribution(temp_base_path, temp_expt_path, temp_image_filename, temp_fl_filename, temp_num_frames,
                              temp_analyzed_scene, temp_num_scenes, temp_threshold, temp_drange, temp_analyzed_frames,
                              temp_label_path):
    # temp_label_path determines whether we use the output from populate_cells_all_scenes_3 or 2.
    directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(temp_analyzed_scene)

    with open(directory + '/cells_scene_{0}_v1.pkl'.format(temp_analyzed_scene), 'rb') as input:
        c = pickle.load(input)
    if not os.path.exists(directory + '/plots'):
        os.makedirs(directory + '/plots')
    temp1 = [[], [], []]  # G1 cells
    temp2 = [[], [], []]  # G2 cells
    for obj in c:
        c1 = [obj.fluor_vals[i0] for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0] and
              obj.frames[i0]<temp_analyzed_frames]
        c2 = [obj.index for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0] and
              obj.frames[i0]<temp_analyzed_frames]
        c3 = [obj.frames[i0] for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0] and
              obj.frames[i0]<temp_analyzed_frames]
        c4 = [obj.fluor_vals[i0] for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0] == 0 and
              obj.frames[i0]<temp_analyzed_frames]
        c5 = [obj.index for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0] == 0 and
              obj.frames[i0]<temp_analyzed_frames]
        c6 = [i0 for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0] == 0 and
              obj.frames[i0]<temp_analyzed_frames]
        temp1[0] += c1
        temp1[1] += c2
        temp1[2] += c3
        temp2[0] += c4
        temp2[1] += c5
        temp2[2] += c6

    lower_bound = 0.01


    v1 = [np.amax(vals) for vals in temp1[0]]
    v2 = [np.amax(vals) for vals in temp2[0]]
    make_dist_plots_1 = 1
    if make_dist_plots_1:
        import seaborn as sns
        fig = plt.figure(figsize=[5, 5])
        sns.distplot(v1, label='Nuclear')
        sns.distplot(v2, label='Cytoplasmic')
        plt.legend()
        plt.xlabel('Max brightness')
        plt.title('Maximum brightness in manually annotated cell populations')
        fig.savefig(directory + '/plots/brightness')
        del fig


    v1 = [scipy.stats.skew(vals) for vals in temp1[0]]
    v2 = [scipy.stats.skew(vals) for vals in temp2[0]]
    make_dist_plots_1 = 1
    if make_dist_plots_1:
        import seaborn as sns
        fig = plt.figure(figsize=[5, 5])
        sns.distplot(v1, label='Nuclear')
        sns.distplot(v2, label='Cytoplasmic')
        plt.legend()
        plt.xlabel('Skewness')
        plt.title('Skewness in manually annotated cell populations')
        fig.savefig(directory + '/plots/skewness')
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
    fig.savefig(directory + '/plots/skewness_bounding')

    # exit()
    # making initial plots of single cell distributions
    make_dist_plots = 0
    if make_dist_plots:
        import seaborn as sns
        inds = [np.random.randint(0, len(temp2[0]), size=10), np.random.randint(0, len(temp1[0]), size=10)]
        labels = ['Whi5 non-nuclear localized', 'Whi5 nuclear localized']
        for i0 in range(2):
            for ind1 in inds[i0]:
                # print ind
                fig = plt.figure(figsize=[5, 5])
                plt.xlabel('Whi5 fluorescence intensity')
                plt.title(labels[i0])
                if i0 == 1:
                    sns.distplot(temp1[0][ind1], label='cell # {0}, frame # {1}'.format(temp1[1][ind1], temp1[2][ind1]))
                    plt.legend()
                    fig.savefig(
                        directory + '/plots/G1_{0}_cell_{1}_frame_{2}'.format(i0, temp1[1][ind1], temp1[2][ind1]))
                else:
                    sns.distplot(temp2[0][ind1], label='cell # {0}, frame # {1}'.format(temp2[1][ind1], temp2[2][ind1]))
                    plt.legend()
                    fig.savefig(
                        directory + '/plots/G1_{0}_cell_{1}_frame_{2}'.format(i0, temp2[1][ind1], temp2[2][ind1]))

        for i0 in range(2):
            for ind1 in inds[i0]:
                # print ind
                fig = plt.figure(figsize=[5, 5])
                plt.xlabel('Whi5 fluorescence intensity')
                plt.title(labels[i0])
                if i0 == 1:
                    sns.distplot(temp1[0][ind1] / np.mean(temp1[0][ind1]),
                                 label='cell # {0}, frame # {1}'.format(temp1[1][ind1], temp1[2][ind1]))
                    plt.legend()
                    fig.savefig(
                        directory + '/plots/norm_G1_{0}_cell_{1}_frame_{2}'.format(i0, temp1[1][ind1], temp1[2][ind1]))
                else:
                    sns.distplot(temp2[0][ind1] / np.mean(temp2[0][ind1]),
                                 label='cell # {0}, frame # {1}'.format(temp2[1][ind1], temp2[2][ind1]))
                    plt.legend()
                    fig.savefig(
                        directory + '/plots/norm_G1_{0}_cell_{1}_frame_{2}'.format(i0, temp2[1][ind1], temp2[2][ind1]))

    # Automatically classifying cells based on skewness cutoff
    automate_analysis = 1
    if automate_analysis:
        for scene in range(1, temp_num_scenes):
            print 'Working on Scene number {0}'.format(scene)
            directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene)
            if temp_label_path is None:  # in this case we didn't do the labeling since we don't have labels
                temp_loading_path = directory + '/cells_fl_scene_{0}.pkl'.format(scene)
            else:
                temp_loading_path = directory + '/cells_fl_lab_scene_{0}.pkl'.format(scene)
            with open(temp_loading_path, 'rb') as input:
                c = pickle.load(input)
            #
            # with open(directory + '/cells_fl_lab_scene_{0}.pkl'.format(scene), 'rb') as input:
            #     c = pickle.load(input)
            c = automated_whi5_assignment_1(c, skew_cutoff)
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
                    temp_base_path + temp_expt_path +temp_image_filename+ temp_fl_filename + 's{0}_t{1}.TIF'.format(str(scene),
                                                                                                str(frame_num)))
                temp_im1 = temp_im/ temp_drange
                if np.sum(temp_im > temp_threshold) > 0:
                    temp_im1 = np.log(np.amax(temp_im1, axis=0) / np.amax(temp_im1))
                    # using a log scale in this case because this image contains the frustrating high energy pixels that
                    # sometimes crop up
                else:
                    temp_im1 = np.amax(temp_im1, axis=0) / np.amax(temp_im1)
                temp_im1 *= outlines[frame_num - 1, :, :] == 0
                temp_im1 += outlines[frame_num - 1, :, :].astype('uint16')
                # print c[ind].frames
                # print [frame_num - c[ind1].frames[0] for ind1 in update_list if
                # (c[ind1].nuclear_whi5[int(frame_num - c[ind1].frames[0])])]
                temp_cell_coords = [c[ind1].position[frame_num - c[ind1].frames[0]] for ind1 in update_list if
                                    (c[ind1].nuclear_whi5[
                                         frame_num - c[ind1].frames[0]])]  # store the centroids of the correct cells.
                # print temp_cell_coords

                # print temp_cell_coords
                # print temp_cell_coords

                # plotting figures
                fig = plt.figure(figsize=[5.12, 5.12], frameon=False)
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
                # exit()


def assign_troublesome_pair(temp_posn, temp_frame_num, temp_dir, temp_cell_num, temp_scene):
    # This will prompt the user to find the position of the partner cell at that timepoint and add it to a .txt tab
    # delimited file for this scene and frame number. The format should be
    # cell  X   Y
    # num   x1  y1
    print 'Assign pair cell for scene {0}, frame {1}, cell number {2}'.format(temp_scene, temp_frame_num,
                                                                              temp_cell_num)
    print 'position: ', temp_posn
    temp_complete = input('Cell partner assigned: ')
    if not os.path.exists(temp_dir+'/cell_coords'):
        os.makedirs(temp_dir+'/cell_coords')
    if temp_complete == 'y':
        temp_posns_pd = pd.read_csv(temp_dir+'/cell_coords_scene_{0}_frame_{1}.csv'.format(temp_scene, temp_frame_num),
                                    sep='\t')
        temp1 = temp_posns_pd.cell == temp_cell_num
        temp_posn_partner = np.array([temp_posns_pd[temp1].X[0], temp_posns_pd.Y[0]])
    else:
        temp_posn_partner = None
    return temp_posn_partner


def assign_lineages(temp_base_path, temp_expt_path, temp_image_filename, temp_fl_filename, temp_num_frames,
                              temp_analyzed_scene, temp_num_scenes, temp_threshold, temp_drange):
    # This function takes as inputs the current cells with specified Whi5 localization, and assigns lineages
    distance_cutoff = 40  # cutoff for centroid distance above which mother-daughter pairings cannot be assigned
    num_divisions = 0
    num_nuclear_events = 0
    for scene in range(1, temp_num_scenes):
        print 'Starting scene', scene
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene)
        with open(directory + '/cells_scene_{0}_v2.pkl'.format(scene), 'rb') as input:
            c = pickle.load(input)
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
            if frame_num%10 == 0:
                print 'Reached Scene {0}, frame {1}'.format(scene, frame_num)

            # Now we check which cells to analyze for each frame
            assigned_inds = [[], []]
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
                    if temp_ind == 0:  # if this cell was segmented for the first time in this frame then we store that,
                        # but we don't do much with it for now.
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
                        temp_started_G1[0].append(ind1)  # these are the cells we have to assign
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

                # We now go through each cell in temp_started_G1 and see what other cell it is closest to.
                for ind1 in range(len(temp_started_G1[0])):
                    # evaluating distances between cells
                    d = 1000 * np.ones(len(temp_started_G1[0]))
                    for ind2 in range(len(temp_started_G1[0])):
                        if ind1 != ind2:
                            # print d, temp_started_G1
                            d[ind2] = np.linalg.norm(temp_started_G1[2][ind1] - temp_started_G1[2][ind2])

                    if temp_mothers[ind1]:  # if this cell is a mother then look within the cells which haven't
                        # previously been classified as mothers.

                        temp_inds = np.nonzero(temp_mothers == 0)
                        temp_ind2 = np.argsort(
                            d[temp_inds])  # ordering the indices of temp_inds based on their distance
                        # print len(temp_inds), temp_ind2, temp_inds[0]
                        if len(temp_inds[0]) > 0:
                            if d[temp_inds[0][temp_ind2[0]]] < distance_cutoff:  # if the lowest distance is less than
                                # the distance cutoff
                                if sum(d[
                                           temp_inds] < distance_cutoff) > 1:  # if there is another cell with low distance in
                                    # the same frame then we neglect this for now. This should ideally be developed later
                                    # with manual correction
                                    # pending_assignment[-1].append(temp_started_G1[0][ind1])
                                    # dir1 = temp_base_path + temp_expt_path
                                    temp_loc = assign_troublesome_pair(temp_started_G1[2][ind1], frame_num,
                                                    directory, temp_started_G1[0][ind1], scene)
                                    # Calculating the distance from this point to each cell that started G1 this
                                    # cell cycle
                                    d1 = np.ones(len(temp_started_G1[0]))
                                    for ind2 in range(len(temp_started_G1[0])):
                                            d1[ind2] = np.linalg.norm(
                                                temp_started_G1[2][ind2] - temp_loc)
                                    if np.amin(d1)<20 and np.amin(d1) in temp_inds[0]:
                                        # if we can assign the point to a cell that has newly started G1
                                        assigned_inds[0].append(temp_started_G1[0][ind1])
                                        assigned_inds[1].append(temp_started_G1[0][np.argmin(d1)])
                                        c = assign_md_pair(c, mother_ind=temp_started_G1[0][ind1],
                                                           daughter_ind=temp_started_G1[0][np.argmin(d1)],
                                                           temp_frame_num=frame_num)
                                        num_divisions += 1
                                    else:
                                        print 'Unable to assign a daughter in scene {0}, frame{1}, cell{2}'.format(
                                            scene, frame_num, temp_started_G1[0][ind1])
                                else:
                                    assigned_inds[0].append(temp_started_G1[0][ind1])
                                    assigned_inds[1].append(temp_started_G1[0][temp_inds[0][temp_ind2[0]]])
                                    c = assign_md_pair(c, mother_ind=temp_started_G1[0][ind1],
                                                         daughter_ind=temp_started_G1[0][temp_inds[0][temp_ind2[0]]],
                                                         temp_frame_num=frame_num)
                                    num_divisions += 1
                        else:
                            pending_assignment[-1].append(temp_started_G1[0][ind1])  # if we can't assign this well,
                            # then we will revisit it in the next timeframe
                    else:
                        temp_ind2 = np.argsort(d)
                        if d[temp_ind2[0]] < distance_cutoff:
                            if sum(d < distance_cutoff) > 1:  # As above for the case of cells which we know to be
                                # mothers, if there is another cell with low distance in
                                # the same frame then we neglect this for now. This should ideally be developed later
                                # with manual correction
                                # pending_assignment[-1].append(temp_started_G1[0][ind1])
                                temp_loc = assign_troublesome_pair(temp_started_G1[2][ind1], frame_num,
                                                                   directory, temp_started_G1[0][ind1], scene)
                                # Calculating the distance from this point to each cell that started G1 this
                                # cell cycle
                                d1 = np.ones(len(temp_started_G1[0]))
                                for ind2 in range(len(temp_started_G1[0])):
                                    d1[ind2] = np.linalg.norm(
                                        temp_started_G1[2][ind2] - temp_loc)
                                if np.amin(d1) < 20:
                                    # if we can assign the point to a cell that has newly started G1
                                    assigned_inds[0].append(temp_started_G1[0][ind1])
                                    assigned_inds[1].append(temp_started_G1[0][np.argmin(d1)])
                                    c = assign_md_pair(c, mother_ind=temp_started_G1[0][ind1],
                                                       daughter_ind=temp_started_G1[0][np.argmin(d1)],
                                                       temp_frame_num=frame_num)
                                    num_divisions += 1
                                else:
                                    print 'Unable to assign a daughter in scene {0}, frame{1}, cell{2}'.format(
                                        scene, frame_num, temp_started_G1[0][ind1])
                            else:  # we add a new pair
                                temp_inds1 = [temp_started_G1[0][ind1], temp_started_G1[0][temp_ind2[0]]]
                                # print temp_inds1
                                temp_val = [c[temp_inds1[0]].ellipse_volume[frame_num - c[temp_inds1[0]].frames[0]],
                                            c[temp_inds1[1]].ellipse_volume[frame_num - c[temp_inds1[1]].frames[0]]]
                                # print temp_val
                                m_ind = np.argmax(temp_val)
                                # find which is the mother cell based on which has the larger cell size.
                                assigned_inds[0].append(temp_inds1[m_ind])
                                assigned_inds[1].append(temp_inds1[m_ind - 1])
                                c = assign_md_pair(c, mother_ind=temp_inds1[m_ind],
                                                     daughter_ind=temp_inds1[m_ind - 1],
                                                     temp_frame_num=frame_num)
                                num_divisions += 1
                        else:
                            pending_assignment[-1].append(temp_started_G1[0][ind1])  # if we can't assign this well,
                            # then we will plan to revisit it in the next timeframe

                    # we take any cell that is pending assignment in the current frame, and compare to the previous frame to
                    # see if there is an appropriate assignment there.
                    for i0 in range(len(pending_assignment[-1])):
                        coord1 = c[pending_assignment[-1][i0]].position[
                            frame_num - c[pending_assignment[-1][i0]].frames[0]]
                        d = []
                        for i1 in range(len(pending_assignment[-2])):
                            if frame_num in c[pending_assignment[-2][i1]].frames:  # if the cell in the previous timept
                                # still exists in the current one we can assign it
                                # compare the positions in the current frame and the previous frame
                                coord2 = c[pending_assignment[-2][i1]].position[frame_num - 1 -
                                                                                c[pending_assignment[-2][i1]].frames[0]]
                                d.append(np.linalg.norm(coord1 - coord2))
                            else:
                                d.append(distance_cutoff+1000)  # we will not risk assigning this lineage since the cell
                                # stopped being tracked in the current timepoint.
                        temp_ind2 = np.argsort(d)
                        if len(temp_ind2) != 0:
                            if d[temp_ind2[0]] < distance_cutoff:
                                if len(d) > 1:
                                    # print d, np.asarray(d) < distance_cutoff
                                    if sum(np.asarray(d) < distance_cutoff) == 1:
                                        temp_inds1 = [pending_assignment[-1][i0], pending_assignment[-2][temp_ind2[0]]]
                                        # print temp_inds1
                                        temp_val = [
                                            c[temp_inds1[0]].ellipse_volume[frame_num - c[temp_inds1[0]].frames[0]],
                                            c[temp_inds1[1]].ellipse_volume[frame_num - 1 - c[temp_inds1[1]].frames[0]]]
                                        m_ind = np.argmax(temp_val)
                                        # find which is the mother cell based on which has the larger cell size.
                                        assigned_inds[0].append(temp_inds1[m_ind])
                                        assigned_inds[1].append(temp_inds1[m_ind - 1])
                                        c = assign_md_pair(c, mother_ind=temp_inds1[m_ind],
                                                             daughter_ind=temp_inds1[m_ind - 1],
                                                             temp_frame_num=frame_num - 1)
                                        num_divisions += 2
                                        # we assign these cells as both being present in the previous frame, since that's when
                                        # the nuclear localization event must have occurred for both
                                        to_remove[0].append(pending_assignment[-1][i0])
                                        to_remove[1].append(pending_assignment[-2][temp_ind2[0]])
                                else:
                                    temp_inds1 = [pending_assignment[-1][i0], pending_assignment[-2][temp_ind2[0]]]
                                    # print temp_inds1
                                    temp_val = [c[temp_inds1[0]].ellipse_volume[frame_num - c[temp_inds1[0]].frames[0]],
                                                c[temp_inds1[1]].ellipse_volume[
                                                    frame_num - 1 - c[temp_inds1[1]].frames[0]]]
                                    m_ind = np.argmax(temp_val)
                                    # find which is the mother cell based on which has the larger cell size.
                                    assigned_inds[0].append(temp_inds1[m_ind])
                                    assigned_inds[1].append(temp_inds1[m_ind - 1])
                                    c = assign_md_pair(c, mother_ind=temp_inds1[m_ind],
                                                         daughter_ind=temp_inds1[m_ind - 1],
                                                         temp_frame_num=frame_num - 1)
                                    num_divisions += 2
                                    # we assign these cells as both being present in the previous frame, since that's when
                                    # the nuclear localization event must have occurred for both. They will appear in the
                                    # output images as being tracked in the current frame.
                                    to_remove[0].append(pending_assignment[-1][i0])
                                    to_remove[1].append(pending_assignment[-2][temp_ind2[0]])
                                    # removing the newly tracked items from pending_assignment. Cells from the previous timestep will still
                                    # appear as not being tracked in the previous frame.
                pending_assignment[-1] = list(set(pending_assignment[-1]) - set(to_remove[0]))
                pending_assignment[-2] = list(set(pending_assignment[-2]) - set(to_remove[1]))
            # producing a figure to track the result of this so far.
            # if frame_num ==50 and scene==4:
            #     for ind1 in assigned_inds[1]:
            #         print ind1, len(c[ind1].position), frame_num-c[ind1].frames[0], c[ind1].nuclear_whi5
            #         print c[ind1].position
            temp_im = io.imread(
                temp_base_path + temp_expt_path + temp_image_filename+temp_fl_filename+
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
    ind = 0  # this indexes the cell we are considering in the current timestep
    # this will be repeated with each new iteration

    for scene in range(1, temp_num_scenes):
        print 'Scene Number: {0}'.format(scene)
        data_index = [temp_base_path + temp_expt_path, scene]
        directory = temp_base_path + temp_expt_path + '/scene_{0}/outputs'.format(scene)
        with open(directory + '/cells_scene_{0}_v3.pkl'.format(scene), 'rb') as input:
            c = pickle.load(input)
        for i0 in range(len(c)):
            if not(i0 == c[i0].index):
                raise ValueError('Cell indexing distinct from cell list indexing')
            else:
                c, cc, ind = create_cycles(c, i0, cc, ind,
                                             temp_data_origin=data_index)  # note that the unique index for each
            # cell cycle is unique across scenes also, but that the
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
    with open(temp_base_path+temp_expt_path+ '/cell_cycles_compiled.pkl', 'rb') as input:
        temp_cycles = pickle.load(input)
    filt_cc = [obj for obj in temp_cycles if obj.complete and not (obj.error) and not (obj.daughter is None)]
    for scene in range(1, temp_num_scenes):
        scene_cycles = [obj for obj in filt_cc if obj.data_origin[1] == scene]  # selecting the cycles in this scene
        frame_list = [obj.frames for obj in scene_cycles]
        for frame_num in range(1, temp_num_frames[scene-1]):
            update_list = [i for i, e in enumerate(frame_list) if frame_num in e]  # gives the number of objects in the
            # current frame
            # load the current image
            filename = temp_image_filename + temp_bf_filename + 's{0}_t{1}.TIF'.format(
                str(scene), str(frame_num).zfill(2))
            temp_im = io.imread(temp_base_path+temp_expt_path+'/scene_{0}'.format(scene)+filename)/65535.0
            # convert this grayscale to rgb
            temp_im1 = np.repeat(temp_im[:, :, np.newaxis], 3, axis=2)
            for temp_ind in update_list:  #  going through to add segmented visualization
                temp_coords = zip(*scene_cycles[temp_ind].segment_coords[frame_num-scene_cycles[temp_ind].frames[0]])
                temp_mask_im = np.zeros(temp_im.shape)
                temp_mask_im[temp_coords] = 1  # creating the mask for the main cell
                if not(scene_cycles[temp_ind].bud_seg[frame_num-scene_cycles[temp_ind].frames[0]] is None):
                    # if we have coords for the bud at this stage we add that data
                    temp_coords = zip(*scene_cycles[temp_ind].bud_seg[frame_num-scene_cycles[temp_ind].frames[0]])
                    temp_mask_im[temp_coords] = 1  # creating the mask for the bud
                temp_mask_im1 = np.repeat(temp_mask_im[:,:,np.newaxis], 3, axis=2)  # now in rgb shape!
                temp_rgb_vals = cm.tab20(scene_cycles[temp_ind].index % 20)  # cyclic color map so that we track cell
                # cycles
                temp_mask_im2 =temp_mask_im1==0
                temp_im1*=temp_mask_im2
                temp_im1+=temp_mask_im1*np.asarray(temp_rgb_vals[:3])
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
            fig.savefig(temp_base_path+temp_expt_path+'/scene_{0}'.format(scene)+'/outputs' +
                        '/images/cell_cycle_validation_frame_{1}.tif'.format(scene, frame_num))
            del fig


def label_boundaries(temp_base_path, temp_expt_path, temp_image_filename, temp_bf_filename, temp_num_scenes,
                     temp_num_frames):
    cells = []
    bdys = set([0, 511])
    for scene in range(1, temp_num_scenes):
        with open(temp_base_path+temp_expt_path + '/scene_{0}/outputs/cells_fl_scene_{0}.pkl'.format(scene), 'rb') as input:
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
    save_object(cells, temp_base_path+temp_expt_path + '/cells_compiled.pkl')
