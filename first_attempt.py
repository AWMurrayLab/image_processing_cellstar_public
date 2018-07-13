import numpy as np
import scipy
import scipy.optimize
import scipy.io as sio
import weakref
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from scipy import stats
import pandas as pd
import cPickle as pickle
# import custom_image_toolkit as I

base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt'
cmap = mpl.cm.get_cmap('tab10')
color_seq = plt.cm.tab10(np.linspace(0.0, 1.0, 11))
# print color_seq

class Cell(object):
    cellCount = 0  # total number of cells

    def __init__(self, birth_parameters):  # birth parameters = [tb, celltype, parent, vb, parent_current]
        self.exists = True
        # these are present for all cells
        self.index = birth_parameters['index']
        self.frames = [birth_parameters['current_frame']]
        self.position = [birth_parameters['centroid']]  # Y then X
        self.index_image = [birth_parameters['index_image']]
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


def single_frame(temp_path, temp_cells, temp_frame):
    # This function takes an existing set of cells, and will append cells from a given frame to them.
    temp_im = sio.loadmat(temp_path)['segments']
    temp1 = np.amax(temp_im)
    # print temp1  # number of segmented cells in this image
    for i0 in range(1,temp1+1):  # this just indexes the cells in the image, not the unique cell identifier
        # print i0
        temp2 = temp_im == i0
        temp3 = tracking_csv[(tracking_csv['Frame_number']==temp_frame+1) & (tracking_csv[' Cell_number']==i0)]
        # note that temp_frame has 1 added here to account for the first frame being the background
        # print c
        temp4 = {}
        temp4['index'] = temp3[' Unique_cell_number'].iloc[0]  # unique cell identifier
        temp4['index_image'] = temp3[' Cell_number'].iloc[0]  # cell identifier in that image
        temp4['centroid'] = np.array([temp3[' Position_X'].iloc[0], temp3[' Position_Y'].iloc[0]])  # cell centroid
        temp4['current_frame'] = temp_frame  # adding the current frame index
        if (len(temp_cells)==0) or (not(temp3[' Unique_cell_number'].iloc[0] in [obj.index for obj in temp_cells])):
            # if this cell has not been added yet
            temp_cells.append(Cell(temp4))  # adding a new cell
        else:
            temp_ind =[obj.index for obj in temp_cells].index(temp3[' Unique_cell_number'].iloc[0])
            if ~(temp_cells[temp_ind].index==temp3[' Unique_cell_number'].iloc[0]):
                raise ValueError('Selected cell is incorrect')
            temp_cells[temp_ind].next_frame(frame_parameters=temp4)
    fig=plt.figure(figsize=[5,5])
    plt.imshow(temp_im)
    for obj in temp_cells:
        # print obj.position[0]
        plt.plot(obj.position[0][0], obj.position[0][1], '.', color=color_seq[(obj.index-1)%11])
    return temp_cells, fig


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


# initiating this over a full timelapse

scene=1
c = []  # this will track all the cells over this timelapse
for frame_num in range(1,31):
    filename = '180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w1Brightfield confocal_s{0}_t{1}_segmentation.mat'.format(str(scene),str(frame_num).zfill(2))
    path = base_path+expt_path+'/timelapse/segments/'+filename
    tracking_csv = pd.DataFrame.from_csv(base_path+expt_path+'/timelapse/segments/tracking/tracking.csv', index_col=None)
    c, im = single_frame(path, c, frame_num)
    # plt.show(im)
    im.savefig('./scene_{0}/frame_{1}.png'.format(scene,frame_num))
    del im
# for obj in c[:50]:
    # print obj.index, obj.frames
save_object(c, './scene_{0}/cells_scene_{1}.pkl'.format(scene, scene))