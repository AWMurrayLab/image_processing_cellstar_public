import numpy as np
import matplotlib.pyplot as plt
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import manual_image_toolkit as M
import os
import scipy
from scipy import stats

# THIS SCRIPT IS UNFINISHED. DO NOT USE.


# this script takes a list of cell tracks with daughter lineages assigned (e.g. the output from assign_lineages_1.py)
# and creates a collection of reliable cell cycles from this by iterating through this list.

base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
cc = []
unique_ind = 0
global ind
ind = 0  # this indexes the cell we are considering in the current timestep
# this will be repeated with each new iteration
for scene in range(1, 4):
    directory = base_path+expt_path+'/scene_{0}/outputs'.format(scene)
    with open(directory + '/cells_scene_{0}_v3.pkl'.format(scene), 'rb') as input:
        c = pickle.load(input)
    lookup_inds = [obj.index for obj in c]




    # by definition, the first cell cycle cannot be tracked fully, but we keep track of it here regardless since it's
    # important for lineages. This is done in such a way that if a daughter starts being tracked at the same time
    # as whi5 enters the nucleus, we assign it as a full cell cycle.
    range_timepoints = range(0, temp_first_timepoints[0]+1)  # this necessarily has length at least 2, since it must
    # exist at least 1 timepoint before the difference variable kicks in.
    temp_params = {'range': range_timepoints, 'cc_num': 0, 'parent': None, 'index': ind, 'unique_ind': unique_ind}

    # create a cell cycle for this event
    cc.append(C.CellCycle(c[ind], temp_params))
    unique_ind += 1  # this is different for each cell cycle.
    # if num_cycles >= 1:  # if there is at least one cell cycle here, otherwise we simply create a
    # placeholder cell which notes that it was not tracked for a full cell cycle
    for i0 in range(num_cycles):  # num_cycles = len(temp_first_timepoints)-1
        range_timepoints = range(temp_first_timepoints[i0]+1, temp_first_timepoints[i0+1]+1)  # note that because this
        # took a difference, the actual timepoints of interest are shifted backward by 1.
        cc.append(C.CellCycle(c[ind], range_timepoints, i0))


def iterate_cells(temp_cells, temp_cycles, temp_ind, temp_lookup):
    # take a cell cycle with a given index, populate it with values, and then output the next generation cells (if there
    # are any).
    # input temp_cycles[temp_ind] is assumed to have a properties 'cell_index' pointing to the cell it came from. It is
    # assumed to have a parent already assigned (either none, or a previous cell cycle). It is assumed to have a param
    # 'cc_num' which tracks the number of the cell cycle to be considered
    temp_obj = temp_cells[temp_lookup.index[temp_cycles[temp_ind].cell_index]]
    # gives the cell that temp_cycles[temp_ind] came from.
    temp = np.diff(np.array(temp_obj.nuclear_whi5))  # 1 gives entering G1, -1 gives leaving G1.
    temp_first_timepoints = np.where(temp == 1)
    temp_first_timepoints += 1  # now this gives the index in the nuclear_whi5 parameter of the appropriate event
    range_points = np.insert(temp_first_timepoints, 0, 0)
    range_points = np.append(range_points, len(temp))  # range_points can now be iterated through to generate cycles of
    # the appropriate length, with the condition that if a cell stopped being tracked just as a new cycle started, this
    # will still work.
    num_cycles = len(temp_first_timepoints) - 1
    # this gives the number of complete cell cycles tracked, ignoring for now
    # whether daughters were assigned at each timepoint.

    temp_gen = temp_cycles[temp_ind].cc_num  # this gives the number of the cell cycle that is being tracked
    range_timepoints = range(range_points[temp_gen], range_points[temp_gen+1]+1)  # gives the indices of the appropriate
    # time points for this cell cycle (relative to the starting tracking time of the cell in question)
    temp_params = {'range': range_timepoints, 'complete': None, 'daughter': None, 'next_gen': None, 'error': None,
                   'start': None}
    if temp_gen == 0 and temp_obj.parent_assignment_frame[0] == temp_obj.frames[0]:
        # if this is a daughter that came into view just as the division was triggered. Deal with this separately since
        # it can potentially be classified as a complete cell cycle and otherwise would not. For now we are not.
        # temp_params['parent'] = c[ind].parent[0]
        if __name__ == '__main__':
            temp_params['complete'] = False  # We say that this cell cycle isn't complete, so they won't be considered.

        # if num_cycles >= 1:  # if we tracked this for a complete cycle
        #     temp_params['complete'] = True
        # else:
        #     temp_params['complete'] = False  # in this case we caught the start but not the end
    elif range_timepoints[0] in temp_first_timepoints and range_timepoints[-1] in temp_first_timepoints:
        temp_params['complete'] = True  # if the start and end points fit the description of temp_first_timepoints this
        # corresponds to a full cell cycle
    else:
        temp_params['complete'] = False  # in this case we missed either the start or the end of the cell cycle.

    # generating the next cell cycles
    if range_timepoints[-1] < len(temp_obj.frames) - 1:  # if there is another generation of the current cell generate
        # that.
        ind += 1
        temp_nextgen_params = {'index': ind, 'cell_ind': temp_obj.index, 'cc_num': temp_gen + 1, 'parent': temp_obj}
        # if the current cell cycle is a bud, then we produce a daughter, otherwise we produce a mother.
        if temp_cycles[temp_ind].celltype == -1:
            temp_nextgen_params['celltype'] = 1
        else:
            temp_nextgen_params['celltype'] = 0
        temp_params['next_gen'] = ind
        temp_cycles.append(C.CellCycle(temp_nextgen_params))

    # Revise this to edit the daughter and generate a bud. We know this will be a bud, since we haven't assigned
    # lineages unless one has a previous timepoint to compare to.
    if range_timepoints[-1] in temp_obj.daughter_assignment_frames:  # if this cell has a daughter assigned at this
        # time point then generate the daughter
        ind += 1
        temp_daughter_params = {'index': ind, 'cell_ind': temp_obj.daughters[temp_obj.daughter_assignment_frames.index(range_timepoints[-1])],
                                'cc_num': 0, 'parent': temp_obj, 'celltype': -1}
        temp_params['bud'] = ind
        temp_cycles.append(C.CellCycle(temp_daughter_params))

    # determining the point of start in this cell cycle:
    temp1 = np.insert(temp, 0, 0)
    if np.sum(temp1[range_timepoints]==-1) > 1:  # testing the number of "start" events in this cell cycle
        temp_params['error'] = True
    elif np.sum(temp1[range_timepoints]==-1) == 1:  # testing whether a single start event was measured
        temp_params['start'] = range_timepoints[np.where(temp1[range_timepoints]==-1)]  # the index of the start event
        if temp1[temp_params['start']]!=-1:
            raise ValueError('Wrong indexing for Start')
    elif np.sum(temp1[range_timepoints]==-1) == 0:  # testing whether no start event was measured
        if temp_params['complete']:
            temp_params['error'] = True  # if we didn't measure a Start event, and yet caught the full cell cycle,
            # then something is wrong.

    # Updating the cell cycle with the appropriate info
    temp_cycles[temp_ind].update_cycle(temp_params)

    return temp_cycles

