import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import manual_image_toolkit as M
import os
import scipy
from scipy import stats

# this script takes as an input the output of "analyze_whi5_distribution.py", i.e. a list of cells with nuclear whi5
# residencies.
# the underlying assumption in this script is that the cell segmentation and tracking is reliable.

# assumes the variable names daughters, parent

base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
fluor_name = '/180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w2515 laser 20_'
drange = 65535.0

# analyzing manually annotated dataset for scene 1
distance_cutoff = 30  # cutoff for centroid distance above which mother-daughter pairings cannot be assigned
num_divisions = 0
num_nuclear_events = 0
frame_nums = [31, 50, 50, 50, 50, 50, 50]
for scene in range(1, 8):
    print 'scene', scene
    directory = base_path+expt_path+'/scene_{0}/outputs'.format(scene)
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
    for frame_num in range(1, frame_nums[scene-1]):  # we go through frame by frame, starting with frame # 2 since we need a previous
        # timepoint to compare it to.
        print 'frame', frame_num
        assigned_inds = [[], []]
        pending_assignment.append([])
        if frame_num>1:

            temp2 = [(frame_num in temp1) for temp1 in frame_list]
            update_list = [i for i, e in enumerate(temp2) if e != 0]  # list of cells in the current frame
            # we now consider only cells for whom Whi5 is in the nucleus. This happens coincidentally for
            # mother and daughter cells, and we use spatial proximity to determine which daughter belongs to which
            # mother.
            update_list1 = [i for i in update_list if c[i].nuclear_whi5[frame_num-c[i].frames[0]] == 1]
            temp_new = [[], [], []]
            temp_started_G1 = [[], [], []]
            for ind1 in update_list1:
                temp_ind = frame_num-c[ind1].frames[0]  # cell age
                if temp_ind == 0:  # if this cell was segmented for the first time in this frame then we store that,
                    # but we don't do much with it for now.
                    print 'Newborn cell {0} with nuclear Whi5 found in timepoint {1}'.format(c[ind1].index, frame_num)
                    temp_new[0].append(ind1)
                    temp_new[1].append(temp_ind)
                    # print c[ind1].position
                    temp_new[2].append(c[ind1].position[temp_ind])
                    num_nuclear_events += 1
                elif c[ind1].nuclear_whi5[temp_ind-1] == 0:
                    # if temp_ind-2==0 or c[ind1].nuclear_whi5[temp_ind-2] == 0:  # we require that either the cell only
                    # started being
                    temp_started_G1[0].append(ind1)  # these are the cells we have to assign
                    temp_started_G1[1].append(temp_ind)
                    temp_started_G1[2].append(c[ind1].position[temp_ind])
                    num_nuclear_events += 1

            temp_mothers = np.array([(len(c[temp_ind1].daughters) > 0) for temp_ind1 in temp_started_G1[0]])
            # tracking whether these cells have already been mothers.
            to_remove = [[], []]
            for ind1 in range(len(temp_started_G1[0])):
                # evaluating distances
                d = 1000*np.ones(len(temp_started_G1[0]))
                for ind2 in range(len(temp_started_G1[0])):
                    if ind1 != ind2:
                        # print d, temp_started_G1
                        d[ind2] = np.linalg.norm(temp_started_G1[2][ind1]-temp_started_G1[2][ind2])

                if temp_mothers[ind1]:  # if this cell is a mother then find the cells which haven't previously been
                    # classified as mothers.
                    temp_inds = np.nonzero(temp_mothers == 0)
                    temp_ind2 = np.argsort(d[temp_inds])  # ordering the indices of temp_inds based on their distance
                    # print len(temp_inds), temp_ind2, temp_inds[0]
                    if len(temp_inds[0]) > 0:
                        if d[temp_inds[0][temp_ind2[0]]] < distance_cutoff:
                            if sum(d[temp_inds] < distance_cutoff) > 1:  # if there is another cell with low distance in
                                # the same frame then we assign this
                                pending_assignment[-1].append(temp_started_G1[0][ind1])
                                c = M.assign_troublesome_pair(c, temp_started_G1[0][ind1], frame_num)
                            else:
                                assigned_inds[0].append(temp_started_G1[0][ind1])
                                assigned_inds[1].append(temp_started_G1[0][temp_inds[0][temp_ind2[0]]])
                                c = C.assign_md_pair(c, mother_ind=temp_started_G1[0][ind1],
                                                     daughter_ind=temp_started_G1[0][temp_inds[0][temp_ind2[0]]],
                                                     temp_frame_num=frame_num)
                                num_divisions+=1
                    else:
                        pending_assignment[-1].append(temp_started_G1[0][ind1])  # if we can't assign this well,
                        # then we will revisit it in the next timeframe
                else:
                    temp_ind2 = np.argsort(d)
                    if d[temp_ind2[0]] < distance_cutoff:
                        if sum(d < distance_cutoff) > 1:
                            pending_assignment[-1].append(temp_started_G1[0][ind1])
                            c = M.assign_troublesome_pair(c, temp_started_G1[0][ind1], frame_num)
                        else:
                            temp_inds1 = [temp_started_G1[0][ind1], temp_started_G1[0][temp_ind2[0]]]
                            # print temp_inds1
                            temp_val = [c[temp_inds1[0]].ellipse_volume[frame_num-c[temp_inds1[0]].frames[0]],
                                               c[temp_inds1[1]].ellipse_volume[frame_num-c[temp_inds1[1]].frames[0]]]
                            # print temp_val
                            m_ind = np.argmax(temp_val)
                            # find which is the mother cell based on which has the larger cell size.
                            assigned_inds[0].append(temp_inds1[m_ind])
                            assigned_inds[1].append(temp_inds1[m_ind-1])
                            c = C.assign_md_pair(c, mother_ind=temp_inds1[m_ind], daughter_ind=temp_inds1[m_ind-1],
                                                 temp_frame_num=frame_num)
                            num_divisions += 1
                    else:
                        pending_assignment[-1].append(temp_started_G1[0][ind1])  # if we can't assign this well,
                        # then we will plan to revisit it in the next timeframe

                # we take any cell that is pending assignment in the current frame, and compare to the previous frame to
                # see if there is an appropriate assignment there.
                for i0 in range(len(pending_assignment[-1])):
                    coord1 = c[pending_assignment[-1][i0]].position[frame_num - c[pending_assignment[-1][i0]].frames[0]]
                    d = []
                    for i1 in range(len(pending_assignment[-2])):
                        # compare the positions in the current frame and the previous frame
                        coord2 = c[pending_assignment[-2][i1]].position[frame_num - 1 -
                                                                        c[pending_assignment[-2][i1]].frames[0]]
                        d.append(np.linalg.norm(coord1 - coord2))
                    temp_ind2 = np.argsort(d)
                    if len(temp_ind2) != 0:
                        if d[temp_ind2[0]] < distance_cutoff:
                            if len(d)>1:
                                # print d, np.asarray(d) < distance_cutoff
                                if sum(np.asarray(d) < distance_cutoff) == 1:
                                    temp_inds1 = [pending_assignment[-1][i0], pending_assignment[-2][temp_ind2[0]]]
                                    # print temp_inds1
                                    temp_val = [c[temp_inds1[0]].ellipse_volume[frame_num - c[temp_inds1[0]].frames[0]],
                                                c[temp_inds1[1]].ellipse_volume[frame_num -1 - c[temp_inds1[1]].frames[0]]]
                                    m_ind = np.argmax(temp_val)
                                    # find which is the mother cell based on which has the larger cell size.
                                    assigned_inds[0].append(temp_inds1[m_ind])
                                    assigned_inds[1].append(temp_inds1[m_ind - 1])
                                    c = C.assign_md_pair(c, mother_ind=temp_inds1[m_ind],
                                                         daughter_ind=temp_inds1[m_ind - 1], temp_frame_num=frame_num-1)
                                    num_divisions += 2
                                    # we assign these cells as both being present in the previous frame, since that's when
                                    # the nuclear localization event must have occurred for both
                                    to_remove[0].append(i0)
                                    to_remove[1].append(temp_ind2[0])
                            else:
                                temp_inds1 = [pending_assignment[-1][i0], pending_assignment[-2][temp_ind2[0]]]
                                # print temp_inds1
                                temp_val = [c[temp_inds1[0]].ellipse_volume[frame_num - c[temp_inds1[0]].frames[0]],
                                            c[temp_inds1[1]].ellipse_volume[frame_num - 1 - c[temp_inds1[1]].frames[0]]]
                                m_ind = np.argmax(temp_val)
                                # find which is the mother cell based on which has the larger cell size.
                                assigned_inds[0].append(temp_inds1[m_ind])
                                assigned_inds[1].append(temp_inds1[m_ind - 1])
                                c = C.assign_md_pair(c, mother_ind=temp_inds1[m_ind],
                                                     daughter_ind=temp_inds1[m_ind - 1], temp_frame_num=frame_num - 1)
                                num_divisions += 2
                                # we assign these cells as both being present in the previous frame, since that's when
                                # the nuclear localization event must have occurred for both. They will appear in the
                                # output images as being tracked in the current frame.
                                to_remove[0].append(i0)
                                to_remove[1].append(temp_ind2[0])
                # removing the newly tracked items from pending_assignment. Cells from the previous timestep will still
                # appear as not being tracked in the previous frame.
            pending_assignment[-1] = list(set(pending_assignment[-1])-set(to_remove[0]))
            pending_assignment[-2] = list(set(pending_assignment[-2]) - set(to_remove[1]))
        # producing a figure to track the result of this so far.
        # if frame_num ==50 and scene==4:
        #     for ind1 in assigned_inds[1]:
        #         print ind1, len(c[ind1].position), frame_num-c[ind1].frames[0], c[ind1].nuclear_whi5
        #         print c[ind1].position
        temp_im = io.imread(
            base_path + expt_path + fluor_name + 's{0}_t{1}.TIF'.format(str(scene), str(frame_num))) / drange
        temp_im1 = np.amax(temp_im, axis=0) / np.amax(temp_im)  # scaling for visualization purposes.
        temp_im1 *= outlines[frame_num - 1, :, :] == 0
        temp_im1 += outlines[frame_num - 1, :, :].astype('uint16')
        m_cell_coords = [c[ind1].position[frame_num - c[ind1].frames[0]] for ind1 in assigned_inds[0]]
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
        if len(m_cell_coords)>0:
            temp_x, temp_y = zip(*m_cell_coords)
            plt.plot(temp_x, temp_y, '.', markersize=6, color='w', linestyle='None')
        if len(d_cell_coords)>0:
            temp_x, temp_y = zip(*d_cell_coords)
            plt.plot(temp_x, temp_y, '.', markersize=4, color='k', linestyle='None')
        if len(bad_cell_coords)>0:
            temp_x, temp_y = zip(*bad_cell_coords)
            plt.plot(temp_x, temp_y, '.', markersize=2, color='r', linestyle='None')
        fig.subplots_adjust(bottom=0)
        fig.subplots_adjust(top=1)
        fig.subplots_adjust(right=1)
        fig.subplots_adjust(left=0)
        fig.savefig(directory + '/images/automated_lineage_assignments_frame_{0}.tif'.format(frame_num))
    C.save_object(c, directory + '/cells_scene_{0}_v3.pkl'.format(scene))
    # print 'I got here', directory + '/images/cells_scene_{0}_v3.pkl'.format(scene)

print num_divisions, num_nuclear_events



