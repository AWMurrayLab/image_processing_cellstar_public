import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage import io
import cPickle as pickle
import custom_image_toolkit as C
import os
import scipy
from scipy import stats

# initial definitions
fluor_name = '/180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w2515 laser 20_'
drange = 65535.0
base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
# analyzing manually annotated dataset for scene 1
scene = 1
directory = base_path+expt_path+'/scene_{0}/outputs'.format(scene)
with open(directory + '/cells_scene_{0}_v1.pkl'.format(scene), 'rb') as input:
    c = pickle.load(input)
if not os.path.exists(directory + '/plots'):
            os.makedirs(directory + '/plots')
temp1 = [[],[],[]]  # G1 cells
temp2 = [[],[],[]]  # G2 cells
for obj in c:
    c1 = [obj.fluor_vals[i0] for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0]]
    c2 = [obj.index for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0]]
    c3 = [obj.frames[i0] for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0]]
    c4 = [obj.fluor_vals[i0] for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0] == 0]
    c5 = [obj.index for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0] == 0]
    c6 = [i0 for i0 in range(len(obj.frames)) if obj.nuclear_whi5[i0] == 0]
    temp1[0] += c1
    temp1[1] += c2
    temp1[2] += c3
    temp2[0] += c4
    temp2[1] += c5
    temp2[2] += c6

lower_bound = 0.01

v1 = [scipy.stats.skew(vals) for vals in temp1[0]]
v2 = [scipy.stats.skew(vals) for vals in temp2[0]]
make_dist_plots_1 = 0
if make_dist_plots_1:
    import seaborn as sns
    fig = plt.figure(figsize=[5,5])
    sns.distplot(v1, label='Nuclear')
    sns.distplot(v2, label='Cytoplasmic')
    plt.legend()
    plt.xlabel('Skewness')
    plt.title('Skewness in manually annotated cell populations')
    fig.savefig(directory+'/plots/skewness')
    del fig
fig = plt.figure(figsize=[5,5])
kde = scipy.stats.gaussian_kde(v1)
kde1 = scipy.stats.gaussian_kde(v2)
xpoints = np.linspace(-2.0,5.0, 100)
vals = [kde.integrate_box_1d(low=-10, high=lim) for lim in xpoints]
sel = [i for i, e in enumerate(vals) if e>lower_bound]
skew_cutoff = xpoints[sel[0]]
plt.axvline(x=xpoints[sel[0]], label='99% bound = {0}'.format(str(np.round(skew_cutoff, decimals=3))))
print kde1.integrate_box_1d(low=skew_cutoff, high=10.0)
plt.plot(xpoints, kde.evaluate(xpoints), label='Nuclear')
plt.plot(xpoints, kde1.evaluate(xpoints), label='Cytoplasmic')
plt.plot(xpoints, vals, label='Nuclear CDF')
plt.legend()
fig.savefig(directory+'/plots/skewness_bounding')

# making initial plots of single cell distributions
make_dist_plots = 0
if make_dist_plots:
    import seaborn as sns
    inds = [np.random.randint(0, len(temp2[0]), size=10), np.random.randint(0, len(temp1[0]), size=10)]
    labels = ['Whi5 non-nuclear localized', 'Whi5 nuclear localized']
    for i0 in range(2):
        for ind1 in inds[i0]:
            # print ind
            fig=plt.figure(figsize=[5, 5])
            plt.xlabel('Whi5 fluorescence intensity')
            plt.title(labels[i0])
            if i0 == 1:
                sns.distplot(temp1[0][ind1], label = 'cell # {0}, frame # {1}'.format(temp1[1][ind1], temp1[2][ind1]))
                plt.legend()
                fig.savefig(directory + '/plots/G1_{0}_cell_{1}_frame_{2}'.format(i0, temp1[1][ind1], temp1[2][ind1]))
            else:
                sns.distplot(temp2[0][ind1], label='cell # {0}, frame # {1}'.format(temp2[1][ind1], temp2[2][ind1]))
                plt.legend()
                fig.savefig(directory + '/plots/G1_{0}_cell_{1}_frame_{2}'.format(i0, temp2[1][ind1], temp2[2][ind1]))

    for i0 in range(2):
        for ind1 in inds[i0]:
            # print ind
            fig=plt.figure(figsize=[5, 5])
            plt.xlabel('Whi5 fluorescence intensity')
            plt.title(labels[i0])
            if i0 == 1:
                sns.distplot(temp1[0][ind1]/np.mean(temp1[0][ind1]), label = 'cell # {0}, frame # {1}'.format(temp1[1][ind1], temp1[2][ind1]))
                plt.legend()
                fig.savefig(directory + '/plots/norm_G1_{0}_cell_{1}_frame_{2}'.format(i0, temp1[1][ind1], temp1[2][ind1]))
            else:
                sns.distplot(temp2[0][ind1]/np.mean(temp2[0][ind1]), label='cell # {0}, frame # {1}'.format(temp2[1][ind1], temp2[2][ind1]))
                plt.legend()
                fig.savefig(directory + '/plots/norm_G1_{0}_cell_{1}_frame_{2}'.format(i0, temp2[1][ind1], temp2[2][ind1]))

# Automatically classifying cells based on skewness cutoff
frames = [31, 51, 51, 51, 51, 51, 51]
automate_analysis = 1
if automate_analysis:
    for scene in range(1, 8):
        directory = base_path + expt_path + '/scene_{0}/outputs'.format(scene)
        with open(directory + '/cells_fl_lab_scene_{0}.pkl'.format(scene), 'rb') as input:
            c = pickle.load(input)
        c = C.automated_whi5_assignment_1(c, skew_cutoff)
        C.save_object(c, directory + '/cells_scene_{0}_v2.pkl'.format(scene))
        # producing plots to ensure that this tracking is adequate
        outlines = np.load(directory + '/cell_outlines_scene_{0}.npy'.format(scene))
        frame_list = [obj.frames for obj in c]

        # organizing figure data
        for frame_num in range(1, frames[scene-1]):
            # checking which cells to look at for this timepoint
            temp2 = [(frame_num in temp1) for temp1 in frame_list]
            update_list = [i for i, e in enumerate(temp2) if e != 0]
            # loading image
            temp_im = io.imread(
                base_path + expt_path + fluor_name + 's{0}_t{1}.TIF'.format(str(scene), str(frame_num))) / drange
            temp_im1 = np.amax(temp_im, axis=0) / np.amax(temp_im)  # scaling for visualization purposes.
            temp_im1 *= outlines[frame_num - 1, :, :] == 0
            temp_im1 += outlines[frame_num - 1, :, :].astype('uint16')
            # print c[ind].frames
            # print [frame_num - c[ind1].frames[0] for ind1 in update_list if (c[ind1].nuclear_whi5[int(frame_num - c[ind1].frames[0])])]
            temp_cell_coords = [c[ind1].position[frame_num - c[ind1].frames[0]] for ind1 in update_list if (c[ind1].nuclear_whi5[
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
            temp_x, temp_y = zip(*temp_cell_coords)
            plt.plot(temp_x, temp_y, '.', markersize=2, color='r', linestyle='None')

            fig.subplots_adjust(bottom=0)
            fig.subplots_adjust(top=1)
            fig.subplots_adjust(right=1)
            fig.subplots_adjust(left=0)
            fig.savefig(directory + '/images/automated_whi5_assignments_frame_{0}.tif'.format(frame_num))
            # exit()

