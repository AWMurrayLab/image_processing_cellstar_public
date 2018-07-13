import numpy as np
import matplotlib.pyplot as plt
from skimage import io
import os


def onclick(event):
    ix, iy = event.xdata, event.ydata
    print 'x = %d, y = %d'%(
        ix, iy)
    x_pts.append(ix)
    y_pts.append(iy)
    line.set_xdata(x_pts)
    line.set_ydata(y_pts)
    fig.canvas.draw()
    global coords
    # coords.append(np.array([np.floor(ix), np.floor(iy)]).astype('int'))
    coords.append((np.floor(ix), np.floor(iy)))
    return coords


def keypress(event):
    # global val
    val = event.key

    if val == 'f':
        # close image
        fig.canvas.mpl_disconnect(cid)
        fig.canvas.mpl_disconnect(cid1)
        plt.close(fig)
        if not os.path.exists(directory+'/fl_loc_centres'):
            os.makedirs(directory+'/fl_loc_centres')
        x_vals, y_vals = zip(*coords)
        temp = np.zeros([len(x_vals),2])
        temp[:, 0] = x_vals[:]
        temp[:, 1] = y_vals[:]
        np.save(directory+'/fl_loc_centres/'+'scene_{0}_frame_{1}'.format(scene, frame_num), temp)
        print 'completed scene {0}, frame {1}, coordinates:'.format(scene, frame_num), temp
    elif val == 'b':
        temp = coords.pop()
        x_pts.pop()
        y_pts.pop()
        line.set_xdata(x_pts)
        line.set_ydatba(y_pts)
        fig.canvas.draw()
        print temp
    # elif val == 'b':
    #     temp = coords.pop()
    #     x_pts.pop()
    #     y_pts.pop()
    #     line.set_xdata(x_pts)
    #     line.set_ydata(y_pts)
    #     fig.canvas.draw()
    #     print temp
    # elif val == 'q':
    return coords

base_path, expt_path = '/home/felix/Documents/image_analysis', '/180531_dye_mix_expt/timelapse'
fluor_name = '/180531_60X_20lp_yFB43_yFB29dyed_yfp_10min_w2515 laser 20_'
drange = 65535.0
coords = []
global ax, temp_im, fig, scene, x_pts, y_pts, directory, frame_num
for scene in range(1, 4):
    directory = base_path+expt_path+'/scene_{0}/outputs'.format(scene)
    outlines = np.load(directory+'/cell_outlines_scene_{0}.npy'.format(scene))
    coord_list = []
    for frame_num in range(1, 31):
        coord_list.append([])
    # print np.amax(outlines)
    # exit()
    for frame_num in range(1, 31):
        temp_im = io.imread(base_path+expt_path+fluor_name+'s{0}_t{1}.TIF'.format(str(scene), str(frame_num)))/drange
        # print temp_im.dtype
        # exit()
        temp_im1 = np.amax(temp_im, axis=0)/np.amax(temp_im)
        # print np.amax(temp_im1)
        # exit()
        temp_im1*= outlines[frame_num-1, :, :]==0
        temp_im1+= outlines[frame_num-1, :, :].astype('uint16')
        coords = []
        # filename = '/scene_{0}/outputs/images/frame_{1}.tif'.format(scene, frame)
        # temp_im = io.imread(base_path+expt_path+filename)
        fig = plt.figure(figsize=[10, 10], frameon=False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        ax.imshow(temp_im1)
        # variable for the points
        x_pts = []
        y_pts = []
        line, = ax.plot(x_pts, y_pts, marker="o", linestyle='None', color='r')
        # ax.imshow(outlines[frame_num-1, :, :])
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        cid1 = fig.canvas.mpl_connect('key_press_event', keypress)
        plt.show(fig)
        # print 'completed scene {0}, frame {1}, coordinates:'.format(scene, frame_num), coords
        coord_list.append(coords)
    coord_vals = np.zeros([sum([len(obj) for obj in coord_list]), 3])
    ind=0
    ind1 = 1
    for obj in coord_list:
        x,y = zip(*obj)
        coord_vals[ind:ind + len(obj), 0] = ind1
        coord_vals[ind:ind + len(obj), 1] = x[:]
        coord_vals[ind:ind + len(obj), 2] = y[:]
        ind1+=1
        ind+=len(obj)
    print len(coord_list), sum([len(obj) for obj in coord_list]), coord_vals
    np.save(directory+'/fl_loc_cells_scene_{0}.npy'.format(scene), coord_vals)
