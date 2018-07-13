import numpy as np
import matplotlib.pyplot as plt


def onclick(event):
    ix, iy = event.xdata, event.ydata
    print 'x = %d, y = %d'%(
        ix, iy)

    global coords
    # coords.append(np.array([np.floor(ix), np.floor(iy)]).astype('int'))
    coords.append((np.floor(ix), np.floor(iy)))
    return coords


def keypress(event):
    global val
    val = event.key

    if val == 'f':
        # close image
        fig.canvas.mpl_disconnect(cid)
        fig.canvas.mpl_disconnect(cid1)
        plt.close(fig)

    elif val == 'b':
        temp = coords.pop()
        print temp
    return coords
