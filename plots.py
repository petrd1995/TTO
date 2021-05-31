import numpy as np
import matplotlib.pyplot as plt


def nodecoords(node, array):
    '''
        Function that returns node's coordinates, based on the inputted node number (label)
    '''
    return array[node]


def plot_res(res, vec, all_nodes, n):
    '''
        Plotting results of the optimization. Each line corresponds to one bar. If single load case is analyzed, then red color means the bar is in tension, and blue color means compression. For multiple load cases this becomes there is no distinction by color.
    '''
    dim = np.shape(vec)[1]
    # btp - bars to plot
    btp = np.empty((np.count_nonzero(res[:, 0]), dim * 2 + 2))
    j = 0
    for i, bar in enumerate(res):
        if res[i, 0]:
            btp[j, 0:dim] = nodecoords(int(bar[2]), all_nodes)
            btp[j, dim:dim * 2] = nodecoords(int(bar[3]), all_nodes)
            btp[j, -2] = n[i]
            btp[j, -1] = bar[0]
            j += 1

    if dim > 2:
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(1, 1, 1, projection='3d')
        for i in range(len(btp)):
            if btp[i, -2] >= 0:
                ax1.plot([btp[i][0], btp[i][dim]], [btp[i][1], btp[i][dim + 1]],
                         [btp[i][2], btp[i][dim + 2]], 'r', linewidth=4 * btp[i, -1] / np.max(btp[:, -1]))
            else:
                ax1.plot([btp[i][0], btp[i][dim]], [btp[i][1], btp[i][dim + 1]],
                         [btp[i][2], btp[i][dim + 2]], 'b', linewidth=4 * btp[i, -1] / np.max(btp[:, -1]))
    else:
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(1, 1, 1)
        for i in range(len(btp)):
            if btp[i, -2] >= 0:
                ax1.plot([btp[i][0], btp[i][dim]], [btp[i][1], btp[i][dim + 1]],
                         'r', linewidth=4 * btp[i, -1] / np.max(btp[:, -1]))
            else:
                ax1.plot([btp[i][0], btp[i][dim]], [btp[i][1], btp[i][dim + 1]],
                         'b', linewidth=4 * btp[i, -1] / np.max(btp[:, -1]))

    ax1.set_xlabel('x [mm]')
    ax1.set_ylabel('y [mm]')
    if dim > 2:
        ax1.set_zlabel('z [mm]')
        ax1.set_zlabel('z [mm]')
        ax1.set_xlabel('x [mm]')
        ax1.set_ylabel('y [mm]')
    plt.show()


def plot_grid(nodes, bars, **plotbars):
    '''
        Function for plotting the grid or the so-called ground structure. Can be modified to show everything, or exclude bar plotting, node plotting etc. I use this function to look and verify the ground structure is correctly set up according to my vision. The other use is to show only node numbers and with that knowledge impose boundary conditions and loads to the corresponding nodes.
    '''
    dim = np.shape(nodes)[1]
    bbox_props = dict(boxstyle="round,pad=0.1,rounding_size=0.2", fc="white", ec="k")
    if dim > 2:
        x, y, z = nodes.T
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1, projection='3d')
        # ax1.scatter(x, y, z, s=50, marker="o", facecolor='white',
        #             edgecolor='k', depthshade=False)
        # if plotbars:
        for i in bars:
            ax1.plot(
                [nodecoords(bars[i[0], 1], nodes)[0], nodecoords(bars[i[0], 2], nodes)[0]],
                [nodecoords(bars[i[0], 1], nodes)[1], nodecoords(bars[i[0], 2], nodes)[1]],
                [nodecoords(bars[i[0], 1], nodes)[2], nodecoords(bars[i[0], 2], nodes)[2]],
                'k', linewidth=0.8)
        for i, txt in enumerate(range(len(nodes))):
            ax1.text(
                    nodes.T[0][i],
                    nodes.T[1][i],
                    nodes.T[2][i],
                    txt, bbox=bbox_props
                     )
    else:
        x, y = nodes.T
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1)
        ax1.scatter(x, y, s=50, marker="o", facecolor='white',
                    edgecolor='k')
        if plotbars:
            for i in bars:
                ax1.plot(
                    [nodecoords(bars[i[0], 1], nodes)[0], nodecoords(bars[i[0], 2], nodes)[0]],
                    [nodecoords(bars[i[0], 1], nodes)[1], nodecoords(bars[i[0], 2], nodes)[1]],
                    'k', linewidth=0.8)
    ax1.set_xlabel('x [mm]')
    ax1.set_ylabel('y [mm]')
    if dim > 2:
        ax1.set_zlabel('z [mm]')
        ax1.set_zlabel('z [mm]')
        ax1.set_xlabel('x [mm]')
        ax1.set_ylabel('y [mm]')
    plt.show()


def plot_bcs(bcs, f, nodes):
    '''
        Function that plots boundary conditions and forces, which can be used to quickly verify, that these were put into the correct nodes.
    '''
    bcs = np.array(bcs)
    f = np.array(f)
    dim = np.shape(bcs)[1] - 1
    if dim > 2:
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1, projection='3d')
        for en, i in enumerate(bcs.T[0]):
            x, y, z = nodecoords(int(i), nodes)
            ax1.scatter(x, y, z, marker='^', c='b')

        for i in f:
            mag = np.sqrt(i[1]**2 + i[2]**2 + i[3]**2)
            x, y, z = nodecoords(int(i[0]), nodes)
            coef = np.max(nodes.T[0]) / 3
            ax1.quiver(x, y, z, i[1] * coef / mag,
                       i[2] * coef / mag, i[3] * coef / mag,
                       color='darkgreen')
    else:
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1)
        for en, i in enumerate(bcs.T[0]):
            x, y = nodecoords(int(i), nodes)
            ax1.scatter(x, y, marker='^', c='b')

        for i in f:
            mag = np.sqrt(i[1]**2 + i[2]**2)
            x, y = nodecoords(int(i[0]), nodes)
            coef = np.max(nodes.T[0]) / 10
            ax1.quiver(x, y, i[1] * coef / mag,
                       i[2] * coef / mag,
                       color='darkgreen')

    ax1.set_xlabel('x [mm]')
    ax1.set_ylabel('y [mm]')
    if dim > 2:
        ax1.set_zlabel('z [mm]')
        ax1.set_zlabel('z [mm]')
        ax1.set_xlabel('x [mm]')
        ax1.set_ylabel('y [mm]')
    plt.show()

def plot_conv(iter, hist_eps):
    '''
        Function for plotting convergence of the optimization.
    '''
    fig, ax = plt.subplots()
    ax.plot(range(iter), hist_eps)
    # fig.title('Convergence')
    ax.set_xlabel('Iteration [-]')
    ax.set_ylabel('Epsilon [%]')

    plt.tight_layout()
    plt.show()
