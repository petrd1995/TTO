import numpy as np
import itertools
import circlecoords


def create_grid(nnodes, x0, y0, z0, nx, ny, nz):

    x = np.linspace(20, x0, nx)
    y = np.linspace(20, y0, ny)

    if z0:
        z = np.linspace(20, z0, nz)
        all_nodes = np.empty((nnodes, 3))
        get_node = itertools.product(x, y, z)

    else:
        all_nodes = np.empty((nnodes, 2))
        get_node = itertools.product(x, y)

    for i, el in enumerate(get_node):
        all_nodes[i] = el

    return all_nodes


def add_one_node(coords, all_nodes):

    # coords musi byt numpy array
    all_nodes = np.block([[all_nodes], [coords]])
    num_nodes = np.shape(all_nodes)[0]

    return all_nodes, num_nodes


def add_circle(dim, all_nodes, radius, plane, n_points, *origin):

    if dim > 2:
        x, y, z = circlecoords.circle(
            dim, plane, radius, n_points, *origin)
        newcoords = np.column_stack((x, np.column_stack((y, z))))
        return add_one_node(newcoords, all_nodes)
    else:
        x, y = circlecoords.circle(dim, plane, radius, n_points, *origin)
        newcoords = np.column_stack((x, y))
        return add_one_node(newcoords, all_nodes)


def node_coords(node_num, all_nodes, **kwargs):

    if 'show' in kwargs:
        print(
            f'the coordinates for node {node_num} are: {all_nodes[node_num]}')
    return all_nodes[node_num]


def rem_node(node_num, all_nodes, num_nodes):

    all_nodes = np.delete(all_nodes, node_num, axis=0)
    num_nodes = np.shape(all_nodes)[0]

    return all_nodes, num_nodes


def create_bars(num_nodes):

    num_bars = int(num_nodes * (num_nodes - 1) / 2)
    node_counter = np.arange(num_nodes)
    bars = np.empty((num_bars, 3), dtype=int)
    comb = itertools.combinations(node_counter, 2)

    for q, i in enumerate(comb):
        bars[q, ] = int(q), *i

    return bars


def rem_bars(bar_num, bars):

    bars = np.delete(bars, bar_num, axis=0)
    num_bars = len(bars)
    bars.T[0] = np.arange(num_bars)

    return bars, num_bars


def rem_long_bars_length(length, num_bars, vec, lens, bars):
    '''remove bars of prescribed length multiple of length of diagonal bar'''

    a = []
    for i in range(num_bars):
        if (lens[i] > length) or (lens[i] < 0.5 * length):
            a.append(i)
    rem_bars(a, bars)  # odstranění prutu/ů?
    vec = np.delete(vec, a, axis=0)
    lens = np.delete(lens, a, axis=0)

    return bars, vec, lens


def vec_len(num_bars, bars, dim, all_nodes):

    lens = np.zeros((num_bars, 1))
    if dim > 2:
        vec = np.zeros((num_bars, 3))
        dim3D = True
    else:
        vec = np.zeros((num_bars, 2))
    for bar in bars:

        start = node_coords(bars[bar[0], 1], all_nodes)
        end = node_coords(bars[bar[0], 2], all_nodes)

        lens[bar[0]] = np.sqrt((end - start).dot(end - start))
        vec[bar[0], 0:2] = ((end - start) / lens[bar[0]])[0:2]
        if dim3D:
            vec[bar[0], 2] = ((end - start) / lens[bar[0]])[2]

    return vec, lens


def matK(E, num_bars, bars, rB, cB, vec, Avec):

    K = np.zeros((int(cB * rB), int(cB * rB)))
    for i in range(num_bars):
        slc1 = cB * bars[i][1]
        slc1end = slc1 + cB
        slc2 = cB * bars[i][2]
        slc2end = slc2 + cB

        K[slc1:slc1end,
            slc1:slc1end] += np.outer(vec[i], vec[i]) * E * Avec[i] / len[i]

        K[slc2:slc2end,
            slc2:slc2end] += np.outer(vec[i], vec[i]) * E * Avec[i] / len[i]

        K[slc1:slc1end, slc2:slc2end] = - \
            np.outer(vec[i], vec[i]) * E * Avec[i] / len[i]

        K[slc2:slc2end, slc1:slc1end] = - \
            np.outer(vec[i], vec[i]) * E * Avec[i] / len[i]

    return K


def forces(F, all_nodes):

    rB, cB = np.shape(all_nodes)
    f = np.zeros((rB * cB, 1))
    for i in range(len(F)):
        slc = cB * int(F[i][0])
        f[slc] = F[i][1]
        f[slc + 1] = F[i][2]
        if cB > 2:
            f[slc + 2] = F[i][3]
    return f


def boundary(bc, K, all_nodes):

    rB, cB = np.shape(all_nodes)
    for i in range(len(bc)):
        slc = cB * int(bc[i, 0])
        if bc[i, 1]:
            K[:, slc] = 0
            K[slc, :] = 0
            K[slc, slc] = 1
        if bc[i, 2]:
            K[:, slc + 1] = 0
            K[slc + 1, :] = 0
            K[slc + 1, slc + 1] = 1
        if cB > 2:
            if bc[i, 3]:
                K[:, slc + 2] = 0
                K[slc + 2, :] = 0
                K[slc + 2, slc + 2] = 1

    return K


def zerocrosssection(K):

    p1 = np.diag(K)
    p1 = np.array(p1)

    findzeros = np.where(p1 < 1e-5)
    for i in findzeros:
        p1[i] = 1
    np.fill_diagonal(K, 0)
    K = K + np.diag(p1)
    return K


def createBCs(node_start, node_end, dim):
    # přidá vetktnutí, dalo by se upravit na posuvnou vazbu změnou 1 na nulu u příslušného styčníku v příslušném směru
    if dim == 3:
        bcs = np.empty((node_end + 1 - node_start, 4))
        for en, i in enumerate(bcs):
            bcs[en][0], bcs[en][1], bcs[en][2], bcs[en][3] = en + \
                node_start, 1, 1, 1
    else:
        bcs = np.empty((node_end + 1 - node_start, 3))
        for en, i in enumerate(bcs):
            bcs[en][0], bcs[en][1], bcs[en][2] = en + node_start, 1, 1
    return bcs


def createForces(node_start, node_end, dim, force):
    # přidá síly do daných uzlů, všechny síly jsou pro jednoduchost stejné
    if dim == 3:
        frcs = np.empty((node_end + 1 - node_start, 4))
        for en, i in enumerate(frcs):
            frcs[en][0], frcs[en][1], frcs[en][2], frcs[en][3] = en + \
                node_start, force[0], force[1], force[2]
    else:
        frcs = np.empty((node_end + 1 - node_start, 3))
        for en, i in enumerate(frcs):
            frcs[en][0], frcs[en][1], frcs[en][2] = en + \
                node_start, force[0], force[1]
    return frcs
