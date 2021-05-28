import numpy as np

'''
    This module includes two functions, which create arrays used for defining boundary conditions and forces, which have too many entries to enter them by hand, thus doing so programmatically.
'''

def createBCs(node_start, node_end, dim):
    '''
        This function takes starting node's number, ending node's number and the dimension of the problem (2D, 3D) as an input and creates a vector of boundary conditions ranging from the starting node to the end node. (Code too difficult - can be easily simplified)
    '''
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
    '''
        Similarly to createBCs this function takes starting node, ending node and dimension of the problem as an input but also the actual force to be inserted in each node. For simplification reasons the force will be the same in all nodes - used as multiple load cases setup, inserting different force in each node would be possible by slightly altering the code (not necessary right now).
    '''
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
