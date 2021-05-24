import new_kostka
import numpy as np
from create_truss import createBCs, createForces
import emailnotify

x0 = 280
y0 = 280
z0 = 280
nx = 11
ny = 11
nz = 11
E = 2.1e5
# E = 1
A0 = 10
ratio = 0.4
Ro = 1
kon = 10
Vol0 = 20640708 

bcnodes = 13
fnodes = 18
f_array = np.array([0, 0, 1])

nnodes = nx * ny * nz
bcnode_start = nnodes
bcnode_end = nnodes + bcnodes - 1
fnode_start = bcnode_end + 1
after_nodes = nnodes + fnodes + bcnodes - 1

bcs = createBCs(bcnode_start, bcnode_end, 3)
f = createForces(fnode_start, after_nodes, 3, f_array)
# print(f)
kostka = new_kostka.Truss('kostkaA', x0, y0, z0, nx, ny, nz, bcs, f, E, A0, Vol0, ratio, Ro, kon)

nds_to_add = np.array([[32.24,32.24,0],
                        [38.4,140,0],
                        [32.24,247.76,0],
                        [89.2,140,0],
                        [140,38.4,0],
                        [140,89.2,0],
                        [140,140,0],
                        [140,190.8,0],
                        [140,241.6,0],
                        [190.8,140,0],
                        [247.76,32.24,0],
                        [241.6,140,0],
                        [247.76,247.76,0],

                        [-20,77,61],
                        [-20,139,70],
                        [-20,144.5,87.5],
                        [-20,133.5,192.5],
                        [-20,141,210],
                        [-20,203,219],
                        [203,61,300],
                        [141,70,300],
                        [135,87.5,300],
                        [77,219,300],
                        [139,210,300],
                        [144.5,192.5,300],
                        [300,192.5,135.5],
                        [300,210,141],
                        [300,219,203],
                        [300,61,77],
                        [300,70,139],
                        [300,87.5,144.5]])

kostka.create_grid()

kostka.add_one_node(nds_to_add)

kostka.create_bars()
kostka.vec_len()
kostka.rem_long_bars(1)


kostka.plot('grid')
# kostka.plot('bcs')

# kostka.opt()
# kostka.out()

## emailnotify.notify()
# kostka.plot('res')
## kostka.plot('bcs')