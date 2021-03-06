'''
    testing file
'''

import TTO
import numpy as np
import emailnotify

x0 = 100
y0 = 100
z0 = 100
nx = 2
ny = 2
nz = 2
E = 2.1e5
r0 = 5
ratio = 0.5
Vol0 = x0*y0*z0
Ro = 1
kon = 1
bcs = np.array([[0, 1, 1, 1], [1, 1, 1, 1], [2, 1, 1, 1], [3, 1, 1, 1]])
f = np.array([[5, 0, -1000, -600]])

nodes = np.array([[0, 0, 0],
                  [10, 0, 0],
                  [0, 10, 0],
                  [0, 0, 10],
                  [0, 30, 30],
                  [50, 0, 50],
                  [50, 20, 50],
                  [50, 50, 50]])


example = TTO.Truss('example', x0, y0, z0, nx, ny, nz,
                    bcs, f, E, r0, Vol0, ratio, Ro, kon)

example.grid_from_list(nodes)
example.create_bars()
example.vec_len()
example.opt()
example.plot('res')
# example.plot('bcs')
# example.plot('grid')
