import TTO
import numpy as np
import emailnotify

x0 = 100
y0 = 100
z0 = 0
nx = 20
ny = 10
nz = 2
E = 2.1e5
r0 = 10
ratio = 0.5
# Vol0 = x0*y0*z0
Vol0 = x0*y0
Ro = 1
kon = 0.1
bcs = np.array([[0, 1, 1], [1, 1, 1],[2, 1, 1], [3, 1, 1], [4, 1, 1], [5, 1, 1], 
               [6, 1, 1], [7, 1, 1], [8, 1, 1], [9, 1, 1]])
f = np.array([[194, 0, -1000]])

example = TTO.Truss('example', x0, y0, z0, nx, ny, nz, bcs, f, E, r0, Vol0, ratio, Ro, kon)

example.default_setup()
# example.rem_long_bars(2)
example.opt()
example.plot('res')
# example.plot('bcs')
# example.plot('grid')
# emailnotify.notify()
# example.plot('conv')
# example.out()

# x0 = 100
# y0 = 100
# z0 = 100
# nx = 2
# ny = 2
# nz = 2
# E = 2.1e5
# r0 = 5
# ratio = 0.5
# Vol0 = x0*y0*z0
# Ro = 1
# kon = 1
# bcs = np.array([[0, 1, 1, 1], [1, 1, 1, 1], [2, 1, 1, 1], [3, 1, 1, 1]])
# f = np.array([[4, 0, -1000, -600]])

# example = TTO.Truss('example', x0, y0, z0, nx, ny, nz,
#                     bcs, f, E, r0, Vol0, ratio, Ro, kon)

# example.default_setup()
# example.opt()
# # example.plot('bcs')
# # example.plot('grid')
# emailnotify.notify()
# # example.plot('conv')
# # example.out()
