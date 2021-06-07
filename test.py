# import TTO
# import numpy as np
# import emailnotify

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
# bcs = np.array([[0, 1, 1, 1], [1, 1, 1, 1],[2, 1, 1, 1], [3, 1, 1, 1]]) 
# f = np.array([[5, 0, -1000, -600]])

# nodes = np.array([[0, 0, 0],
#                   [10, 0, 0], 
#                   [0, 10, 0], 
#                   [0, 0, 10], 
#                   [0, 30, 30], 
#                   [50, 0, 50],
#                   [50, 20, 50],
#                   [50, 50, 50]])


# example = TTO.Truss('example', x0, y0, z0, nx, ny, nz, bcs, f, E, r0, Vol0, ratio, Ro, kon)

# example.grid_from_list(nodes)
# example.create_bars()
# example.vec_len()
# example.opt()
# example.plot('res')
# # example.plot('bcs')
# # example.plot('grid')

import nlopt
from numpy import *


def myfunc(x, grad):
    if grad.size > 0:
        grad[0] = 0.0
        grad[1] = 0.5 / sqrt(x[1])
    return sqrt(x[1])


def myconstraint(x, grad, a, b):
    if grad.size > 0:
        grad[0] = 3 * a * (a*x[0] + b)**2
        grad[1] = -1.0
    return (a*x[0] + b)**3 - x[1]


opt = nlopt.opt(nlopt.LD_MMA, 2)
opt.set_lower_bounds([-float('inf'), 0])
opt.set_min_objective(myfunc)
opt.add_inequality_constraint(
    lambda x, grad: myconstraint(x, grad, 2, 0), 1e-8)
opt.add_inequality_constraint(
    lambda x, grad: myconstraint(x, grad, -1, 1), 1e-8)
opt.set_xtol_rel(1e-4)
x = opt.optimize([1.234, 5.678])
minf = opt.last_optimum_value()
print("optimum at ", x[0], x[1])
print("minimum value = ", minf)
print("result code = ", opt.last_optimize_result())
