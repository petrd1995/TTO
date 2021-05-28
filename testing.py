'''
    testing file
'''

import TTO
import numpy as np

example = TTO.Truss('example', 100, 100, 0, 2, 2, 2, np.array([[0,1,1,1],[1,1,1,1]]), np.array([[2,0,100, 0]]), 2.1e5, 10, 0, 0.2, 1, 1)

example.default_setup()
example.opt()


# example.Avec = np.ones_like(example.len) * example.A0

# example.MK()
# example.forces()
# K0 = example.K
# example.boundary()
# K1 = example.Ksys
# print(K0)
# print(K1)

# F = example.f
# Fsys = example.fsys

# u = Fsys
# u[:,1] = np.linalg.inv(K1) @ Fsys[:,1]
