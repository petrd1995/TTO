import TTO
import numpy as np
import emailnotify
'''
    Simple example of cube-like domain with 8 nodes in the corners - four of which are bound and one is loaded by single force. 
    Run with the following uncommented at the end of this file:
        example.default_setup()
        example.plot('bcs')
        example.opt()
        example.plot('res')

'''
x0 = 100 # length of domain in x direction [mm]
y0 = 100 # length of domain in y direction [mm]
z0 = 100 # length of domain in z direction [mm]
nx = 2  # number of nodes in x direction [-]
ny = 2 # number of nodes in y direction [-]
nz = 2  # number of nodes in z direction [-]
E = 2.1e5 # Young's modulus E[MPa]
r0 = 5 # initial guess for bar radii [mm]
ratio = 0.5 # ratio between final volume of material and initial V/V0 [-]
Vol0 = x0*y0*z0 # Initial volume from which we optimize to the desired ratio [mm^3]
Ro = 1 # Density of the material (would be used in dynamic analysis, which is not default)
kon = 1 # convergence of 1 % is set to end the optimization loop (when only one percentage is being made, the loop exits)
bcs = np.array([[0, 1, 1, 1], [1, 1, 1, 1],[2, 1, 1, 1], [3, 1, 1, 1]]) # 4 boundary conditions
f = np.array([[4, 0, -1000, -600]]) # one force

example = TTO.Truss('example', x0, y0, z0, nx, ny, nz, bcs, f, E, r0, Vol0, ratio, Ro, kon)


example.default_setup() # performs default setup of a problem - this means running:
    # example.create_grid() # creates ground structure's nodes
    # example.create_bars() # creates bars connecting these nodes
    # example.vec_len() # computes lengths and directions of these bars
example.plot('bcs') # plotting boundary conditions and forces
example.opt() # optimization
example.plot('res') # plotting results
