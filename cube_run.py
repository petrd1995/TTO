import TTO_cube
import numpy as np
from create_bc_f import createBCs, createForces
import emailnotify
'''
    A problem of finding optimal topology based multiple load cases in particular setup.
    To run the following must be present at the end of the file(after nds_to_add definition):

        cube.create_grid() # creates grid
        cube.add_one_node(nds_to_add) # adds set of nodes in particular coordinates, which will be used as nodes for boundary conditions and forces
        cube.create_bars() # creates bars
        cube.vec_len() # computes directions and lengths of bars
        cube.rem_long_bars(1) # so that the optimization doesn't take too long
            # cube.plot('grid') # insert this line if you want to see the grid from which we will optimize
            # cube.plot('bcs') # insert this if curious about boundary conditions
        cube.opt() # optimization
        cube.plot('res') # plotting resulting structure
'''
x0 = 100  # length of domain in x direction [mm]
y0 = 100  # length of domain in y direction [mm]
z0 = 100  # length of domain in z direction [mm]
nx = 2  # number of nodes in x direction [-]
ny = 2  # number of nodes in y direction [-]
nz = 2  # number of nodes in z direction [-]
E = 2.1e5  # Young's modulus E[MPa]
r0 = 5  # initial guess for bar radii [mm]
ratio = 0.5  # ratio between final volume of material and initial V/V0 [-]
# Initial volume from which we optimize to the desired ratio [mm^3]
Vol0 = x0*y0*z0
# Density of the material (would be used in dynamic analysis, which is not default)
Ro = 1
# convergence of 1 % is set to end the optimization loop (when only one percentage is being made, the loop exits)

x0 = 280 # length of domain in x direction [mm]
y0 = 280 # length of domain in y direction [mm]
z0 = 280 # length of domain in z direction [mm]
nx = 11 # number of nodes in x direction [-]
ny = 11 # number of nodes in y direction [-]
nz = 11 # number of nodes in z direction [-]
E = 2.1e5  # Young's modulus E[MPa]
# E = 1
A0 = 10  # initial guess for bar Area [mm^2]
ratio = 0.4  # ratio between final volume of material and initial V/V0 [-]
Ro = 1 # Density of the material (would be used in dynamic analysis, which is not default)

kon = 10 # convergence of 1 % is set to end the optimization loop (when only one percentage is being made, the loop exits)
Vol0 = 20640708 # Initial volume [mm^3]

bcnodes = 13 # number of nodes to be added for boundary conditions
fnodes = 18 # number of nodes to be added for forces
f_array = np.array([0, 0, 1]) # aribtrary vector for forces - i later changed it in the TTO_cube to compute multiple load cases in all directions by itself, so this step is obsolete

nnodes = nx * ny * nz # total number of nodes
bcnode_start = nnodes # boundary condition's start node
bcnode_end = nnodes + bcnodes - 1  # boundary condition's end node
fnode_start = bcnode_end + 1 # starting node for forces
after_nodes = nnodes + fnodes + bcnodes - 1  # ending node for forces

bcs = createBCs(bcnode_start, bcnode_end, 3) # creating boundary conditions
f = createForces(fnode_start, after_nodes, 3, f_array) # creating forces

cube = TTO_cube.Truss('cube1', x0, y0, z0, nx, ny, nz, bcs, f, E, A0, Vol0, ratio, Ro, kon)

# adding nodes for boundary conditions and forces
# the first 13 are for BCs
# the last 18 are for forces
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

cube.create_grid()  # creates grid

cube.add_one_node(nds_to_add) # adds set of nodes in particular coordinates, which will be used as nodes for boundary conditions and forces
cube.create_bars()  # creates bars
cube.vec_len()  # computes directions and lengths of bars
cube.rem_long_bars(1)  # so that the optimization doesn't take too long
# cube.plot('grid') # insert this line if you want to see the grid from which we will optimize
# cube.plot('bcs') # insert this if curious about boundary conditions
cube.opt()  # optimization
cube.plot('res')  # plotting resulting structure
