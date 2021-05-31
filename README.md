# Truss topology optimization

Programs in this repo contain functions used to create and solve problems in truss topology optimization written in Python.

The core program is TTO.py, where the class Truss is defined along with its various methods.

*example of problem creation:*

    f = np.array([[4, 0, -1000, -600]])
    bcs = np.array([[0, 1, 1, 1], [1, 1, 1, 1],[2, 1, 1, 1], [3, 1, 1, 1]])

    benchmark = Truss('benchmark', 100, 100, 100, 2, 2, 2, bcs, f, 2.1e5, 5, 0, 0.2, 1, 1)

where the inputs are:

    Truss(name, x0, y0, z0, nx, ny, nz, bc, F, E, r0, Vol0, ratio, Ro, kon)

which means that would create problem with dimensions of 100x100x100 in x, y and z directions. Number of nodes would be 2 in x, y and z. Force with components [0, -1000, -600] would be put into node 4, and boundary conditions would be imposed onto nodes 0, 1, 2, 3. Young's modulus would be = 2.1e5 MPa. The initial guess for bar radii r0 = 5, Vol0 would be set to 0, which would simply mean the design domain itself would be Vol0, so that Vol0 = x0\*y0\*z0. Ratio = 0.2, density of material = 1 and covergence criteria = 1.

The individual inputs for creating the instance are further explained in [Variables used](#Variables_used) section and in example_run() or cube_run().

To run the problem we have to call some methods from the class:

    benchmark.create_grid()
    benchmark.create_bars()
    benchmark.vec_len()
    benchmark.opt()

## Quick use

To simply run some examples use files cube_run.py and example_run.py, which are set up and should run out of the box. These files provide a look into how a problem is defined and which functions in what order are called (not all are used).
### My general workflow
I usually create a separate file from which I call individual methods, as is the case in the two example files cube_run.py and example_run.py.

First I create instance of the class Truss, which stands as the problem itself. Initially I set the force vector and boundary conditions vector as empty or e.g.:

    bc = [[0,1,1]]
    f = [[1,0,0]]

where I just want to make sure the instance is created, despite having wrong bcs and forces. Then I proceed with calling individual methods to set up nodes, bars etc. with following functions:

    name_of_instance.create_grid()
    name_of_instance.create_bars()
    name_of_instance.vec_len()

after these, I can plot the grid with individual node labels, in order to find out in which nodes I want to put bcs/forces.

    name_of_instance.plot('bcs')

where one just has to make sure, that in plots.py the corresponding function plot_grid() has uncommented the blocks of code responsible for showing labels. If one uses large number of nodes (let's say more than 10 in each direction) the preferred way would be to simply add nodes by using the method add_one_node() as is shown in cube_run.py:

    name_of_instance.create_grid()
    name_of_instance.add_one_node(nds_to_add)

note that running create_grid() before this method is necessary. The array nds_to_add should contain all desired nodes to be added, in the following way:

    nds_to_add = np.array[[x1,y1,z1], 
                          [x2,y2,z2]]

where one adds nodes based on their coordinates (depending on problem dimension also z coordinate). Here are shown two nodes to be added. These are appended to the rest of the nodes and so to access the first added node (x1, y1, z1) its number(label) will be {total number of nodes + 1} but since Python counts from 0 it will be just equal to total number of nodes. Let's say I created a problem which has 8 nodes and added one node. In order to put force (e.g. of 100 N in x and 200 N in y direction) into this node I would create array f like so:

    f = [[8, 100, 200]]

This is thus done before creating the instance, but after doing so the add_one_node() part has to follow after the create_grid() method otherwise we would get error, that we are trying to access node number 8(9th node) but we only have 8 nodes.

To check if bcs and forces were created as expected we can call:

    name_of_instance.plot('bcs')

If everything is set up the way we want we can finally call:

    name_of_instance.opt()

after which we can do:

    name_of_instance.plot('res')

to plot the results, or:

    name_of_instance.out()

to save the results into csv files in a direction of the name of the instance (the first argument to the instance creation).

## Files and their description:

TTO.py - the main file containing the class Truss, which has several methods used to create and 
    solve basic truss topology optimization problems in 2D and 3D. 

TTO_cube.py - modified version of TTO.py which is more suited to solve problems with multiple load cases (experimental).

plots.py - file containing functions for plotting grid, boundary conditions, results etc. Was separated from the main file for simplicity of the class Truss.

example_run.py - file with simple problem setup to get an idea how to setup problems yourself.

emailnotify.py - function for sending e-mails programmatically. I use this with complex problems, where there is a lot of computation time to notify me that the optimization has ended. Works only for gmail and one has to perform setup as shown in the file itself.

dynamics.py - file containing functions for solving modal analysis(eigenproblems). Not included in optimization yet but TBD...

cube_run.py - another example problem. This one is more complex.

create_bc_f.py - file with two functions used to create arrays of bcs and forces programmatically - for larger problems with a lot of bcs and forces.

circlecoords.py - older file used to create nodes on a circle. Was used for specific boundary condition shapes.

testing.py - just a testing file

## <a name="Variables_used"></a>Variables used:

**name** - name of the problem (for saving purposes)

**x0** - length of the design domain in x direction (mm)

**y0** - length of the design domain in y direction (mm)

**z0** - length of the design domain in z direction (mm)

**nx** - number of nodes in x direction

**ny** - number of nodes in y direction

**nz** - number of nodes in z direction

**bc** - vector(np.array) containing boundary conditions, which are explained by example:

    bc = [[0,1,1],[3,1,0]]

    [node_nr, bcx, bcy, (bcz)]

where two vectors mean we have two bcs (or bc in two nodes) the nodes in which the bcs are present are defined by the firs column of each vector, so the first bc will be set in node number 0 and the second in node number 3. The remaining numbers define if corresponding degree of freedom is constrained or not - 1 means constrained, 0 means free. Node 0 is therefore constrained in both x and y, whereas node 3 is only constrained in x, and free in y. For 3D problems it is needed to add one more column for z direction.

**F** - vector(np.array) containing forces (N), similar to bcs:

    F = [[2, 0, -1000, 0],[4, 10, -1000, 20]]
    
    [node_nr, Fx, Fy, (Fz)]

where two vectors mean two forces are defined. First number in each vector (row) again means node number in which the force is present, and the other numbers mean components of the individual forces.

**E** - Young's modulus E of the material (MPa)

**r0** - initial radius of all bars (mm)

**Vol0** - particular volume from which optimal design is to be found in further defined ratio, if Vol0 is to be defined by the volume of the design domain leave Vol0 = 0

**ratio** - ratio of volume of final design to initial Vol0, that is ratio = Vol/Vol0

**R0** - density of the material (kg/mm^3) for dynamic problems

**kon** - convergence criteria, which is later defined as difference of two consecutive
    design of areas of individual bars. kon = 1 means that when only improvement
    of one percent is achieved, the optimization loop ends

**num_nodes** - total number of nodes 

**all_nodes** -  contains info about individual nodes and their coordinates, that means it
    is a vector, where every row corresponds to coordinates of particular node

**get_node** - helper variable for creating nodes with unique coordinates

**num_bars** - total number of bars

**node_counter** - helper variable used to label individual bars

**bars** - vector containing info about all bars. Each row corresponds to one bar, where
    the first column contains the bar's number (label), second column is the number (label)
    of starting node, the third column the number (label) of ending node

**comb** - helper used to create unique combinations

**epsilon** - initial difference between two consecutive designs
    (must be initially set as bigger than variable self.kon
    afterwards get recomputed automatically)

**maxit** - maximum number of iterations 

**rmax** - maximum radius of bars (further info README)

**rB** - number of rows of array all_nodes ~ final number of nodes

**cB** - number of columns of array all_nodes ~ problem dimension (2D, 3D)

**iteration** - the current iteration number

**hist_A** - array containing areas of all iterations

**hist_epsilon** - array containing all previous epsilons

**u** - array of displacements

**n** - array of bar's inner forces

**res** - array containing results, first column contains bar areas, and further columns are as in variable bars

**nonzero_res** - same as res, but only for non-zero bars