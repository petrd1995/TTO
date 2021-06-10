import numpy as np
import itertools

from numpy.core.numeric import zeros_like
import plots
import os
import pandas as pd
import create_truss as ct

"""
    This program contains class Truss, which contains several methods used to create and 
    solve basic truss topology optimization problems in 2D and 3D.
    Created by Petr David and Tomáš Mareš.
"""

class Truss:

    def __init__(self, name, x0, y0, z0, nx, ny, nz, bc, F, E, r0, Vol0, ratio, Ro, kon):
        '''
            Initialization of class instance with variables defining the particular problem.
            The variables are as follows (for further info see README.md):
            
            name - name of the problem (for saving purposes)
            x0 - length of the design domain in x direction
            y0 - length of the design domain in y direction
            z0 - length of the design domain in z direction
            nx - number of nodes in x direction
            ny - number of nodes in y direction
            nz - number of nodes in z direction
            bc - vector(np.array) containing boundary conditions, which are explained by 
            F - vector(np.array) containing forces, similar to bcs:
            E - Young's modulus E of the material (MPa)
            r0 - initial radius of all bars (mm)
            Vol0 - particular volume from which optimal design is to be found
            ratio - ratio of volume of final design to initial Vol0, that is ratio = Vol/Vol0
            R0 - density of the material (kg/mm^3) for dynamic problems
            kon - konvergence criteria
        '''
        self.name = name
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.A0 = np.pi * r0**2
        self.E = E
        self.Ro = Ro
        self.F = np.array(F)
        self.bc = np.array(bc)
        self.ratio = ratio
        self.Ro = Ro
        self.kon = kon
        self.Vol0 = Vol0
        if not self.Vol0:
            if self.z0:
                self.Vol0 = self.x0 * self.y0 * self.z0
            else:
                self.Vol0 = self.x0 * self.y0 

    def create_grid(self):
        '''
            A sort of 'grid' creation from defined dimensions and number of nodes
            by means of numpy's linspace method.
            Variables here defined are:
                num_nodes - total number of nodes
                all_nodes -  contains info about individual nodes and their coordinates
                get_node - helper variable for creating nodes with unique coordinates
        '''
        self.x = np.linspace(0, self.x0, self.nx)
        self.y = np.linspace(0, self.y0, self.ny)
        self.z = np.linspace(0, self.z0, self.nz)

        if self.z0:
            self.num_nodes = self.nx * self.ny * self.nz
            self.all_nodes = np.empty((self.num_nodes, 3))
            self.get_node = itertools.product(self.x, self.y, self.z)

        else:
            self.num_nodes = self.nx * self.ny
            self.all_nodes = np.empty((self.num_nodes, 2))
            self.get_node = itertools.product(self.x, self.y)

        for i, el in enumerate(self.get_node):
            self.all_nodes[i] = el

    def grid_from_list(self, nodes):
        '''
            Method used, when one wishes to create ground structure from known set of nodes
            and wishes to evade creating linspace grid(from create_grid method) alltogether.
        '''
        self.num_nodes = nodes.shape[0]
        self.all_nodes = nodes

    def add_one_node(self, coords):
        '''
            Method to add node to the existing grid by passing the node's
            coordinates.
        '''
        # coords musi byt numpy array
        self.all_nodes = np.block([[self.all_nodes], [coords]])
        self.num_nodes = len(self.all_nodes)

    def plot(self, plot_type):
        '''
            Method for calling plot functions from the plots.py module.
            The plots are called as all other methods, that is if we have
            instance 'example', we can plot the initial grid, as well as
            optionally the connecting bars by calling:
            
                example.plot('grid')

            where we previously had to run the following methods:
                example.create_grid()
                example.create_bars()
                example.vec_len()
            
            to plot boundary conditions and forces
            (for e.g. verification purposes, that we created the correct problem)
            we call
                
                example.plot('bcs')

            to plot the resulting structure we call

                example.plot('res')
            
            where we previously had to run the following methods:
                example.create_grid()
                example.create_bars()
                example.vec_len()
                example.opt()

            to plot convergence we simply call:
            
                example.plot('conv')
            
            which has the same dependecies as plotting 'res'
        '''
        if plot_type == 'grid':
            plots.plot_grid(self.all_nodes, self.bars)

        if plot_type == 'bcs':
            plots.plot_bcs(self.bc, self.F, self.all_nodes)

        if plot_type == 'res':
            plots.plot_res(self.res, self.vec, self.all_nodes, self.n)

        if plot_type == 'conv':
            plots.plot_conv(self.iteration, self.hist_epsilon)

    def node_coords(self, node_num, **kwargs):
        '''
            Method used to return coordinates of wanted node by calling the node
            by its number
        '''
        if 'show' in kwargs:
            print(f'the coordinates for node {node_num} are: {self.all_nodes[node_num]}')
        return self.all_nodes[node_num]

    def rem_node(self, node_num):
        '''
            Method used to remove particular node(s) which is(are)
            deleted according to the inputed node number(s) - 
            their labels
        '''
        self.all_nodes = np.delete(self.all_nodes, node_num, axis=0)
        self.num_nodes = len(self.all_nodes)

    def create_bars(self):
        '''
            Method which creates bars from unique combinations of 
            individual nodes. This creates All unique combinations,
            as such overlapping bars are present(but can be removed later).
            Variables defined here are:
                num_bars - total number of bars
                node_counter - helper variable used to label individual bars
                bars - vector containing info about all bars (see README)
                comb - helper used to create unique combinations
        '''
        self.num_bars = int(self.num_nodes * (self.num_nodes - 1) / 2)
        self.node_counter = np.arange(self.num_nodes)
        self.bars = np.empty((self.num_bars, 3), dtype=int)
        comb = itertools.combinations(self.node_counter, 2)

        for q, i in enumerate(comb):
            self.bars[q,]= int(q), *i

    def rem_bars(self, bar_num):
        '''
            Method used to remove particular bars. Input is(are) bar(s) number(s).
        '''
        self.bars = np.delete(self.bars, bar_num, axis=0)
        self.num_bars = len(self.bars)
        self.bars.T[0] = np.arange(self.num_bars)

    def rem_long_bars(self, lenm):
        '''
            Remove bars of prescribed length multiple of length of diagonal bar
            of the smallest possible cell inside grid. For example if you set
            lenm = 1.1, then you remove all long bars and thus are left with
            a sort of cube-like lattice.
        '''

        a = []
        if self.z0:
            diag = np.sqrt((self.x0 / (self.nx - 1))**2 + (self.y0 / (self.ny - 1))**2 + (self.z0 / (self.nz - 1))**2)
        else:
            diag = np.sqrt((self.x0 / (self.nx - 1))**2 + (self.y0 / (self.ny - 1))**2)
        for i in range(self.num_bars):
            if (self.len[i] > lenm * diag) or (self.len[i] < 0.5 * lenm * diag):
                a.append(i)
        self.rem_bars(a)  # odstranění prutu/ů?
        self.vec = np.delete(self.vec, a, axis=0)
        self.len = np.delete(self.len, a, axis=0)

    def rem_long_bars_length(self, length):
        '''
            Similar to method rem_long_bars, but here the input is 
            a particular length, where longer bars than this are removed.
        '''

        a = []
        for i in range(self.num_bars):
            if (self.len[i] > length) or (self.len[i] < 0.5 * length):
                a.append(i)
        self.rem_bars(a)  # odstranění prutu/ů?
        self.vec = np.delete(self.vec, a, axis=0)
        self.len = np.delete(self.len, a, axis=0)

    def vec_len(self):
        '''
            Method for computing directional unit vectors and lengths
            of all bars in the ground structure(base truss).
            Variables defined here are:
                len - vector holding lengths of all bars
                vec - vector containing directional info about all bars,
                    every row corresponds to one bar's directional
                    unit vector's components x, y, (z). As with len
                    vector, the row number corresponds to the particular 
                    bar number (label)
        '''
        self.len = np.zeros((self.num_bars, 1))
        if self.z0:
            self.vec = np.zeros((self.num_bars, 3))
        else:
            self.vec = np.zeros((self.num_bars, 2))
        for bar in self.bars:

            start = self.node_coords(self.bars[bar[0], 1])
            end = self.node_coords(self.bars[bar[0], 2])

            self.len[bar[0]] = np.sqrt((end - start).dot(end - start))
            self.vec[bar[0], 0:2] = ((end - start) / self.len[bar[0]])[0:2]
            if self.z0:
                self.vec[bar[0], 2] = ((end - start) / self.len[bar[0]])[2]

    def matK(self):
        '''
            Method which creates the stiffness matrix K by performing
            the so called assembly operation (mapping of local bar's
            stiffness to the global frame of reference).
        '''
        self.gradK = np.zeros((self.num_bars, int(self.cB * self.rB), int(self.cB * self.rB)))
        self.lokK = zeros_like(self.gradK)
        
        for i in range(self.num_bars):
            slc1 = self.cB * self.bars[i][1]
            slc1end = slc1 + self.cB
            slc2 = self.cB * self.bars[i][2]
            slc2end = slc2 + self.cB

            self.gradK[i, slc1:slc1end, slc1:slc1end] += np.outer(self.vec[i], self.vec[i]) * self.E / self.len[i]

            self.gradK[i, slc2:slc2end, slc2:slc2end] += np.outer(self.vec[i], self.vec[i]) * self.E / self.len[i]

            self.gradK[i, slc1:slc1end, slc2:slc2end] = - np.outer(self.vec[i], self.vec[i]) * self.E / self.len[i]

            self.gradK[i, slc2:slc2end, slc1:slc1end] = - np.outer(self.vec[i], self.vec[i]) * self.E / self.len[i]

            
            self.lokK[i] = self.gradK[i] * self.Avec[i]
        self.K = np.sum(self.lokK, axis = 0)



    def forces(self):
        '''
            Function mapping forces and their components
            to corresponding degrees of freedom
        '''
        self.f = np.zeros((self.rB * self.cB, 1))
        for i in range(len(self.F)):
            slc = self.cB * int(self.F[i][0])
            self.f[slc] = self.F[i][1]
            self.f[slc + 1] = self.F[i][2]
            if self.z0:
                self.f[slc + 2] = self.F[i][3]

    def boundary(self):
        '''
            Function which assings boundary conditions by modifying
            the global stiffness matrix K. 
        '''
        for i in range(len(self.bc)):
            slc = self.cB * int(self.bc[i, 0])
            if self.bc[i, 1]:
                self.K[:, slc] = 0
                self.K[slc, :] = 0
                self.K[slc, slc] = 1
            if self.bc[i, 2]:
                self.K[:, slc + 1] = 0
                self.K[slc + 1, :] = 0
                self.K[slc + 1, slc + 1] = 1
            if self.z0:
                if self.bc[i, 3]:
                    self.K[:, slc + 2] = 0
                    self.K[slc + 2, :] = 0
                    self.K[slc + 2, slc + 2] = 1

    def zerocrosssection(self):
        '''
            Function which finds bars, whose cross-section
            got throughout optimization too small, and removes them
            entirely to prevent singularities. 
            (This may effect the final design!)
        '''
        p1 = np.diag(self.K)
        p1 = np.array(p1)

        findzeros = np.where(p1 < 1e-5)
        for i in findzeros:
            p1[i] = 1
        np.fill_diagonal(self.K, 0)
        self.K = self.K + np.diag(p1)

    def opt(self):

        force = self.forces()


        c_history = []
        history = []
        lokK, dKdA = maticeK()


        def cf(x, grad):

            for i in range(num_bars):
                lokK[i] = dKdA[i] * x[i]
            stiffness = np.sum(lokK, axis=0)
            Kres = boundary(stiffness)

            u = np.linalg.inv(Kres) @ force
            for i in range(num_bars):
                grad[i] = -u.T @ dKdA[i] @ u
            cost = u.T @ Kres @ u
            history.append(x)
            c_history.append(float(cost))
            # cost = u.T @ u
            # cost = force.T @ np.linalg.inv(Kres) @ force
            return float(cost)


        def myconstraint(x, grad):
            for i in range(num_bars):
                grad[i] = lengths[i]
            return float(x.T @ lengths - Vol0)


        lb = 0.0001*np.ones(num_bars)
        ub = 10000*np.ones(num_bars)
        ini_x = np.ones(num_bars)

        opt = nlopt.opt(nlopt.LD_MMA, num_bars)
        opt.set_lower_bounds(lb)
        opt.set_upper_bounds(ub)
        opt.set_min_objective(cf)
        opt.add_inequality_constraint(
            lambda x, grad: myconstraint(x, grad), 1e-8)
        opt.set_xtol_rel(1e-6)
        # x = opt.optimize(10*np.ones(num_bars))
        x = opt.optimize(ini_x)
        minf = opt.last_optimum_value()
        # print("optimum at ", x)
        print("minimum value = ", minf)
        # print("result code = ", opt.last_optimize_result())

        res = np.column_stack(
            (np.around(np.sqrt(x/np.pi), decimals=1), bars))

        n = np.zeros((num_bars, 1))

        for i in range(num_bars):
            lokK[i] = dKdA[i] * x[i]
        stiffness = np.sum(lokK, axis=0)
        Kres = boundary(stiffness)

        u = np.linalg.inv(Kres) @ force
        for i in range(num_bars):
            slc2 = cB * bars[i][2]
            slc1 = cB * bars[i][1]
            n[i] = E * x[i] * float(vec[i] @ (
                u[slc2:slc2 + cB] - u[slc1:slc1 + cB])) / lengths[i]

    def out(self):
        '''
            Function for saving the results into csv's. A folder of same name as the name
            inputed into instance initiation is created and into that folder the results are
            saved in csv format.
        '''
        def write(array, name, arrname):
            pd.DataFrame(array).to_csv(os.path.join(name, arrname + '.csv'), header=None, index=None)

        os.mkdir(self.name)
        write(self.nonzero_res, self.name, 'nonzero_res')
        write(self.res, self.name, 'res')
        write(self.all_nodes, self.name, 'all_nodes')
        write(self.vec, self.name, 'vec')
        write(self.len, self.name, 'lengths')

    def __repr__(self):

        return f'Truss({self.x0}, {self.y0}, {self.z0}, {self.nx}, {self.ny}, {self.nz})'

    def default_setup(self):
        '''
            Function which contains methods which are necessary to run every problem.
            It is just here for convenience purposes. (orderof the functions in which 
            they are called matters!)
        '''
        self.create_grid()
        self.create_bars()
        self.vec_len()
        # self.rem_long_bars_length(14)

def cube_run():
    '''
        Particular example along with a more 'difficult' setup. This can give you
        idea how the program is run..
        Easier example is shown at the end of the file.
    '''
    x0 = 240
    y0 = 240
    z0 = 240
    nx = 11
    ny = 11
    nz = 11
    E = 2.1e5
    A0 = 100
    ratio = 0.2
    Ro = 1
    kon = 1

    nnodes = nx * ny * nz
    fnodes = 10
    bcnodes = 8
    radius = x0 / 4 + 7
    after_nodes = nnodes + fnodes * 5 + bcnodes - 1
    f_array = np.array([0, 0, 100])
    bcnode_start = nnodes
    bcnode_end = nnodes + bcnodes - 1
    fnode_start = bcnode_end + 1
    bcs = ct.createBCs(bcnode_start, bcnode_end, 3)
    f = ct.createForces(fnode_start, after_nodes, 3, f_array)

    cube = Truss('cube2403', x0, y0, z0, nx, ny, nz, bcs, f, E, A0, ratio, Ro, kon)

    cube.create_grid()

    a = (cube.x0 + 20)
    b = (cube.y0 + 20)
    c = (cube.z0 + 20)

    cube.add_circle(radius, 'z', bcnodes, a / 2, b / 2, 0)
    cube.add_circle(radius, 'x', fnodes, 0, b / 2, c / 2)
    cube.add_circle(radius, 'x', fnodes, a, b / 2, c / 2)
    cube.add_circle(radius, 'y', fnodes, a / 2, 0, c / 2)
    cube.add_circle(radius, 'y', fnodes, a / 2, b, c / 2)
    cube.add_circle(radius, 'z', fnodes, a / 2, b / 2, c)
    cube.after_nodes = np.shape(cube.all_nodes)[0]

    cube.create_bars()
    cube.vec_len()
    cube.rem_long_bars(1)

    cube.opt()
    # cube.out()
    # emailnotify.notify()
    cube.plot('res')

    # cube.plot('bcs')

def default_run(example):
    '''
        Function for running individual examples. Uncomment as necessary(order 
        of the functions in which they are called matters!).
    '''
    example.default_setup()
    # example.opt()
    # example.plot('bcs')
    # example.plot('grid')
    # example.plot('res')
    # example.plot('conv')
    # example.out()


# b2D = Truss('b2D', 100, 100, 000, 3, 3, 3, np.array([[0, 1, 1], [1, 1, 1], [2, 1, 1]]), np.array([[6, 0, -1000]]), 2.1e5, 10, 0.2, 10000, 0.005)
# b3D = Truss('b3D', 100, 100, 100, 3, 3, 3, np.array([[0, 1, 1, 1], [1, 1, 1, 1], [2, 1, 1, 1],[3, 1, 1, 1], [4, 1, 1, 1], [5, 1, 1, 1]]), np.array([[12, 0, 0, -1000],[18, 0, 0, -1000]]), 2.1e5, 10, 0.1, 10000, 1)

# default_run(b3D)

# benchmark = Truss('convergence', 100, 100, 100, 10, 9, 2, np.array([[1, 1, 1, 1], [3, 1, 1, 1], [5, 1, 1, 1],[7, 1, 1, 1], [9, 1, 1, 1], [11, 1, 1, 1],
#     [13, 1, 1, 1], [15, 1, 1, 1], [17, 1, 1, 1]]), np.array([[171, 0, -1000, 0]]), 2.1e5, 10, 0.1, 10000, 1)

if __name__ == "__main__":
    '''
        I do not usually run this program itself, but create a different one, where i just
        call methods from this file, but here is an example anyways.
        This is a more easier example..
    '''
    benchmark = Truss('benchmark', 100, 100, 100, 2, 2, 2, np.array([[0, 1, 1, 1], [1, 1, 1, 1],[2, 1, 1, 1], [3, 1, 1, 1]]), np.array([[4, 0, -1000, -600]]), 2.1e5, 5, 0, 0.2, 100, 1)
    default_run(benchmark)
