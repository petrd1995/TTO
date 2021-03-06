import numpy as np
import circlecoords
import itertools
import plots
import os
import pandas as pd
import create_bc_f as ct

"""
    This program contains class Truss, which contains several methods used to create and 
    solve basic truss topology optimization problems in 2D and 3D.
    Created by Petr David and Tomáš Mareš.
"""

class Truss:

    def __init__(self, name, x0, y0, z0, nx, ny, nz, bc, F, E, r0, Vol0, ratio, Ro, kon):

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

    def add_one_node(self, coords):

        # coords musi byt numpy array
        self.all_nodes = np.block([[self.all_nodes], [coords]])
        self.num_nodes = len(self.all_nodes)

    def add_circle(self, radius, plane, n_points, *origin):

        if self.z0:
            dim = 3
            x, y, z = circlecoords.circle(
                dim, plane, radius, n_points, *origin)
            newcoords = np.column_stack((x, np.column_stack((y, z))))
            self.add_one_node(newcoords)
        else:
            dim = 2
            x, y = circlecoords.circle(dim, plane, radius, n_points, *origin)
            newcoords = np.column_stack((x, y))
            self.add_one_node(newcoords)

    def plot(self, plot_type):

        if plot_type == 'grid':
            plots.plot_grid(self.all_nodes, self.bars)

        if plot_type == 'bcs':
            plots.plot_bcs(self.bc, self.F, self.all_nodes)

        if plot_type == 'res':
            plots.plot_res(self.res, self.vec, self.all_nodes, self.n)

        if plot_type == 'conv':
            plots.plot_conv(self.iteration, self.hist_epsilon)

    def node_coords(self, node_num, **kwargs):

        if 'show' in kwargs:
            print(f'the coordinates for node {node_num} are: {self.all_nodes[node_num]}')
        return self.all_nodes[node_num]

    def rem_node(self, node_num):

        self.all_nodes = np.delete(self.all_nodes, node_num, axis=0)
        self.num_nodes = len(self.all_nodes)

    def create_bars(self):

        self.num_bars = int(self.num_nodes * (self.num_nodes - 1) / 2)
        self.node_counter = np.arange(self.num_nodes)
        self.bars = np.empty((self.num_bars, 3), dtype=int)
        comb = itertools.combinations(self.node_counter, 2)

        for q, i in enumerate(comb):
            self.bars[q,]= int(q), *i

    def rem_bars(self, bar_num):

        self.bars = np.delete(self.bars, bar_num, axis=0)
        self.num_bars = len(self.bars)
        self.bars.T[0] = np.arange(self.num_bars)

    def rem_long_bars(self, lenm):
        '''remove bars of prescribed length multiple of length of diagonal bar'''

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
        '''remove bars of prescribed length multiple of length of diagonal bar'''

        a = []
        for i in range(self.num_bars):
            if (self.len[i] > length) or (self.len[i] < 0.5 * length):
                a.append(i)
        self.rem_bars(a)  # odstranění prutu/ů?
        self.vec = np.delete(self.vec, a, axis=0)
        self.len = np.delete(self.len, a, axis=0)

    def vec_len(self):

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

        self.K = np.zeros((int(self.cB * self.rB), int(self.cB * self.rB)))

        for i in range(self.num_bars):
            slc1 = self.cB * self.bars[i][1]
            slc1end = slc1 + self.cB
            slc2 = self.cB * self.bars[i][2]
            slc2end = slc2 + self.cB

            self.K[slc1:slc1end, slc1:slc1end] += np.outer(self.vec[i], self.vec[i]) * self.E * self.Avec[i] / self.len[i]

            self.K[slc2:slc2end, slc2:slc2end] += np.outer(self.vec[i], self.vec[i]) * self.E * self.Avec[i] / self.len[i]

            self.K[slc1:slc1end, slc2:slc2end] = - np.outer(self.vec[i], self.vec[i]) * self.E * self.Avec[i] / self.len[i]

            self.K[slc2:slc2end, slc1:slc1end] = - np.outer(self.vec[i], self.vec[i]) * self.E * self.Avec[i] / self.len[i]

    def forces(self):

        self.f = np.zeros((self.rB * self.cB, 1))
        for i in range(len(self.F)):
            slc = self.cB * int(self.F[i][0])
            self.f[slc] = self.F[i][1]
            self.f[slc + 1] = self.F[i][2]
            if self.z0:
                self.f[slc + 2] = self.F[i][3]

    def unit_forces(self):

        self.f = np.zeros((self.rB * self.cB, 3))

        for i in range(len(self.F)):
            slc = self.cB * int(self.F[i][0])
            self.f[slc,0] = 1
            self.f[slc,1] = 0
            self.f[slc,2] = 0
            self.f[slc + 1,0] = 0
            self.f[slc + 1,1] = 1
            self.f[slc + 1,2] = 0
            if self.z0:
                self.f[slc + 2,0] = 0
                self.f[slc + 2,1] = 0
                self.f[slc + 2,2] = 1

    def boundary(self):

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

        p1 = np.diag(self.K)
        p1 = np.array(p1)

        findzeros = np.where(p1 < 1e-5)
        for i in findzeros:
            p1[i] = 1
        np.fill_diagonal(self.K, 0)
        self.K = self.K + np.diag(p1)

    def opt(self):

        epsilon = 100
        maxit = 2000
        # rmax = 0.5 * np.min([self.x0 / self.nx, self.y0 / self.ny, self.z0 / self.nz])
        rmax = 15

        self.rB, self.cB = np.shape(self.all_nodes)
        self.Avec = np.ones_like(self.len) * self.A0
        self.unit_forces()

        self.Vol = self.ratio * self.Vol0

        self.iteration = 0

        self.hist_A = self.Avec
        '''Proměnná ukládající průřezy ve všech iteracích'''

        self.hist_epsilon = []
        '''Proměnná ukládající rozdíl mezi průřezy prutů dvou po sobě jdoucích iteracích'''


        # cf = []
        # cfcurrent = 0
        while epsilon > self.kon:
            self.iteration += 1
            # tvorba matice K
            self.matK()
            # zahrnutí OP
            self.boundary()
            # zohlednění nulových průřezů
            self.zerocrosssection()

            self.u = np.linalg.inv(self.K) @ self.f
            '''Vektor posuvů'''

            self.n = np.zeros((self.num_bars, 1))
            self.n1 = np.zeros((self.num_bars, 1))
            self.n2 = np.zeros((self.num_bars, 1))
            self.n3 = np.zeros((self.num_bars, 1))
            '''Vektor vnitřních sil'''

            # cfit = 0
            # počítání "lepšího" odhadu pomocí Lagrangeovy metody
            for i in range(self.num_bars):
                slc2 = self.cB * self.bars[i][2]
                slc1 = self.cB * self.bars[i][1]
                # self.n[i] = 0.25*(self.E * self.Avec[i] * float(self.vec[i] @ (self.u[slc2:slc2 + self.cB,0] - self.u[slc1:slc1 + self.cB,0])) / self.len[i]) + \
                #             0.25*(self.E * self.Avec[i] * float(self.vec[i] @ (self.u[slc2:slc2 + self.cB,1] - self.u[slc1:slc1 + self.cB,1])) / self.len[i]) + \
                #             0.5*(self.E * self.Avec[i] * float(self.vec[i] @ (self.u[slc2:slc2 + self.cB,2] - self.u[slc1:slc1 + self.cB,2])) / self.len[i])
                self.n1[i] = 0.25*(self.E * self.Avec[i] * float(self.vec[i] @ (self.u[slc2:slc2 + self.cB,0] - self.u[slc1:slc1 + self.cB,0])) / self.len[i])
                self.n2[i] = 0.25*(self.E * self.Avec[i] * float(self.vec[i] @ (self.u[slc2:slc2 + self.cB,1] - self.u[slc1:slc1 + self.cB,1])) / self.len[i])
                self.n3[i] = 0.5*(self.E * self.Avec[i] * float(self.vec[i] @ (self.u[slc2:slc2 + self.cB,2] - self.u[slc1:slc1 + self.cB,2])) / self.len[i])
            #     if self.Avec[i]:
            #         cfit += np.sqrt(self.n1[i]**2+self.n2[i]**2+self.n3[i]**2)**2*self.len[i]/(2*self.E*self.Avec[i])
            # cf.append(cfit)
            # print(cf[-1])
            Afrak1 = (self.n1 ** 2) / (2 * self.E)
            Afrak2 = (self.n2 ** 2) / (2 * self.E)
            Afrak3 = (self.n3 ** 2) / (2 * self.E)
            Afrak = 0.25*Afrak1 + 0.25*Afrak2 + 0.5*Afrak3
            Acurrent = (self.Vol * Afrak ** 0.5) / float(Afrak.T ** 0.5 @ self.len.reshape((self.num_bars, 1)))
            # odstranění prutů s příliš malým průřezem a zajištění maximálního možného průřezu
            # if self.iteration > 10:
            for i in range(len(Acurrent)):
                radius = np.sqrt(Acurrent[i]/np.pi)
                if radius < 1:
                    Acurrent[i] = 0.000000001
                elif radius > rmax:
                    Acurrent[i] = rmax**2*np.pi

            # epsilon = float(np.abs(cfcurrent-cf[-1])/cf[-1])*1000
            # cfcurrent = cf[-1]

            # computing epsilon (for convergence - difference between norms of two consecutive vectors of bar Areas)
            epsilon = np.linalg.norm(Acurrent - self.Avec.reshape((self.num_bars, 1)))
            # another way of computing epsilon
            # epsilon = np.linalg.norm(Acurrent - self.Avec.reshape((self.num_bars, 1)))/np.linalg.norm(self.Avec.reshape((self.num_bars, 1)))



            print(f"it: {self.iteration}, cfdiff = {epsilon}")
            self.Avec = Acurrent
            # ukončení smyčky v případě nekonvergence
            if self.iteration == maxit:
                epsilon = 0
                print("Maximum number of iterations reached")

            # ukládání průřezů z konkrétní iterace do do matice
            self.hist_A = np.column_stack((self.hist_A, self.Avec))
            # ukládání odchylky z dané iterace
            self.hist_epsilon.append(epsilon)

        print(np.dot(self.len.T,np.around(Acurrent, decimals=1))/self.Vol0, self.Vol)
        # ukládání výsledků
        self.res = np.column_stack((np.around(np.sqrt(Acurrent/np.pi), decimals=1), self.bars))
        # print(self.res)
        self.nonzero_res = np.empty((np.count_nonzero(self.res[:,0]),4))
        self.nonzero_res[:,] = self.res[np.nonzero(self.res[:,0]),:]
        # self.nonzero_len = np.empty_like(nonzero_res)


        self.u = np.linalg.inv(self.K) @ self.f[:,2]
        for i in range(18):
            print(i)
            print(self.u[(1344+i)*3:(1345+i)*3])

    def out(self):

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

        self.create_grid()
        self.create_bars()
        self.vec_len()
        # self.rem_long_bars_length(14)

def kostka_run():

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

    kostka = Truss('kostka2403', x0, y0, z0, nx, ny, nz, bcs, f, E, A0, ratio, Ro, kon)

    kostka.create_grid()

    a = (kostka.x0 + 20)
    b = (kostka.y0 + 20)
    c = (kostka.z0 + 20)

    kostka.add_circle(radius, 'z', bcnodes, a / 2, b / 2, 0)
    kostka.add_circle(radius, 'x', fnodes, 0, b / 2, c / 2)
    kostka.add_circle(radius, 'x', fnodes, a, b / 2, c / 2)
    kostka.add_circle(radius, 'y', fnodes, a / 2, 0, c / 2)
    kostka.add_circle(radius, 'y', fnodes, a / 2, b, c / 2)
    kostka.add_circle(radius, 'z', fnodes, a / 2, b / 2, c)
    kostka.after_nodes = np.shape(kostka.all_nodes)[0]

    kostka.create_bars()
    kostka.vec_len()
    kostka.rem_long_bars(1)

    kostka.opt()
    # kostka.out()
    # emailnotify.notify()
    kostka.plot('res')

    # kostka.plot('bcs')

def default_run(example):

    example.default_setup()
    example.opt()
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

# benchmark = Truss('benchmark', 100, 100, 100, 2, 2, 2, np.array([[0, 1, 1, 1], [1, 1, 1, 1],[2, 1, 1, 1], [3, 1, 1, 1]]), np.array([[4, 0, -1000, -600]]), 2.1e5, 100, 0.1, 100, 1)

# default_run(benchmark)

if __name__ == "__main__":
    kostka_run()
