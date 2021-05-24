import numpy as np
import circlecoords
import itertools
import plots
import os
import pandas as pd
# import emailnotify

class Truss:

    def __init__(self, name, x0, y0, z0, nx, ny, nz, bc, F, E, r0, ratio, Ro, kon):

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

    def create_grid(self):

        self.x = np.linspace(20, self.x0, self.nx)
        self.y = np.linspace(20, self.y0, self.ny)
        self.z = np.linspace(20, self.z0, self.nz)

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
        self.num_nodes = np.shape(self.all_nodes)[0]

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
            plots.plot_grid(self.all_nodes, self.bars, plotbars=True)

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
        self.num_nodes = np.shape(self.all_nodes)[0]

    def create_bars(self):

        self.num_bars = int(self.num_nodes * (self.num_nodes - 1) / 2)
        self.node_counter = np.arange(self.num_nodes)
        self.bars = np.empty((self.num_bars, 3), dtype=int)
        comb = itertools.combinations(self.node_counter, 2)

        for q, i in enumerate(comb):
            self.bars[q, 0], self.bars[q, 1], self.bars[q, 2] = int(q), *i

    def rem_bars(self, bar_num):

        self.bars = np.delete(self.bars, bar_num, axis=0)
        self.num_bars = np.shape(self.bars)[0]
        self.bars.T[0] = np.arange(self.num_bars)

    def rem_long_bars(self, lenm):
        '''remove bars of prescribed length multiple of length of diagonal bar'''

        a = []
        if self.z0:
            diag = np.sqrt((self.x0 / (self.nx - 1))**2 + (self.y0 / (self.ny - 1))**2 + (self.z0 / (self.nz - 1))**2)
        else:
            diag = np.sqrt((self.x0 / (self.nx - 1))**2 + (self.y0 / (self.ny - 1))**2)
        for i in np.arange(self.num_bars):
            if (self.len[i] > lenm * diag) or (self.len[i] < 0.5 * lenm * diag):
                a.append(i)
        self.rem_bars(a)  # odstranění prutu/ů?
        self.vec = np.delete(self.vec, a, axis=0)
        self.len = np.delete(self.len, a, axis=0)

    def rem_long_bars_length(self, length):
        '''remove bars of prescribed length multiple of length of diagonal bar'''

        a = []
        for i in np.arange(self.num_bars):
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

        for i in np.arange(self.num_bars):
            self.K[self.cB * self.bars[i][1]:self.cB * self.bars[i][1] + self.cB, self.cB * self.bars[i][1]:self.cB * self.bars[i][1] + self.cB] =\
                self.K[self.cB * self.bars[i][1]:self.cB * self.bars[i][1] + self.cB,
                       self.cB * self.bars[i][1]:self.cB * self.bars[i][1] + self.cB] + np.outer(self.vec[i], self.vec[i]) * self.E * self.Avec[i] / self.len[i]

            self.K[self.cB * self.bars[i][2]:self.cB * self.bars[i][2] + self.cB, self.cB * self.bars[i][2]:self.cB * self.bars[i][2] + self.cB] =\
                self.K[self.cB * self.bars[i][2]:self.cB * self.bars[i][2] + self.cB,
                       self.cB * self.bars[i][2]:self.cB * self.bars[i][2] + self.cB] + np.outer(self.vec[i], self.vec[i]) * self.E * self.Avec[i] / self.len[i]

            self.K[self.cB * self.bars[i][1]:self.cB * self.bars[i][1] + self.cB, self.cB * self.bars[i][2]:self.cB * self.bars[i][2] + self.cB] = \
                - np.outer(self.vec[i], self.vec[i]) * self.E * self.Avec[i] / self.len[i]

            self.K[self.cB * self.bars[i][2]:self.cB * self.bars[i][2] + self.cB, self.cB * self.bars[i][1]:self.cB * self.bars[i][1] + self.cB] = \
                - np.outer(self.vec[i], self.vec[i]) * self.E * self.Avec[i] / self.len[i]

    def forces(self):

        self.f = np.zeros((self.rB * self.cB, 1))
        for i in np.arange(np.shape(self.F)[0]):
            self.f[self.cB * int(self.F[i][0])] = self.F[i][1]
            self.f[self.cB * int(self.F[i][0]) + 1] = self.F[i][2]
            if self.z0:
                self.f[self.cB * int(self.F[i][0]) + 2] = self.F[i][3]

    def boundary(self):

        for i in np.arange(np.shape(self.bc)[0]):
            if self.bc[i, 1]:
                self.K[:, self.cB * int(self.bc[i, 0])] = 0
                self.K[self.cB * int(self.bc[i, 0]), :] = 0
                self.K[self.cB * int(self.bc[i, 0]), self.cB * int(self.bc[i, 0])] = 1
            if self.bc[i, 2]:
                self.K[:, self.cB * int(self.bc[i, 0]) + 1] = 0
                self.K[self.cB * int(self.bc[i, 0]) + 1, :] = 0
                self.K[self.cB * int(self.bc[i, 0]) + 1, self.cB * int(self.bc[i, 0]) + 1] = 1
            if self.z0:
                if self.bc[i, 3]:
                    self.K[:, self.cB * int(self.bc[i, 0]) + 2] = 0
                    self.K[self.cB * int(self.bc[i, 0]) + 2, :] = 0
                    self.K[self.cB * int(self.bc[i, 0]) + 2, self.cB * int(self.bc[i, 0]) + 2] = 1

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
        Amax = 0.2 * np.min([self.x0 / self.nx, self.y0 / self.ny, self.z0 / self.nz])

        self.rB, self.cB = np.shape(self.all_nodes)
        self.Avec = np.ones_like(self.len) * self.A0
        self.forces()

        if self.z0:
            self.Vol0 = self.x0 * self.y0 * self.z0
        else:
            self.Vol0 = self.x0 * self.y0
        self.Vol = self.ratio * self.Vol0

        self.iteration = 0

        self.hist_A = self.Avec
        '''Proměnná ukládající průřezy ve všech iteracích'''

        self.hist_epsilon = []
        '''Proměnná ukládající rozdíl mezi průřezy prutů dvou po sobě jdoucích iteracích'''


        cf = []
        cfcurrent = 0
        while epsilon > self.kon:
            # print(f"it: {self.iteration}, eps: {np.around(epsilon, decimals=3)}")
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
            '''Vektor vnitřních sil'''

            cfit = 0
            # počítání "lepšího" odhadu pomocí Lagrangeovy metody
            for i in np.arange(self.num_bars):
                self.n[i] = self.E * self.Avec[i] * float(self.vec[i] @ (self.u[self.cB * self.bars[i][2]:self.cB * self.bars[i][2] + self.cB]
                                                                         - self.u[self.cB * self.bars[i][1]:self.cB * self.bars[i][1] + self.cB])) / self.len[i]
                if self.Avec[i]:
                    cfit += self.n[i]**2*self.len[i]/(2*self.E*self.Avec[i])
            cf.append(cfit)
            # print(cf[-1])

            Afrak = (self.n ** 2) / (2 * self.E)

            Acurrent = (self.Vol * Afrak ** 0.5) / float(Afrak.T ** 0.5 @ self.len.reshape((self.num_bars, 1)))

            epsilon = float(np.abs(cfcurrent-cf[-1])/cf[-1])*1000
            cfcurrent = cf[-1]
            

            # výpočet rozdílu průřezů po sobě jdoucích iterací
            # Acurr_sum = np.sum(Acurrent)
            # Aprev_sum = np.sum(self.Avec[:, -1])
            # epsilon_prev = epsilon
            # epsilon = np.abs(Acurr_sum - Aprev_sum) / Aprev_sum
            # print(f"it: {self.iteration}, eps: {np.around(epsilon*100, decimals=5)}, diffA = {np.abs(np.sum(Acurr_sum) - np.sum(Aprev_sum))}")
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
        # odstranění prutů s příliš malým průřezem
        for i in np.arange(np.shape(Acurrent)[0]):

            diam = np.sqrt(Acurrent[i]*4/np.pi)
            if diam < 2:
                Acurrent[i] = 0

        # ukládání výsledků
        self.res = np.column_stack((np.around(Acurrent, decimals=1), self.bars))
        # print(self.res)
        self.nonzero_res = np.empty((np.count_nonzero(self.res[:,0]),4))
        self.nonzero_res[:,] = self.res[np.nonzero(self.res[:,0]),:]

        print(self.nonzero_res)

    def max_bar_rad():
        unique_nodes = np.unique([self.nonzero_res[:,2],self.nonzero_res[:,3]])
        for i in np.arange(np.shape(unique_nodes)[0]):
            un_indices = np.isin(self.nonzero_res[:,2:4],unique_nodes[i]) # unique_nodes[i]
            # print(unique_nodes[i])
            nonzero_un_indices = np.nonzero(un_indices)
            max_in_node = np.max(self.nonzero_res[nonzero_un_indices[0],0])
            # print(max_in_node)
            print(np.sqrt(max_in_node*4/np.pi), self.node_coords(int(unique_nodes[i])))

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
        self.rem_long_bars_length(14)


def createBCs(node_start, node_end, dim):
    # přidá vetktnutí, dalo by se upravit na posuvnou vazbu změnou 1 na nulu u příslušného styčníku v příslušném směru
    if dim == 3:
        bcs = np.empty((node_end + 1 - node_start, 4))
        for en, i in enumerate(bcs):
            bcs[en][0], bcs[en][1], bcs[en][2], bcs[en][3] = en + node_start, 1, 1, 1
    else:
        bcs = np.empty((node_end + 1 - node_start, 3))
        for en, i in enumerate(bcs):
            bcs[en][0], bcs[en][1], bcs[en][2] = en + node_start, 1, 1
    return bcs

def createForces(node_start, node_end, dim, force):
    # přidá síly do daných uzlů, všechny síly jsou pro jednoduchost stejné
    if dim == 3:
        frcs = np.empty((node_end + 1 - node_start, 4))
        for en, i in enumerate(frcs):
            frcs[en][0], frcs[en][1], frcs[en][2], frcs[en][3] = en + node_start, force[0], force[1], force[2]
    else:
        frcs = np.empty((node_end + 1 - node_start, 3))
        for en, i in enumerate(frcs):
            frcs[en][0], frcs[en][1], frcs[en][2] = en + node_start, force[0], force[1]
    return frcs

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
    bcnodes = 7
    radius = x0 / 4 + 7
    after_nodes = nnodes + fnodes * 5 + bcnodes - 1
    f_array = np.array([0, 0, 100])
    bcnode_start = nnodes
    bcnode_end = nnodes + bcnodes - 1
    fnode_start = bcnode_end + 1
    bcs = createBCs(bcnode_start, bcnode_end, 3)
    f = createForces(fnode_start, after_nodes, 3, f_array)

    kostka = Truss('kostka2203_01', x0, y0, z0, nx, ny, nz, bcs, f, E, A0, ratio, Ro, kon)

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
    # kostka.plot('res')

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

benchmark = Truss('convergence', 100, 100, 100, 10, 9, 2, np.array([[1, 1, 1, 1], [3, 1, 1, 1], [5, 1, 1, 1],[7, 1, 1, 1], [9, 1, 1, 1], [11, 1, 1, 1],
    [13, 1, 1, 1], [15, 1, 1, 1], [17, 1, 1, 1]]), np.array([[171, 0, -1000, 0]]), 2.1e5, 10, 0.1, 10000, 1)

# default_run(benchmark)

kostka_run()
