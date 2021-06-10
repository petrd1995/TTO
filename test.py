import numpy as np
import nlopt
import itertools
import plots

x0 = 100
y0 = 100
z0 = 0
nx = 6
ny = 5
nz = 2
E = 2.1e5
r0 = 5
ratio = 0.5
Vol0 = 10000*ratio
Ro = 1
kon = 1
bcs = np.array([[0, 1, 1], [1, 1, 1], [2, 1, 1],
               [3, 1, 1], [4, 1, 1]])
forces = np.array([[27, 0, -1000]])
# bcs = np.array([[0, 1, 1], [1, 1, 1], [2, 1, 1],
#                [3, 1, 1], [4, 1, 1], [5, 1, 1], [6, 1, 1], [7, 1, 1], [8, 1, 1]])
# forces = np.array([[112, 0, -1000]])

x = np.linspace(0, x0, nx)
y = np.linspace(0, y0, ny)
z = np.linspace(0, z0, nz)

if z0:
    num_nodes = nx * ny * nz
    all_nodes = np.empty((num_nodes, 3))
    get_node = itertools.product(x, y, z)

else:
    num_nodes = nx * ny
    all_nodes = np.empty((num_nodes, 2))
    get_node = itertools.product(x, y)

for i, el in enumerate(get_node):
    all_nodes[i] = el




# all_nodes = np.array([[0, 0],
#                       [0, 100],
#                       [100, 0],
#                       [100, 100]])

num_nodes = len(all_nodes)
num_bars = int(num_nodes * (num_nodes - 1) / 2)
node_counter = np.arange(num_nodes)
bars = np.empty((num_bars, 3), dtype=int)
comb = itertools.combinations(node_counter, 2)
for q, i in enumerate(comb):
    bars[q, ] = int(q), *i


def node_coords(node_num, all_nodes):
    return all_nodes[node_num]


lengths = np.zeros((num_bars, 1))
vec = np.zeros((num_bars, 2))
for bar in bars:
    start = node_coords(bars[bar[0], 1], all_nodes)
    end = node_coords(bars[bar[0], 2], all_nodes)

    lengths[bar[0]] = np.sqrt((end - start).dot(end - start))
    vec[bar[0], 0:2] = ((end - start) / lengths[bar[0]])[0:2]
    if z0:
        vec[bar[0], 2] = ((end - start) / lengths[bar[0]])[2]

def rem_bars(bar_num, barx):

    barx = np.delete(barx, bar_num, axis=0)
    num_barx = len(barx)
    barx.T[0] = np.arange(num_barx)

    return barx

def rem_long_bars(lenm, vecs, lens, bar_arr):
    a = []
    if z0:
        diag = np.sqrt((x0 / (nx - 1))**2 + (y0 /
                        (ny - 1))**2 + (z0 / (nz - 1))**2)
    else:
        diag = np.sqrt((x0 / (nx - 1))**2 +
                        (y0 / (ny - 1))**2)
    for i in range(num_bars):
        if (lens[i] > lenm * diag) or (lens[i] < 0.5 * lenm * diag):
            a.append(i)
    bar_arr = rem_bars(a, bar_arr)  # odstranění prutu/ů?
    vecs = np.delete(vecs, a, axis=0)
    lens = np.delete(lens, a, axis=0)

    return vecs, lens, bar_arr


print(num_bars)
vec, lengths, bars = rem_long_bars(1.1, vec, lengths, bars)
print(num_bars)
rB, cB = np.shape(all_nodes)


Avec = r0**2*np.pi*np.ones(num_bars)

def maticeK():
    gradK = np.zeros((num_bars, int(cB * rB), int(cB * rB)))

    lokK = np.zeros_like(gradK)
    for i in range(num_bars):
        slc1 = cB * bars[i][1]
        slc1end = slc1 + cB
        slc2 = cB * bars[i][2]
        slc2end = slc2 + cB

        gradK[i, slc1:slc1end, slc1:slc1end] += np.outer(
            vec[i], vec[i]) * E / lengths[i]

        gradK[i, slc2:slc2end, slc2:slc2end] += np.outer(
            vec[i], vec[i]) * E / lengths[i]

        gradK[i, slc1:slc1end, slc2:slc2end] = - \
            np.outer(vec[i], vec[i]) * E / lengths[i]

        gradK[i, slc2:slc2end, slc1:slc1end] = - \
            np.outer(vec[i], vec[i]) * E / lengths[i]

    return lokK, gradK


def force():
    f = np.zeros((rB * cB, 1))
    for i in range(len(forces)):
        slc = cB * int(forces[i][0])
        f[slc] = forces[i][1]
        f[slc + 1] = forces[i][2]
        if z0:
            f[slc + 2] = forces[i][3]
    return f


def boundary(mK):
    for i in range(len(bcs)):
        slc = cB * int(bcs[i, 0])
        if bcs[i, 1]:
            mK[:, slc] = 0
            mK[slc, :] = 0
            mK[slc, slc] = 1
        if bcs[i, 2]:
            mK[:, slc + 1] = 0
            mK[slc + 1, :] = 0
            mK[slc + 1, slc + 1] = 1
        if z0:
            if bcs[i, 3]:
                mK[:, slc + 2] = 0
                mK[slc + 2, :] = 0
                mK[slc + 2, slc + 2] = 1
    return mK


force = force()
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

# plots.plot_grid(all_nodes, bars)
# plots.plot_bcs(bcs, forces, all_nodes)
plots.plot_res(res, vec, all_nodes, n)
# plots.plot_conv(iteration, hist_epsilon)
