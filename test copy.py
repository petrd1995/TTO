import numpy as np
import nlopt
import itertools
from scipy.optimize import minimize

x0 = 100
y0 = 100
z0 = 0
nx = 2
ny = 2
nz = 2
E = 2.1e5
r0 = 5
ratio = 0.5
Vol0 = x0*y0*ratio
Ro = 1
kon = 1
bcs = np.array([[0, 1, 1], [1, 1, 1]])
forces = np.array([[2, 0, -1000]])

all_nodes = np.array([[0, 0],
                        [0, 100],
                        [100, 0],
                        [100, 100]])

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

rB, cB = np.shape(all_nodes)



Avec = r0**2*np.pi*np.ones(num_bars)

def maticeK(A):
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

        lokK[i] = gradK[i] * A[i]
    stiffness = np.sum(lokK, axis=0)
    return stiffness, gradK


def force():
    f = np.zeros((rB * cB, 1))
    for i in range(len(forces)):
        slc = cB * int(forces[i][0])
        f[slc] = forces[i][1]
        f[slc + 1] = forces[i][2]
        if z0:
            f[slc + 2] = forces[i][3]
    return f

def boundary(mK, gK):
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
        for j in range(num_bars):
            if bcs[i, 1]:
                gK[j,:, slc] = 0
                gK[j,slc, :] = 0
                gK[j,slc, slc] = 1
            if bcs[i, 2]:
                gK[j,:, slc + 1] = 0
                gK[j,slc + 1, :] = 0
                gK[j,slc + 1, slc + 1] = 1
            if z0:
                if bcs[i, 3]:
                    gK[j,:, slc + 2] = 0
                    gK[j,slc + 2, :] = 0
                    gK[j,slc + 2, slc + 2] = 1
    return mK, gK

force = force()



# def cf(x, grad):

#     K, dKdA = maticeK(x)
#     Kres, dKdAres = boundary(K, dKdA)

#     u = np.linalg.inv(Kres) @ force
#     grad = np.array([None, None, None, None, None, None])
#     for i in range(num_bars):
#         grad[i] = -u.T @ dKdAres[i] @ u
#     cost = u.T @ Kres @ u
#     return cost


# opt = nlopt.opt(nlopt.LD_MMA, num_bars)
# opt.set_lower_bounds([0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
# opt.set_upper_bounds([100, 100, 100, 100, 100, 100])
# opt.set_min_objective(cf)
# opt.set_xtol_rel(1e-4)
# # x = opt.optimize(10*np.ones(num_bars))
# x = opt.optimize([10, 10, 10, 10, 10, 10])
# minf = opt.last_optimum_value()
# print("optimum at ", x)
# print("minimum value = ", minf)
# print("result code = ", opt.last_optimize_result())


# K, dKdA = maticeK(Avec)
# Kres, dKdAres = boundary(K, dKdA)
# u = np.linalg.inv(Kres) @ force

# print(float(u.T @ Kres @ u))

# print(Avec.T@lengths)


def cf(x):

    K, dKdA = maticeK(x)
    Kres, dKdAres = boundary(K, dKdA)
    u = np.linalg.inv(Kres) @ force
    
    return float(u.T @ Kres @ u)


cons = ({'type': 'ineq', 'fun': lambda x:  x.T@lengths - Vol0, 'jac':  })

bnds = ((1,100),
        (1, 100),
        (1, 100),
        (1, 100),
        (1, 100),
        (1, 100),)

initial_guess = 10*np.ones(num_bars)

res = minimize(cf, initial_guess, method='SLSQP', bounds=bnds,
               constraints=cons)

print(res)
