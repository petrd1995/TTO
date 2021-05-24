import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import plotly.express as px
import plots

path = 'C:\\Users\\Petr\\OneDrive - České vysoké učení technické v Praze\\Dizertace\\TTO\\TTO\\f'


res = np.genfromtxt(path + '\\res.csv', delimiter=',')
vec = np.genfromtxt(path + '\\vec.csv', delimiter=',')
lengths = np.genfromtxt(path + '\\lengths.csv', delimiter=',')
all_nodes = np.genfromtxt(path + '\\all_nodes.csv', delimiter=',')

# a = plots.plot_res(res, vec, all_nodes, n)
bc = [[0, 1, 1, 1], [1, 1, 1, 1], [2, 1, 1, 1], [3, 1, 1, 1]]
f = [[5, 0, -1000, 0]]
