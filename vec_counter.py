import numpy as np
import os
import csv
import json

path = 'C:\\Users\\david\\OneDrive - České vysoké učení technické v Praze\\Dizertace\\TTO\\TTO\\kostka0104_1_aktualni'

nonzero_res = np.genfromtxt(path + '\\nonzero_res.csv', delimiter=',')
vec = np.genfromtxt(path + '\\vec.csv', delimiter=',')
lengths = np.genfromtxt(path + '\\lengths.csv', delimiter=',')
all_nodes = np.genfromtxt(path + '\\all_nodes.csv', delimiter=',')
nonzero_vecs = vec[nonzero_res[:, 1].astype(int)]

vec_dict = {}
nnonzerobars = np.shape(nonzero_vecs)[0]

for i in range(0, nnonzerobars):
    counter = 0
    for j in range(i, nnonzerobars):
        if (nonzero_vecs[i, 0] == nonzero_vecs[j, 0]) and (nonzero_vecs[i, 1] == nonzero_vecs[j, 1]) and (nonzero_vecs[i, 2] == nonzero_vecs[j, 2]) and (f"{nonzero_vecs[i]}" not in vec_dict):
            counter += 1
    if counter != 0:
        vec_dict[f"{nonzero_vecs[i]}"] = counter

with open('vektory.csv', 'w') as f:
    for key in vec_dict.keys():
        f.write(f"{key}, {vec_dict[key]}\n")

with open("vektory.json", "w") as outfile:
    json.dump(vec_dict, outfile)
