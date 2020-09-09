# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:42:21 2019

@author: isphording
"""

import os
import pickle
import parameters as p
import numpy as np
from math import sqrt
from statistics import stdev
from matplotlib import pyplot as plt
from main import dist_paneth_center

pars = set()
for file in os.listdir():
    if 'sim_' in file:
        pars.add(file.split('_')[-1])

result = [[],[],[]]
for par in pars:
    shortrange = []
    for file in os.listdir():
        if file.split('_')[-1] == par:
            p.init_from_file(file + '/parameters.pickle')
            
            end = 0
            for x in os.listdir(file + '/data'):
                if '.pickle' in x:
                    end = max(int(x[:-7]), end)
            end -= end%20
            
            p_cells = []
            path = file + '/data/0.pickle'
            pkl = open(path, "rb")
            vertices, edges, cells, rand = pickle.load(pkl)
            pkl.close()
            for c in cells:
                if cells[c].celltype == 'paneth':
                    p_cells.append(c)
            
            path = file + '/data/' + str(end) + '.pickle'
            pkl = open(path, "rb")
            vertices, edges, cells, rand = pickle.load(pkl)
            pkl.close()
            
            for i in range(len(p_cells)):
                if p_cells[i] in cells and dist_paneth_center(cells[p_cells[i]]) < 2:
                    shortrange.append(1)
                else:
                    shortrange.append(0)
    result[0].append(float(par))
    result[1].append(np.mean(shortrange))
    result[2].append(stdev(shortrange)/sqrt(len(shortrange)))

fig, ax = plt.subplots()
ax.errorbar(result[0], result[1], yerr = result[2], fmt = 'o')
#ax.scatter(result[0], result[1], c = 'black')
plt.xlabel('gradient force')
plt.ylabel('fraction of paneth cells')

plt.savefig('surfaces_areas.png')
plt.close()