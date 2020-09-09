# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 15:17:18 2019

@author: isphording
"""

import os
import pickle
import parameters as p
from matplotlib import pyplot as plt
from main import dist_paneth_center

pars = set()
for file in os.listdir():
    if 'sim_' in file:
        pars.add(file.split('_')[-1])

pars = [float(x) for x in pars]
pars.sort()
pars = [str(x) for x in pars]

full = [[],[],[]]
for par in pars:
    result = [[],[]]
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
            
            for i in range(len(p_cells)):
                result[0].append(dist_paneth_center(cells[p_cells[i]]))
                full[0].append(dist_paneth_center(cells[p_cells[i]]))
                full[2].append(par)
            
            path = file + '/data/' + str(end) + '.pickle'
            pkl = open(path, "rb")
            vertices, edges, cells, rand = pickle.load(pkl)
            pkl.close()
            
            for i in range(len(p_cells)):
                if p_cells[i] in cells:
                    result[1].append(dist_paneth_center(cells[p_cells[i]]))
                    full[1].append(dist_paneth_center(cells[p_cells[i]]))
                else:
                    result[1].append(p.radius)
                    full[1].append(p.radius)
    
    plt.xlabel("time (seconds)")
    plt.ylabel("distance to center")
    plt.title("Distance of each paneth cell over time")
    plt.scatter(result[0], result[1], c = 'black')
    plt.savefig(par + '_pdist_time.png')
    plt.close()

import matplotlib.colors

cmap = plt.cm.rainbow
norm = matplotlib.colors.Normalize(vmin = .8, vmax = 2.2)

from matplotlib.lines import Line2D
custom_lines = []

for par in pars:
    custom_lines.append(Line2D([0], [0], marker='o', color='w', markerfacecolor= cmap(norm(float(par))), markersize = 8))

fig, ax = plt.subplots()
ax.scatter(full[0], full[1], c = [cmap(norm(float(x))) for x in full[2]])
ax.legend(custom_lines, pars)
plt.xlabel('initial distance to center')
plt.ylabel('final distance to center')

plt.savefig('init_v_final.png')
plt.close()