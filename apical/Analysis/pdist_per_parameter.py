# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 12:01:31 2019

@author: isphording
"""

import os
import pickle
import parameters as p
import numpy as np
from matplotlib import pyplot as plt
from main import dist_paneth_center



pars = set()
for file in os.listdir():
    if 'sim_' in file:
        pars.add(file.split('_')[-1])
pars = list(pars)

for par in pars:
    result = []
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
            
            t = []
            pdist = [[] for n in p_cells]
            
            for seq in range(0, end + 20, 20):
                path = file + '/data/' + str(seq) + '.pickle'
                pkl = open(path, "rb")
                vertices, edges, cells, rand = pickle.load(pkl)
                pkl.close()
                for i in range(len(p_cells)):
                    if p_cells[i] in cells:
                        pdist[i].append(dist_paneth_center(cells[p_cells[i]]))
                    else:
                        pdist[i].append(None)
                t.append(seq * p.stepsize)
            result.append(pdist)
    
    
    
    plt.xlabel("time (seconds)")
    plt.ylabel("distance to center")
    plt.title("Distance of each paneth cell over time")
    colors = ["crimson", "limegreen", "fuchsia", "navy", "darkorange"]
    for x in range(len(result)):
        pdist = result[x]
        for i in range(len(pdist)):
            plt.plot(t, [y for y in pdist[i]], label = 'paneth cell %d'%i, color = colors[x])
    plt.savefig(par + '_pdist_time.png')
    plt.close()