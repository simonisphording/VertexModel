# -*- coding: utf-8 -*-
"""
Created on Tue May 28 12:25:52 2019

@author: isphording
"""

import os
import pickle
import time
import matplotlib as mpl
from matplotlib import pyplot as plt
from main import n_stem_cells, dist_paneth_center,\
avg_radial_distribution, gif_from_data
import parameters as p

mpl.rcParams.update({'figure.max_open_warning': 0})

# For every folder:
#   Make a video
#   Plot #Stemcells over time
#   Plot

end = 0
for x in os.listdir('data'):
    if '.pickle' in x:
        end = max(int(x[:-7]), end)
end -= end%20

p.init_from_file('parameters.pickle')

## Timelapse
gif_from_data('data')

## Radial distribution
file = 'data/' + str(end) + '.pickle'
pkl = open(file, "rb")
vertices, edges, cells, rand = pickle.load(pkl)
pkl.close()

data = avg_radial_distribution(cells, .1)
plt.plot(range(len(data)),data)
plt.savefig(os.getcwd().split(os.sep)[-1] + '_dist_correlation.png')
plt.close()

## Plots over time

p_cells = []
file = 'data/0.pickle'
pkl = open(file, "rb")
vertices, edges, cells, rand = pickle.load(pkl)
pkl.close()
for c in cells:
    if cells[c].celltype == 'paneth':
        p_cells.append(c)

nstem = []
t = []
pdist = [[] for n in p_cells]

for seq in range(0, end + 20, 20):
    file = 'data/' + str(seq) + '.pickle'
    pkl = open(file, "rb")
    vertices, edges, cells, rand = pickle.load(pkl)
    pkl.close()
    nstem.append(n_stem_cells(cells))
    for i in range(len(p_cells)):
        if p_cells[i] in cells:
            pdist[i].append(dist_paneth_center(cells[p_cells[i]]))
        else:
            pdist[i].append(None)
    t.append(seq * p.stepsize)

# N stem cells over time
plt.xlabel("time (seconds)")
plt.ylabel("N Stem cells")
plt.title("Number of stem cells over time")
plt.plot(t, nstem)
plt.savefig(os.getcwd().split(os.sep)[-1] + '_nstem_time.png')
plt.close()

# Distance to distance for each paneth cell over time
plt.xlabel("time (seconds)")
plt.ylabel("distance to center")
plt.title("Distance of each paneth cell over time")
for i in range(len(pdist)):
    plt.plot(t, [y for y in pdist[i]], label = 'paneth cell %d'%i)
plt.savefig(os.getcwd().split(os.sep)[-1] + '_pdist_time.png')
plt.close

time.sleep(0.1)
