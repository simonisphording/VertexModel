# -*- coding: utf-8 -*-

import pickle

seconds = 1000
stepsize = .1
nsteps = seconds/stepsize

### Initiation parameters ###

# Voronoi
radius = 6
center = [0., 0., 0.]
ncells = int(3.14 * radius ** 2)

### Force parameters ###

# Line tension
tension = .08
# Perimeter contractibility
contract = .1
# Surface tension
init_surface = 1
elastic = 4
# Gradient pull
gradient = 2

### Update parameters ###

growthrate = .01 * stepsize # increase in cell size with every step

# T1 parameters
transition_boundary = .2
new_edge_length = .3

# T3 parameters
paneth_dist = 3

def init_from_file(path):
    # Reads parameters from a .pickle file
    global nsteps, stepsize, radius, center, ncells, tension, contract
    global init_surface, elastic, gradient, growthrate, transition_boundary
    global new_edge_length, paneth_dist
    pkl = open(path, "rb")
    pars = pickle.load(pkl)
    pkl.close()
    seconds = pars['seconds']
    stepsize = pars['stepsize']
    nsteps = int(seconds/stepsize)
    radius = pars['radius']
    center = pars['center']
    ncells = int(3.14 * radius ** 2)
    tension = pars['tension']
    contract = pars['contract']
    init_surface = pars['init_surface']
    elastic = pars['elastic']
    gradient = pars['gradient']
    growthrate = .01 * stepsize
    transition_boundary = pars['transition_boundary']
    new_edge_length = pars['new_edge_length']
    paneth_dist = pars['paneth_dist']