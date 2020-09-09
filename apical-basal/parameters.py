# -*- coding: utf-8 -*-

import pickle

seconds = 1000
stepsize = .1
nsteps = int(seconds/stepsize)

### Force parameters ###

# Line tension
tension = .1
# Surface tension
basal_tension = .25
apical_tension = .25
lateral_tension = 0
# Volume
compression = 4
init_volume = 2.6

# Apical pressure
pressure = 0

### Update parameters ###

growthrate = .001 * stepsize

# T1 parameters
transition_boundary = .1 # Of both an edge and its sister edge summed
new_edge_length = .1

# T3 parameters
paneth_dist = 3

def init_from_file(path):
    global nsteps, stepsize, tension, basal_tension, apical_tension
    global lateral_tension, compression, init_volume, pressure, growthrate
    global growthrate, transition_boundary, new_edge_length, paneth_dist
    pkl = open(path, "rb")
    pars = pickle.load(pkl)
    pkl.close()
    seconds = pars['seconds']
    stepsize = pars['stepsize']
    nsteps = int(seconds/stepsize)
    tension = pars['tension']
    basal_tension = pars['basal_tension']
    apical_tension = pars['apical_tension']
    lateral_tension = pars['lateral_tension']
    compression = pars['compression']
    init_volume = pars['init_volume']
    pressure = pars['pressure']
    growthrate = pars['growthrate'] * stepsize
    transition_boundary = pars['transition_boundary']
    new_edge_length = pars['new_edge_length']
    paneth_dist = pars['paneth_dist']