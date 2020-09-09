# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 14:29:39 2019

@author: isphording
"""

import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from vertex_objects import Vertex, Cell, Edge
from voronoi import initiate_cells
from transitions import t1, t3
import parameters as p
import math
from numpy.random import random, get_state
import pickle
import imageio
import sys
import os
sys.setrecursionlimit(100000)

### Model functions ###
    
def system_energy(vertices, edges, cells):
    # Calculating the free energy of the system
    surface = 0
    perimeter = 0
    line = 0
    curvature = 0
    for c in cells:
        cell = cells[c]
        surface += p.elastic / 2 * (norm(cell.surface()) - cell.target_volume)**2
        perimeter += p.contract / 2 * cell.perimeter()**2
        if cell.celltype == 'paneth':
            curvature += p.gradient * norm(cell.center() - p.center)
    for edge in edges:
        line += p.tension * edges[edge].length()
    return surface + perimeter + line + curvature

def update(vertices, edges, cells, step = .1):
    # Updates the system to the next time step
    force = {}
    for v in vertices:
        if not vertices[v].bound:
            force[v] = vertices[v].force()
        else:
            force[v] = np.array([0])
    for v in vertices:
        vertices[v].location += force[v] * step
    for e in edges:
        if edges[e].length() < p.transition_boundary:
            t1(edges[e], vertices, edges)
    # Division of stem cells
    for cell in list(cells):
        if cells[cell].celltype == 'stem':
            cells[cell].target_volume += cells[cell].target_volume * p.growthrate
            if cells[cell].surface() > 2 * p.init_surface:
                new_cells, vertices, edges, cells = t3(cells[cell], vertices, edges, cells)
                # Differentiation of newly formed stem cells
                for c in new_cells:
                    if c.dist_to_paneth(cells) > p.paneth_dist:
                        c.celltype = None
                        c.random_color()
    # Delete cells with centers outside of the perimeter
    for c in list(cells):
        if norm(cells[c].center() - p.center) > p.radius:
            vertices, edges, cells = remove_outer_cell(cells[c], vertices, edges, cells)
        elif cells[c].celltype == 'paneth' and any([len(e.cells) < 2 for e in cells[c].edges]):
            vertices, edges, cells = remove_outer_cell(cells[c], vertices, edges, cells)
    
    return(vertices, edges, cells)

def remove_outer_cell(cell, vertices, edges, cells):
    # Removes one outer cell and its dependencies
    for vertex in cell.vertices():
        vertex.remove_cell(cell)
        if not vertex.cells:
            del vertices[vertex.name]
    for edge in list(cell.edges):
        edge.remove_cell(cell)
        if not edge.cells:
            del edges[edge.name]
    del cells[cell.name]
    return vertices, edges, cells

### Parsing and plotting ###

def parse_setup(file, silent = False):
    # Initiates the model from a .txt file
    vertices = {}; cells = {}; edges = {}
    file = open(file, 'r')
    for line in file:
        line = line.split()
        if not line:
            if silent == False: 
                print('line is empty')
        elif line[0] == 'vertex':
            vertices[line[1]] = Vertex([float(line[2]), float(line[3]), float(line[4])], name = line[1])
            vertices[line[1]].bound = (line[5] == 'bound')
        elif line[0] == 'edge':
            edges[line[1]] = Edge([vertices[line[2]], vertices[line[3]]], name = line[1])
            vertices[line[2]].add_edge(edges[line[1]])
            vertices[line[3]].add_edge(edges[line[1]])
        elif line[0] == 'cell':
            cells[line[1]] = Cell(name = line[1])
            cell = cells[line[1]]
            cell.random_color()
            e = []
            d = []
            for i in range(2, len(line)):
                if line[i][0] == 'e':
                    e.append(edges[line[i]])
                    edges[line[i]].add_cell(cell)
                else:
                    d.append(int(line[i]))
            cell.direction = d
            cell.edges = e
            # Adding cell to vertices
            for edge in e:
                for v in edge.vertices:
                    v.add_cell(cell)
        else:
             if silent == False: 
                print('line does not start with vertex, edge or cell')
    return (vertices, edges, cells)

def model_to_txt(file, vertices, edges, cells):
    # Writes the current system to a .txt file
    file = open(file, 'w')
    for v in vertices:
        file.write('vertex ' + v)
        for dim in vertices[v].location:
            file.write(' ' + str(dim))
        if vertices[v].bound:
            file.write(' bound\n')
        else:
            file.write(' unbound\n')
    for e in edges:
        file.write('edge ' + e)
        for v in edges[e].vertices:
            file.write(' ' + v.name)
        file.write('\n')
    for c in cells:
        file.write('cell ' + c)
        for e in cells[c].edges:
            file.write(' ' + e.name + ' ' + str(cells[c].edge_direction(e)))
        file.write('\n')
    file.close()

def plot_space(verts, es, cs, elev = 10, azim = 120, radius = None, center = None,\
               vnumbers = False, enumbers = False, cnumbers = False, video = False):
    # Plotting the system
    if video:
        fig = plt.figure(figsize = (20,20))
    else:
        fig = plt.figure(figsize = (10,10))
    ax = Axes3D(fig)

    xs = [verts[v].location[0] for v in verts]
    ys = [verts[v].location[1] for v in verts]
    zs = [verts[v].location[2] for v in verts]
    
    patches = []; colors = []
    for c in cs:
        locs = [vertex.location for vertex in cs[c].ordered_vertices()]
        patches.append(locs)
        colors.append(cs[c].color)
    pc = Poly3DCollection(patches, facecolors=colors)
    #pc.set_facecolor('c')
    #pc.set_array(colors)
    ax.add_collection3d(pc)
    
    if vnumbers:
        name = [verts[v].name for v in verts]
        for name, x, y, z in zip(name, xs, ys, zs):
            ax.text(x, y, z, name, None)
            
    if enumbers:
        e_xs = [es[e].center()[0] for e in es]
        e_ys = [es[e].center()[1] for e in es]
        e_zs = [es[e].center()[2] for e in es]
        e_name = [es[e].name for e in es]        
        for name, x, y, z in zip(e_name, e_xs, e_ys, e_zs):
            ax.text(x, y, z, name, None)
    
    if cnumbers:
        c_xs = [cs[c].center()[0] for c in cs]
        c_ys = [cs[c].center()[1] for c in cs]
        c_zs = [cs[c].center()[2] for c in cs]
        c_name = [cs[c].name for c in cs]        
        for name, x, y, z in zip(c_name, c_xs, c_ys, c_zs):
            ax.text(x, y, z, name, None, 
                    horizontalalignment='center',
                    verticalalignment='center')

    if not radius:
        ax.set_xlim(min(xs)-.1*np.mean(xs), max(xs)+.1*np.mean(xs))
        ax.set_ylim(min(ys)-.1*np.mean(ys), max(ys)+.1*np.mean(ys))
        ax.set_zlim(min(zs)-.1*np.mean(zs) - 0.0001, max(zs)+.1*np.mean(zs) + 0.0001)
    else:
        ax.set_xlim(center[0] - radius, center[0] + radius)
        ax.set_ylim(center[1] - radius, center[1] + radius)
        ax.set_zlim(-0.0001, 0.0001)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.scatter(xs, ys, zs)
    ax.view_init(elev = elev, azim = azim)
    
    if video:
        # Used to return the plot as an image array
        plt.axis('off')
        fig.canvas.draw()       # draw the canvas, cache the renderer
        image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        plt.close()
        return image
    else:
        plt.show()

def gif_from_data(path):
    # takes path to data folder, creates a gif
    plots = []
    end = 0
    for x in os.listdir(path):
        if '.pickle' in x:
            end = max(end, int(x[:-7]))
    end -= end%20
    
    for seq in range(0, end + 20, 20):
        file = path + '/' + str(seq) + '.pickle'
        pkl = open(file, "rb")
        vertices, edges, cells, rand = pickle.load(pkl)
        pkl.close()
        plots.append(plot_space(vertices, edges, cells, elev = 90, azim = 0,\
                                cnumbers = True, video = True,\
                                center = p.center, radius = p.radius))
    
    kwargs_write = {'fps':1/(p.stepsize*20), 'quantizer':'nq'}
    imageio.mimsave(os.getcwd().split(os.sep)[-1] + '_result.gif', plots, fps=10)

### Quantification functions ###

def n_stem_cells(cells):
    n = 0
    for c in cells:
        if cells[c].celltype == 'stem':
            n += 1
    return n

def dist_paneth_center(cell):
    return norm(cell.center() - p.center)

def avg_dist_paneth_center(cells):
    dist = []
    for c in cells:
        if cells[c].celltype == 'paneth':
            dist.append(norm(cells[c].center() - p.center))
    return np.mean(dist)

def avg_radial_distribution(cells, deltar):
    p_cells = []
    for c in cells:
        if cells[c].celltype == 'paneth':
            p_cells.append(cells[c])
    x = 0
    result = [[] for n in p_cells]
    while(x < p.radius):
        for i in range(len(p_cells)):
            n = 0
            for cell in p_cells:
                if cell is not p_cells[i]:
                    dist = norm(p_cells[i].center() - cell.center())
                    if x < dist and dist < x + deltar:
                        n += 1
            result[i].append(n)
        x += deltar
    return np.array(np.mean(np.matrix(result), axis = 0))[0]

### Main ###
if __name__ == '__main__':
    if not os.path.exists('data'):
        os.mkdir('data')
    else:
        print('output folder already exists, overwriting previous experiment')
    
    # Set the parameters
    p.init_from_file('parameters.pickle')
    
    # Initial configuration
    vertices, edges, cells = initiate_cells(p.radius, p.center, p.ncells)
    
    # get the system to steady state before doing any experiments
    # stops after the system reaches equilibrium or 500 timesteps have passed
    diff = - math.inf
    prev_energy = math.inf
    step = 0
    while(diff < -.1 and step < 500):
        vertices, edges, cells = update(vertices, edges, cells, step = p.stepsize)
        if step % 10 == 0:
            diff = system_energy(vertices, edges, cells) - prev_energy
            prev_energy = system_energy(vertices, edges, cells)
        step += 1
    
    # Make random cells within the middle third of the tissue into paneth or stem cells
    
    # Making sure there are at least 2 paneth cells
    paneth = []
    while(len(paneth) < 2):
        paneth = []
        for c in cells:
            if norm(cells[c].center() - p.center) < p.radius/2:
                if random() < .1:
                    paneth.append(c)
        
    for c in paneth:
        cells[c].celltype = 'paneth'
        cells[c].color = 'black'
        for cell in cells[c].adjacent_cells():
            cell.celltype = 'stem'
            cell.color = 'pink'
    
    for c in cells:
        if cells[c].dist_to_paneth(cells) > p.paneth_dist:
            cells[c].celltype = None
            cells[c].random_color()    
    
    print('initialisation completed')
    plot_space(vertices, edges, cells, elev = 90, azim = 0, cnumbers = True)
    
    for step in range(p.nsteps):
        vertices, edges, cells = update(vertices, edges, cells, step = p.stepsize)
        if step % 20 == 0:
            pkl = open('data/%d.pickle' % (step),'wb')
            pickle.dump([vertices, edges, cells, get_state()], pkl)
            pkl.close()
        if step % 1 == 0:
            print('step %d of %d' % (step, p.nsteps))
