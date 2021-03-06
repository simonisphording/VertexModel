# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:50:38 2019

@author: isphording
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d
from numpy.random import random, get_state
from vertex_objects import Vertex, Cell, Edge
from transitions import t1, t3
import parameters as p
import sys, os
import imageio
import pickle
sys.setrecursionlimit(100000)

### Model functions ###
    
def system_energy(vertices, edges, cells):
    # Returns the free energy of the system
    surface = 0
    volume = 0
    line = 0
    
    for c in cells:
        cell = cells[c]
        surface += cell.apical_tension * cell.surface(cell.apical_edges)
        surface += cell.basal_tension * cell.surface(cell.basal_edges)
        for surface_edges in cell.lateral_surfaces():
            surface += p.lateral_tension * cell.surface(surface_edges)
        volume += p.compression * (cell.volume() - cell.target_volume)**2
    
    for e in edges:
        line += p.tension * edges[e].length()
    
    return surface + volume + line

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
        edge = edges[e]
        if edge.is_apical() and edge.length() + edge.sister_edge().length() < p.transition_boundary:
            vertices, edges, cells = t1(edges[e], vertices, edges, cells)
    for c in list(cells):
        if cells[c].celltype == 'stem':
            cells[c].target_volume += cells[c].target_volume * p.growthrate
            if cells[c].volume() > 2 * p.init_volume:
                new_cells, vertices, edges, cells = t3(cells[c], vertices, edges, cells)
                # Differentiation of newly formed stem cells
                for c in new_cells:
                    if not any([x.celltype == 'paneth' for x in c.adjacent_cells()]): #c.dist_to_paneth(cells) > p.paneth_dist:
                        c.celltype = None
                        c.random_color()
    return(vertices, edges, cells)

### Parsing and plotting ###

def parse_setup(file, silent = False):
    # Parses initial configuration from a text file, ie. generated by round_sheet
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
            i = 2
            while(i < len(line)):
                edges[line[i]].add_cell(cell)
                if line[i][0] == 'a':
                    cell.apical_edges.append(edges[line[i]])
                    cell.apical_direction.append(int(line[i+1]))
                    i += 2
                elif line[i][0] == 'b':
                    cell.basal_edges.append(edges[line[i]])
                    cell.basal_direction.append(int(line[i+1]))
                    i += 2
                elif line[i][0] == 'l':
                    cell.lateral_edges.append(edges[line[i]])
                    i += 1
                else:
                    print('incorrect entry!!!')
                    i += 1
            # Adding cell to vertices
            for site in [cell.apical_edges, cell.basal_edges, cell.lateral_edges]:
                for edge in site:
                    for v in edge.vertices:
                        v.add_cell(cell)
        else:
             if silent == False: 
                print('line does not start with vertex, edge or cell')
    return (vertices, edges, cells)

def plot_space(verts, es, cs, elev = 10, azim = 120, video = False,\
               vnumbers = False, enumbers = False, cnumbers = False, forces = False):
    # Plots the simulation space
    
    if video:
        fig = plt.figure(figsize = (20,20))
    else:
        fig = plt.figure(figsize = (10,10))
    
    ax = fig.add_subplot(111, projection = '3d')

    xs = [verts[v].location[0] for v in verts]
    ys = [verts[v].location[1] for v in verts]
    zs = [verts[v].location[2] for v in verts]
    
    patches = []; colors = []
    for c in cs:
        patches.append([vertex.location for vertex in cs[c].ordered_vertices(cs[c].apical_edges)])
        colors.append(cs[c].color)
        patches.append([vertex.location for vertex in cs[c].ordered_vertices(cs[c].basal_edges)])
        colors.append(cs[c].color)
    pc = Poly3DCollection(patches, facecolors=colors)
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
        c_xs = [cs[c].center(cs[c].ordered_vertices(cs[c].apical_edges))[0] for c in cs]
        c_ys = [cs[c].center(cs[c].ordered_vertices(cs[c].apical_edges))[1] for c in cs]
        c_zs = [cs[c].center(cs[c].ordered_vertices(cs[c].apical_edges))[2] for c in cs]
        c_name = [cs[c].name for c in cs]        
        for name, x, y, z in zip(c_name, c_xs, c_ys, c_zs):
            ax.text(x, y, z, name, None)
    
    ax.set_xlim(min(xs)-.1*np.mean(xs), max(xs)+.1*np.mean(xs))
    ax.set_ylim(min(ys)-.1*np.mean(ys), max(ys)+.1*np.mean(ys))
    #ax.set_zlim(min(zs)-.1*np.mean(zs), max(zs)+.1*np.mean(zs))
    ax.set_zlim(-2, 4)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.scatter(xs, ys, zs)
    ax.view_init(elev = elev, azim = azim)
    plt.draw()
    
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

def gif_from_data(path, elev = 10, azim = 120):
    # Takes path to data folder, creates a gif
    stepsPerFrame = int(2/p.stepsize)
    plots = []
    end = 0
    for x in os.listdir(path):
        if '.pickle' in x:
            end = max(end, int(x[:-7]))
    end -= end%stepsPerFrame
    
    for seq in range(0, end + stepsPerFrame, stepsPerFrame):
        file = path + '/' + str(seq) + '.pickle'
        pkl = open(file, "rb")
        vertices, edges, cells, rand = pickle.load(pkl)
        pkl.close()
        plots.append(plot_space(vertices, edges, cells, \
                                elev = elev, azim = azim, video = True))
    
    kwargs_write = {'fps': 10, 'duration': p.stepsize}
    imageio.mimsave(os.getcwd().split(os.sep)[-1] + '_result.gif', plots, **kwargs_write)

### Main ###
if __name__ == '__main__':
    
    # Create a folder for .pickle outputs
    if not os.path.exists('data'):
        os.mkdir('data')
    else:
        print('output folder already exists, overwriting previous experiment')
    
    # Set the parameters
    p.init_from_file('parameters.pickle')
      
    # Set initial model topology    
    vertices, edges, cells = parse_setup('4_rings.txt', silent = True)
    
    # Set some cells to be paneth or stem cells
    for c in ['c1']:
        cells[c].celltype = 'paneth'
        cells[c].color = 'black'
        cells[c].apical_tension = .5
        cells[c].basal_tension = .1
        for cell in cells[c].adjacent_cells():
            cell.celltype = 'stem'
            cell.color = 'pink'
    
    # The outer vertices are bound to their location
    for v in vertices:
        if len(vertices[v].cells) < 3:
            vertices[v].bound = True
    
    for step in range(p.nsteps):
        vertices, edges, cells = update(vertices, edges, cells, step = p.stepsize)
        if step % 2/p.stepsize == 0:
            pkl = open('data/%d.pickle' % (step),'wb')
            pickle.dump([vertices, edges, cells, get_state()], pkl)
            pkl.close()

        if(step % 10 == 0):
            print(step, system_energy(vertices, edges, cells))