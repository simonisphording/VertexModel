# -*- coding: utf-8 -*-
"""
Created on Tue May  7 10:15:26 2019

@author: isphording
"""

from vertex_objects import Vertex, Edge, Cell
from scipy.spatial import Voronoi
from numpy.linalg import norm
import numpy as np
import random
import math

def rand_point_in_circle(radius, center):
    # random angle
    alpha = 2 * math.pi * random.random()
    # random radius    
    u = random.random() + random.random()
    r = radius * (2 - u if u > 1 else u)
    # return coordinates
    return np.array([r * math.cos(alpha) + center[0], r * math.sin(alpha) + center[1]])

def convert_to_dict(l):
    result = {}
    for item in l:
        result[item.name] = item
    return result

def initiate_cells(radius, center, ncells):
    # Create cells, edges and vertices using Voronoi Tessellation
    points =  [rand_point_in_circle(radius, center) for i in range(ncells)]
    vor = Voronoi(points)
    
    vertices = []
    i = 1
    for vertex in vor.vertices:
        vnew = Vertex(vertex, name = 'v' + str(i))
        vnew.location = np.append(vnew.location, 0)
        
        vertices.append(vnew)
        i += 1
    
    edges = []
    i = 1
    for edge in vor.ridge_vertices:
        if -1 not in edge:
            v1 = vertices[edge[0]]; v2 = vertices[edge[1]]
            # Only create edges within the radius
            if norm(center - v1.location) < radius and norm(center - v2.location) < radius:
                enew = Edge([v1, v2], name = 'e' + str(i))
                v1.add_edge(enew); v2.add_edge(enew)
                edges.append(enew)
                i += 1
    
    cells = []
    i = 1
    for cell in vor.regions:
        if -1 not in cell and cell:
            # Only create cells with all vertices within the radius
            if all([norm(center - vertices[v].location) < radius for v in cell]):
                cnew = Cell(name = 'c' + str(i))
                cnew.random_color()
                for x in range(len(cell)):
                    vs = np.roll(cell, x)
                    v1 = vertices[vs[0]]; v2 = vertices[vs[1]]
                    for edge in edges:
                        if edge.vertices[0] == v1 and edge.vertices[1] == v2:
                            cnew.edges.append(edge)
                            cnew.direction.append(1)
                            edge.add_cell(cnew)
                        if edge.vertices[1] == v1 and edge.vertices[0] == v2:
                            cnew.edges.append(edge)
                            cnew.direction.append(-1)
                            edge.add_cell(cnew)
                    v1.add_cell(cnew)
                cells.append(cnew)
                i +=1
    
    # Remove vertices that are not connected to the rest of the system
    for v in list(vertices):
        if norm(center - v.location) > radius or not v.cells:
            vertices.remove(v)
    
    vertices = convert_to_dict(vertices)
    edges = convert_to_dict(edges)
    cells = convert_to_dict(cells)
    
    return vertices, edges, cells