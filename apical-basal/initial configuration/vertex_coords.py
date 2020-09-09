# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 09:35:51 2019

@author: simon
"""

import numpy as np
from vertex_objects import Vertex, Cell, Edge

# Get cell centers
cells = {}
w = 12
h = 6
z_offset = 0
n = 0
for x in range(w):
    for y in range(h):
        n += 1
        cells['c'+str(n)] = [x, y * 2, 0]
        n += 1
        cells['c'+str(n)] = [.5 + x, 1 + y * 2, 0]

vertices = {}
n = 1
for c in cells:
    a1 = cells[c] + np.array([-.5,-.5, 1])
    a2 = cells[c] + np.array([0,-.5, 1])
    a3 = cells[c] + np.array([.5,-.5, 1])
    a4 = cells[c] + np.array([.5,.5, 1])
    a5 = cells[c] + np.array([0,.5, 1])
    a6 = cells[c] + np.array([-.5,.5, 1])
    b1 = cells[c] + np.array([-.5,-.5, 0])
    b2 = cells[c] + np.array([0,-.5, 0])
    b3 = cells[c] + np.array([.5,-.5, 0])
    b4 = cells[c] + np.array([.5,.5, 0])
    b5 = cells[c] + np.array([0,.5, 0])
    b6 = cells[c] + np.array([-.5,.5, 0])
    if cells[c][0] == 0:
        new = [a2, a3, a4, a5, b2, b3, b4, b5]
    elif cells[c][0] == w-.5:
        new = [a1, a2, a5, a6, b1, b2, b5, b6]
    else:
        new = [a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, b6]
    cells[c] = Cell()
    for v in new:
        if any(v is x for x in [a1, a2, a3, a4, a5, a6]):
            vertex = Vertex(v, site='apical')
        else:
            vertex = Vertex(v, site='basal')
        # Check whether vertex already exists in vertices.
        if list(vertex.location) not in [list(v.location) for v in vertices.values()]:
            vertex.name = 'v' + str(n)
            vertices['v' + str(n)] = vertex
            cells[c].init_vertices.append(vertex)
            n += 1
        else:
            for v in vertices:
                if list(vertices[v].location) == list(vertex.location):
                    cells[c].init_vertices.append(vertices[v])

file = str(int(w/2)) + 'x' + str(h) + '.txt'
file = open(file, 'w')

for v in vertices:
    loc = vertices[v].location
    if loc[0] in [0, w - .5] or loc[1] in [-.5, h * 2 - .5]:
        bound = 'bound'
    else:
        bound = 'unbound'
        loc[2] -= z_offset
    file.write('vertex ' + v + ' ' 
               + str(loc[0]) + ' ' 
               + str(loc[1]) + ' ' 
               + str(loc[2]) + ' '+ bound +'\n')

edges = {}
n = 1
for c in cells:
    vs = [v.name for v in cells[c].init_vertices]
    cells[c].init_vertices = vs
    a_vs = []; b_vs = []
    for v in vs:
        if vertices[v].site == 'apical':
            a_vs.append(v)
        if vertices[v].site == 'basal':
            b_vs.append(v)
    # Apical edges
    for i in range(len(a_vs)):
        rolled = np.roll(a_vs, i)
        v1 = rolled[0]
        v2 = rolled[1]
        if [v2, v1] not in [edges[e].vertices for e in edges]:
            edges['a' + str(n)] = Edge([v1, v2])
            n += 1
    # Basal edges
    for i in range(len(b_vs)):
        rolled = np.roll(b_vs, i)
        v1 = rolled[0]
        v2 = rolled[1]
        if [v2, v1] not in [edges[e].vertices for e in edges]:
            edges['b' + str(n)] = Edge([v1, v2])
            n += 1
    # Lateral edges
    for i in range(int(len(vs)/2)):
        if [vs[i],vs[i+int(len(vs)/2)]] not in [edges[e].vertices for e in edges]:
            edges['l' + str(n)] = Edge([vs[i], vs[i+int(len(vs)/2)]])
            n += 1

n = 1
for edge in edges:
    file.write('edge ' + edge + ' ' + edges[edge].vertices[0] + ' ' + edges[edge].vertices[1] + '\n')
    n += 1

for c in cells:
    file.write('cell ' + c)
    vs = cells[c].init_vertices
    a_vs = []; b_vs = []
    for v in vs:
        if vertices[v].site == 'apical':
            a_vs.append(v)
        if vertices[v].site == 'basal':
            b_vs.append(v)
    for i in range(len(a_vs)):
            rolled = np.roll(a_vs, i)
            for e in edges:
                if [rolled[0], rolled[1]] == edges[e].vertices:
                    file.write(' ' + e + ' ' + '1')
                if [rolled[1], rolled[0]] == edges[e].vertices:
                    file.write(' ' + e + ' ' + '-1')
    # basal are in the opposite orientation
    for i in range(len(b_vs)):
            rolled = np.roll(b_vs, i)
            for e in edges:
                if [rolled[0], rolled[1]] == edges[e].vertices:
                    file.write(' ' + e + ' ' + '-1')
                if [rolled[1], rolled[0]] == edges[e].vertices:
                    file.write(' ' + e + ' ' + '1')
    
    for i in range(int(len(vs)/2)):
        for e in edges:
            if [vs[i], vs[i + int(len(vs)/2)]] == edges[e].vertices:
                file.write(' ' + e)
    file.write('\n')

file.close()