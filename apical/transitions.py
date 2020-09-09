# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:05:24 2019

@author: isphording
"""

import numpy as np
from numpy.linalg import norm
from numpy.random import random
from scipy.optimize import brentq
import parameters as p
from vertex_objects import Vertex, Cell, Edge, gs

def t1(e, vertices, edges):
    v1 = e.vertices[0]
    v2 = e.vertices[1]
    if len(v1.cells) > 2 and len(v2.cells) > 2:
        c1 = e.cells[0]
        c2 = e.cells[1]
        if len(c1.edges) > 3 and len(c2.edges) > 3:
            fourfold = [] # List of 4 edges connected to v1 and v2
            for edge in v1.edges:
                if edge != e:
                    fourfold.append(edge)
            for edge in v2.edges:
                if edge != e:
                    fourfold.append(edge)
        
            # Get enew, perpendicular to e (see methods for explanation)
            c1_edges = []
            c2_edges = []
            for edge in fourfold:
                if c1 in edge.cells:
                    c1_edges.append(edge)
                if c2 in edge.cells:
                    c2_edges.append(edge)
            
            e1 = np.array([0. for dim in v1.location])
            e2 = np.array([0. for dim in v1.location])
            
            for edge in c1_edges:
                if v1 in edge.vertices:
                    e1 += edge.other_vertex(v1).location - v1.location
                else:
                    e1 += edge.other_vertex(v2).location - v2.location
            for edge in c2_edges:
                if v1 in edge.vertices:
                    e2 += edge.other_vertex(v1).location - v1.location
                else:
                    e2 += edge.other_vertex(v2).location - v2.location
            
            u = v1.location - v2.location
            e1 = gs(u, e1)
            e2 = gs(u, e2)
            u_ = np.cross(u, e1 + e2)
            if np.array_equal(np.around(e1, decimals = 5), np.around(e2, decimals = 5))\
            or np.array_equal(np.around(e1, decimals = 5), np.around(-e2, decimals = 5)):
                u_ = e1/norm(e1)
            else:
                u_ = u_/norm(u_)
            r1 = e.center() + .5 * p.new_edge_length * u_
            r2 = e.center() - .5 * p.new_edge_length * u_
            
            # List of cells involved but not connected to e
            other_cells = set()
            for edge in fourfold:
                for cell in edge.cells:
                    if cell not in [c1, c2]:
                        other_cells.add(cell)
            other_cells = list(other_cells)
            
            # Creating new vertices
            index = 'v' + str(int(list(vertices.keys())[-1][1:])+1)
            v1_new = Vertex(r1, name = index)
            vertices[index] = v1_new
            index = 'v' + str(int(list(vertices.keys())[-1][1:])+1)
            v2_new = Vertex(r2, name = index)
            vertices[index] = v2_new
            # Creating new edge
            index = 'e' + str(int(list(edges.keys())[-1][1:])+1)
            enew = Edge([v1_new, v2_new], cells = other_cells, name = index)
            edges[index] = enew
                
            # Connect new edge with 4 neighbouring edges, checking the orientation of enew
                
            L1 = 0; L2 = 0
            
            for edge in c1_edges:
                if v1 in edge.vertices:
                    L1 += norm(v1_new.location - edge.other_vertex(v1).location)
                    L2 += norm(v2_new.location - edge.other_vertex(v1).location)
                else:
                    L1 += norm(v1_new.location - edge.other_vertex(v2).location)
                    L2 += norm(v2_new.location - edge.other_vertex(v2).location)
            for edge in c2_edges:
                if v1 in edge.vertices:
                    L1 += norm(v2_new.location - edge.other_vertex(v1).location)
                    L2 += norm(v1_new.location - edge.other_vertex(v1).location)
                else:
                    L1 += norm(v2_new.location - edge.other_vertex(v2).location)
                    L2 += norm(v1_new.location - edge.other_vertex(v2).location)
            
            if L1 < L2:
                va = v1_new
                vb = v2_new
            else:
                va = v2_new
                vb = v1_new
            
            for edge in c1_edges:
                if v1 in edge.vertices:
                    edge.vertices[edge.vertices.index(v1)] = va
                else:
                    edge.vertices[edge.vertices.index(v2)] = va
                va.add_edge(edge)
            for edge in c2_edges:
                if v1 in edge.vertices:
                    edge.vertices[edge.vertices.index(v1)] = vb
                else:
                    edge.vertices[edge.vertices.index(v2)] = vb
                vb.add_edge(edge)
            
            v1_new.add_edge(enew)
            v2_new.add_edge(enew)
            
            # Set cell index for new vertices
            for cell in other_cells:
                va.add_cell(cell)
                vb.add_cell(cell)
            va.add_cell(c1)
            vb.add_cell(c2)
            
            # Changing edge indeces of cells
            for cell in [c1, c2]:
                cell.remove_edge(e)
            for cell in other_cells:
                cell.add_edge(enew)
        
            # Delete old objects
            del vertices[v1.name]
            del vertices[v2.name]
            del edges[e.name]

def t2(c, vertices, edges, cells):
    # Apoptosis
    
    # Create new vertex at the center
    index = 'v' + str(int(list(vertices.keys())[-1][1:])+1)
    vnew = Vertex(c.center(), name = index)
    vertices[index] = vnew
    
    # Redirect all edges connected to c's vertices to vnew and delete vertices
    for vertex in c.vertices():
        for edge in vertex.edges:
            if edge not in c.edges:
                i = edge.vertices.index(vertex)
                edge.vertices[i] = vnew
                vnew.add_edge(edge)
                for cell in edge.cells:
                    vnew.add_cell(cell)
                
    # Delete old vertices and the old cell
    for vertex in c.vertices():
        del vertices[vertex.name]
    
    # Remove old edges
    for edge in c.edges:
        for cell in edge.cells:
            if cell is not c:
                cell.remove_edge(edge)
        del edges[edge.name]
    
    del cells[c.name]

def split_edge(e, location, vertices, edges):
    # Create a new vertex in the middle of the edge
    index = 'v' + str(int(list(vertices.keys())[-1][1:])+1)
    vnew = Vertex(location, name = index, cells = list(e.cells))
    vertices[index] = vnew
    
    # Create two new edges from the original
    index = 'e' + str(int(list(edges.keys())[-1][1:])+1)
    e1 = Edge([e.vertices[0], vnew], cells = list(e.cells), name = index)
    edges[index] = e1
    for vertex in e1.vertices:
        vertex.add_edge(e1)
        vertex.remove_edge(e)
    
    index = 'e' + str(int(list(edges.keys())[-1][1:])+1)
    e2 = Edge([vnew, e.vertices[1]], cells = list(e.cells), name = index)
    edges[index] = e2
    for vertex in e2.vertices:
        vertex.add_edge(e2)
        vertex.remove_edge(e)
    
    # Replace original edge with new edges
    for cell in e.cells:
        cell.remove_edge(e)
        cell.add_edge(e1)
        cell.add_edge(e2)
    
    del edges[e.name]

    return vnew, vertices, edges

def edge_from_length(cell, length):
    # Given a length along the perimeter of the cell, 
    # return an edge and a location on that edge accordingly
    edges = cell.ordered_edges()
    l = 0
    i = 0
    edge = None
    while edge == None:
        e = edges[i]
        if l + e.length() <= length:
            l += e.length()
        else:
            edge = e
            a = cell.directed_edge(edge)[0].location
            b = cell.directed_edge(edge)[1].location
            loc = a + (length - l) * (b-a)/norm(b-a)
        i += 1
    return edge, loc

def vol_diff(x, cell, e1, loc1, len1, vertices, edges, cells):
    # Given a location on the cell perimeter, calculate the difference in the daughter cell's volume
    x = x * cell.perimeter()
    if x < cell.perimeter() - len1:
        y = x + len1
    else:
        y = x + len1 - cell.perimeter()
    e2, loc2 = edge_from_length(cell, y)
    if e1 == e2:
        # Dependent on locations, return - vol cell or + vol cell
        if norm(loc2 - cell.directed_edge(e1)[0].location) < norm(loc1 - cell.directed_edge(e1)[0].location):
            return cell.surface()
        else:
            return -1 * cell.surface()
    else:
        return mock_divide(cell, e1, e2, loc1, loc2)

def mock_divide(cell, e1, e2, loc1, loc2):
    # Calculates volume difference of cells without actually changing the system
    verts = [loc1]
    e = e1
    while(e is not e2):
        v = cell.directed_edge(e)[1]
        verts.append(v.location)
        for edge in cell.edges:
            if cell.directed_edge(edge)[0] == v:
                e = edge
    verts.append(loc2)
    area = np.array([0. for dim in verts[0]])
    n = len(verts)
    for i in range(n):
        verts = np.roll(verts, 1, axis = 0)
        area += np.cross(verts[0], verts[1])
    return norm(area/2) - cell.surface()/2
    

def divide_cell(c, e1, e2, loc1, loc2, vertices, edges, cells):
    v1, vertices, edges = split_edge(e1, loc1, vertices, edges)
    v2, vertices, edges = split_edge(e2, loc2, vertices, edges)
    
    # Create a new edge between v1 and v2
    index = 'e' + str(int(list(edges.keys())[-1][1:])+1)
    enew = Edge([v1, v2], name = index)
    edges[index] = enew
    for vertex in enew.vertices:
        vertex.add_edge(enew)

    # Get the edges connecting v1 to v2 and directions. Same for v2 v1.
    # Get next edge from v1, find vnext, repeat until v2 is reached...
    new_cells = []
    for a in [v1, v2]:
        index = 'c' + str(int(list(cells.keys())[-1][1:])+1)
        cnew = Cell(name = index, color = c.color, celltype = c.celltype, target_volume = c.target_volume/2)
        vnext = a
        if vnext == v1:
            b = v2
        else:
            b = v1
        while vnext is not b:
            for edge in vnext.edges:
                # Check whether edge points away from vnext
                if edge in c.edges:
                    if c.directed_edge(edge)[0] == vnext:
                        cnew.edges.append(edge)
                        cnew.direction.append(c.edge_direction(edge))
                        edge.add_cell(cnew)
                        edge.remove_cell(c)
                        vnext.add_cell(cnew)
                        vnext.remove_cell(c)
                        vnext = c.directed_edge(edge)[1]
        b.add_cell(cnew)
        b.remove_cell(c)
        # Finally, add the new edge
        cnew.add_edge(enew)
        enew.add_cell(cnew)
        cells[index] = cnew
        new_cells.append(cnew)
    del cells[c.name]
    
    return new_cells, vertices, edges, cells

def t3(c, vertices, edges, cells):
    # Cell division with exact halfs of surface
    peri = c.perimeter()
    rpoint = random() * peri
    e1, loc1 = edge_from_length(c, rpoint)
    # Determine best second point on edge to split the cells at
    sol = brentq(vol_diff, .1, .9, args = (c, e1, loc1, rpoint, vertices, edges, cells))
    x = sol * peri
    if x < peri - rpoint:
        y = x + rpoint
    else:
        y = x + rpoint - peri
    e2, loc2 = edge_from_length(c, y)
    # Split the mother cell into daughter cells
    new_cells, vertices, edges, cells = divide_cell(c, e1, e2, loc1, loc2, vertices, edges, cells)
    return new_cells, vertices, edges, cells
