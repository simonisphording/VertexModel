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

def t1(e, vertices, edges, cells):
    # Do a t1 transition, replacing edge e and its sister edge
    a = e.vertices[0]
    b = e.vertices[1]
    if len(a.cells) > 2 and len(b.cells) > 2:
        c1 = e.cells[0]
        c2 = e.cells[1]
        if len(c1.apical_edges) > 3 and len(c2.apical_edges) > 3 and len(a.cells) > 2 and len(b.cells) > 2:
            sister = e.sister_edge()
            v1, v2, o1 = t1_single_edge(e, vertices, edges)
            v3, v4, o2 = t1_single_edge(sister, vertices, edges)
            # Delete old lateral edges
            for edge in o1:
                for cell in edge.cells:
                    cell.remove_edge(edge, 'lateral')
                del edges[edge.name]
            
            # Check which pairs of vertices are connected laterally
            if all([c in v1.cells for c in v3.cells]):
                pairs = [[v1,v3], [v2, v4]]
            else:
                pairs = [[v1,v4], [v2, v3]]
            for pair in pairs:
                # Create new lateral edge
                index = 'l' + str(int(list(edges.keys())[-1][1:])+1)
                enew = Edge(pair, cells = pair[0].cells, name = index)
                edges[index] = enew
                # Add lateral edge to new vertices
                for v in pair:
                    v.add_edge(enew)
                # Add lateral edge to cells
                for c in enew.cells:
                    c.add_edge(enew, 'lateral')
    return vertices, edges, cells

def t1_single_edge(e, vertices, edges):
    # Do a t1 transition, but only on either the apical or the basal side
    v1 = e.vertices[0]
    v2 = e.vertices[1]
    
    old_edges = []
    for vert in [v1, v2]:
        for edge in vert.edges:
            if edge.is_lateral():
                old_edges.append(edge)
    
    c1 = e.cells[0]
    c2 = e.cells[1]
    fourfold = [] # List of 4 edges connected to v1 and v2
    for edge in v1.edges:
        if edge != e and not edge.is_lateral():
            fourfold.append(edge)
    for edge in v2.edges:
        if edge != e and not edge.is_lateral():
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
    if e.is_apical():
        pro = 'a'
    else:
        pro = 'b'
    index = pro + str(int(list(edges.keys())[-1][1:])+1)
    enew = Edge([v1_new, v2_new], cells = other_cells, name = index)
    edges[index] = enew
    # Connect new edge with 4 neighbouring edges, checking the orientation of enew
    # Calc distance between vertices and new edge in both orientations to see which oune has the shortest distances
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
    if e.is_apical():
        site = 'apical'
    else:
        site = 'basal'
    for cell in [c1, c2]:
        cell.remove_edge(e, site)
    for cell in other_cells:
        cell.add_edge(enew, site)
    
    # Deleting old lateral edge
    old_edges = []
    for vert in [v1, v2]:
        for edge in vert.edges:
            if edge.is_lateral():
                old_edges.append(edge)
        
    # Delete old objects
    del vertices[v1.name]
    del vertices[v2.name]
    del edges[e.name]
    return [v1_new, v2_new, old_edges]

def split_edge(e, location, site, vertices, edges):
    # Create a new vertex on given location, splitting edge e
    index = 'v' + str(int(list(vertices.keys())[-1][1:])+1)
    vnew = Vertex(location, name = index, cells = list(e.cells))
    vertices[index] = vnew
    
    if e.is_apical():
        pre = 'a'
    else:
        pre = 'b'
    
    # Create two new edges from the original
    index = pre + str(int(list(edges.keys())[-1][1:])+1)
    e1 = Edge([e.vertices[0], vnew], cells = list(e.cells), name = index)
    edges[index] = e1
    for vertex in e1.vertices:
        vertex.add_edge(e1)
        vertex.remove_edge(e)
    
    index = pre + str(int(list(edges.keys())[-1][1:])+1)
    e2 = Edge([vnew, e.vertices[1]], cells = list(e.cells), name = index)
    edges[index] = e2
    for vertex in e2.vertices:
        vertex.add_edge(e2)
        vertex.remove_edge(e)
    
    # Replace original edge with new edges
    for cell in e.cells:
        cell.remove_edge(e, site)
        cell.add_edge(e1, site)
        cell.add_edge(e2, site)
    
    del edges[e.name]

    return vnew, vertices, edges

def edge_from_length(cell, edges, length):
    # Given a length, return a matching location on the perimeter of a cell and the edge in lies on
    edges = cell.ordered_edges(edges)
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

def vol_diff(x, cell, e1, loc1, len1):
    x = x * cell.perimeter(cell.apical_edges)
    if x < cell.perimeter(cell.apical_edges) - len1:
        y = x + len1
    else:
        y = x + len1 - cell.perimeter(cell.apical_edges)
    e2, loc2 = edge_from_length(cell, cell.apical_edges, y)
    # If e1 and e2 are the same edge, one of the two resulting cells has all the cell volume
    if e1 == e2:
        # Dependent on locations, return - vol cell or + vol cell
        if norm(loc2 - cell.directed_edge(e1)[0].location) < norm(loc1 - cell.directed_edge(e1)[0].location):
            return cell.surface(cell.apical_edges)
        else:
            return -1 * cell.surface(cell.apical_edges)
    else:
        return mock_divide(cell, e1, e2, loc1, loc2)

def mock_divide(cell, e1, e2, loc1, loc2):
    # Tries out a split of the cell surface and returns the resulting cell surfaces
    verts = [loc1]
    e = e1
    while(e is not e2):
        v = cell.directed_edge(e)[1]
        verts.append(v.location)
        for edge in cell.apical_edges:
            if cell.directed_edge(edge)[0] == v:
                e = edge
    verts.append(loc2)
    area = np.array([0. for dim in verts[0]])
    n = len(verts)
    for i in range(n):
        verts = np.roll(verts, 1, axis = 0)
        area += np.cross(verts[0], verts[1])
    return norm(area/2) - cell.surface(cell.apical_edges)/2
    

def divide_cell(c, es, locs, vertices, edges, cells):
    # Divides a cell c given the locations of new vertices (locs) and the edges they lie on (es)
    v1, vertices, edges = split_edge(es[0], locs[0], 'apical', vertices, edges)
    v2, vertices, edges = split_edge(es[1], locs[1], 'apical', vertices, edges)
    v3, vertices, edges = split_edge(es[2], locs[2], 'basal', vertices, edges)
    v4, vertices, edges = split_edge(es[3], locs[3], 'basal', vertices, edges)
    
    # Create a new edge between v1 and v2
    index = 'a' + str(int(list(edges.keys())[-1][1:])+1)
    e12 = Edge([v1, v2], name = index)
    edges[index] = e12
    for vertex in e12.vertices:
        vertex.add_edge(e12)
    
    # Create a new edge between v3 and v4
    index = 'b' + str(int(list(edges.keys())[-1][1:])+1)
    e34 = Edge([v3, v4], name = index)
    edges[index] = e34
    for vertex in e34.vertices:
        vertex.add_edge(e34)   
    
    # Create new lateral edges
    pairs = [[v1, v3],[v2, v4]]
    new_laterals = []
    for pair in pairs:
        # Create new lateral edge
        index = 'l' + str(int(list(edges.keys())[-1][1:])+1)
        enew = Edge(pair, cells = list(pair[0].cells), name = index)
        edges[index] = enew
        # Add lateral edge to new vertices
        for v in pair:
            v.add_edge(enew)
        # Add lateral edge to cells
        for cell in enew.cells:
            cell.add_edge(enew, 'lateral')
        new_laterals.append(enew)
    
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
                if edge in c.apical_edges:
                    if c.directed_edge(edge)[0] == vnext:
                        cnew.apical_edges.append(edge)
                        cnew.apical_direction.append(c.edge_direction(edge))
                        edge.add_cell(cnew)
                        edge.remove_cell(c)
                        vnext.add_cell(cnew)
                        vnext.remove_cell(c)
                        vnext = c.directed_edge(edge)[1]
        b.add_cell(cnew)
        b.remove_cell(c)
        
        # Add the new edge
        cnew.add_edge(e12, 'apical')
        e12.add_cell(cnew)
        cells[index] = cnew
        
        # Add lateral edges
        for vertex in cnew.vertices(cnew.apical_edges):
            for edge in vertex.edges:
                if edge.is_lateral() or edge in new_laterals:
                    cnew.add_edge(edge, 'lateral')
                    edge.cells = list(edge.vertices[0].cells)
        
                        
        # Add basal edges
        apical_verts = cnew.ordered_vertices(cnew.apical_edges)
        basal_verts = []
        for vertex in apical_verts:
            for edge in vertex.edges:
                if edge.is_lateral():
                    basal_verts.append(edge.other_vertex(vertex))
        basal_verts.reverse()
        
        for i in range(len(basal_verts)):
            basal_verts = np.roll(basal_verts, 1)
            va, vb = basal_verts[0:2]
            for edge in va.edges:
                if edge.vertices[0] == va and edge.vertices[1] == vb:
                    cnew.basal_edges.append(edge)
                    cnew.basal_direction.append(1)
                    edge.add_cell(cnew)
                    edge.remove_cell(c)
                if edge.vertices[0] == vb and edge.vertices[1] == va:
                    cnew.basal_edges.append(edge)
                    cnew.basal_direction.append(-1)
                    edge.add_cell(cnew)
                    edge.remove_cell(c)
        
        for vertex in basal_verts:
            vertex.add_cell(cnew)
            vertex.remove_cell(c)

        new_cells.append(cnew)
    
    del cells[c.name]
    
    return new_cells, vertices, edges, cells

def t3(c, vertices, edges, cells):
    # Cell division, resulting in daughter cells with exact halves of the mother cell's apical surface
    peri = c.perimeter(c.apical_edges)
    rpoint = random() * peri
    e1, loc1 = edge_from_length(c, c.apical_edges, rpoint)
    sol = brentq(vol_diff, .1, .9, args = (c, e1, loc1, rpoint))
    x = sol * peri
    if x < peri - rpoint:
        y = x + rpoint
    else:
        y = x + rpoint - peri
    e2, loc2 = edge_from_length(c, c.apical_edges, y)
    
    # Get e3, loc3, e4 and loc4
    e3 = e1.sister_edge()
    e4 = e2.sister_edge()
    dist = norm(c.directed_edge(e1)[0].location - loc1) / e1.length()
    b, a = [v.location for v in c.directed_edge(e3)]
    loc3 = a + dist * e3.length() * (b-a)/norm(b-a)
    dist = norm(c.directed_edge(e2)[0].location - loc2) / e2.length()
    b, a = [v.location for v in c.directed_edge(e4)]    
    loc4 = a + dist * e4.length() * (b-a)/norm(b-a)
    
    new_cells, vertices, edges, cells = divide_cell(c, [e1,e2,e3,e4], [loc1, loc2, loc3, loc4], vertices, edges, cells)
    return new_cells, vertices, edges, cells