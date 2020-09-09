# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 10:30:27 2019

@author: isphording
"""

import numpy as np
from numpy.linalg import norm
import parameters as p
import random
import math

class Vertex:
    
    def __init__(self, location, edges = None, cells = None,\
                 site = None, bound = False, name = None):
        self.location = np.array(location)
        if edges is None:
            self.edges = []
        else:
            self.edges = edges
        if cells is None:
            self.cells = []
        else:
            self.cells = cells
        self.name = name
        self.bound = bound
        self.site = site
    
    def force(self):
        # Calculates the force on a vertex
        # The force consists of surface tensions, volume conservation and line tensions
        surface = [0 for dimension in self.location]
        volume = [0 for dimension in self.location]
        pressure = [0 for dimension in self.location]
        for cell in self.cells:
            if self in cell.vertices(cell.apical_edges):
                surface += -1 * cell.apical_tension * self.surface_gradient(cell, cell.apical_edges)
                pressure += -1 * p.pressure * cell.surface_vector(cell.apical_edges) / len(cell.apical_edges)
            else:
                surface += -1 * cell.basal_tension * self.surface_gradient(cell, cell.basal_edges)
            for surface_edges in cell.lateral_surfaces():
                if self in cell.vertices(surface_edges):
                    surface += -1 * p.lateral_tension * self.surface_gradient(cell, surface_edges)
            volume += 1 * p.compression * (cell.volume() - cell.target_volume) * self.volume_gradient(cell)
        
        line = [0 for dimension in self.location]
        for edge in self.edges:
            line += -1 * p.tension * edge.gradient(vertex = self)
        return surface + line + volume + pressure
    
    def volume_gradient(self, cell, debug = False):
        # Returns the gradient of the volume of a cell
        result = [0 for dim in self.location]
        # for each surface, get ordered_vertices
        surfaces = cell.lateral_surfaces()
        surfaces.append(cell.apical_edges)
        surfaces.append(cell.basal_edges)
        for surface in surfaces:
            vertices = cell.ordered_vertices(surface)
            c = cell.center(vertices)
            n = len(vertices)
            for vertex in vertices:
                i = vertices.index(vertex)
                rolled = np.roll(vertices, 1-i)
                r0 = rolled[0].location
                r1 = rolled[1].location
                r2 = rolled[2].location
                a1 = triangle_area_vector(r0, r1, c)
                a2 = triangle_area_vector(r1, r2, c)
                if vertex == self:
                    # A_1
                    x = np.dot([0, .5 * (- (r0[2] - c[2]) - (r1[2] - r0[2])/n), .5 * ((r0[1] - c[1]) + (r1[1] - r0[1])/n)], c) + a1[0]/n
                    y = np.dot([.5 * ((r0[2] - c[2]) + (r1[2] - r0[2])/n), 0, .5 * (- (r0[0] - c[0]) - (r1[0] - r0[0])/n)], c) + a1[1]/n
                    z = np.dot([.5 * (- (r0[1] - c[1]) - (r1[1] - r0[1])/n), .5 * ((r0[0] - c[0]) + (r1[0] - r0[0])/n), 0], c) + a1[2]/n
                    result += np.array([x,y,z])/3
                    # A_2
                    x = np.dot([0, .5 * ((r2[2] - c[2]) + (r1[2] - r2[2])/n), .5 * (- (r2[1] - c[1]) - (r1[1] - r2[1])/n)], c) + a2[0]/n
                    y = np.dot([.5 * (- (r2[2] - c[2]) - (r1[2] - r2[2])/n), 0, .5 * ((r2[0] - c[0]) + (r1[0] - r2[0])/n)], c) + a2[1]/n
                    z = np.dot([.5 * ((r2[1] - c[1]) + (r1[1] - r2[1])/n), .5 * (- (r2[0] - c[0]) - (r1[0] - r2[0])/n), 0], c) + a2[2]/n
                    result += np.array([x,y,z])/3
                elif rolled[0] != self:
                    # A_{1,2,j}
                    x = np.dot([0, .5/n * - (r1[2] - r0[2]), .5/n * (r1[1] - r0[1])], c) + a1[0]/n
                    y = np.dot([.5/n * (r1[2] - r0[2]), 0, .5/n * - (r1[0] - r0[0])], c) + a1[1]/n
                    z = np.dot([.5/n * - (r1[1] - r0[1]), .5/n * (r1[0] - r0[0]), 0], c) + a1[2]/n
                    result += np.array([x,y,z])/3
            if debug:
                print([edge.name for edge in surface], np.array([x,y,z])/3)
        return result
    
    def surface_gradient(self, cell, edges):
        # Returns the gradient of a given surface
        vertices = cell.ordered_vertices(edges)
        result = [0 for dim in self.location]
        c = cell.center(vertices)
        n = len(vertices)
        for vertex in vertices:
            i = vertices.index(vertex)
            rolled = np.roll(vertices, 1-i)
            r0 = rolled[0].location
            r1 = rolled[1].location
            r2 = rolled[2].location
            a1 = triangle_area_vector(r0, r1, c)
            a2 = triangle_area_vector(r1, r2, c)
            if vertex == self:
                # A_1
                x = .5 * - (r0[2] - c[2]) * a1[1] + .5 * (r0[1] - c[1]) * a1[2]
                y = .5 * (r0[2] - c[2]) * a1[0] + .5 * - (r0[0] - c[0]) * a1[2]
                z = .5 * - (r0[1] - c[1]) * a1[0] + .5 * (r0[0] - c[0]) * a1[1]
                result += np.array([x,y,z])/norm(a1)
                # A_2
                x = .5 * (r2[2] - c[2]) * a2[1] + .5 * - (r2[1] - c[1]) * a2[2]
                y = .5 * - (r2[2] - c[2]) * a2[0] + .5 * (r2[0] - c[0]) * a2[2]
                z = .5 * (r2[1] - c[1]) * a2[0] + .5 * - (r2[0] - c[0]) * a2[1]
                result += np.array([x,y,z])/norm(a2)
            # A_{1,2,j}
            x = .5/n * - (r1[2] - r0[2]) * a1[1] + .5/n * (r1[1] - r0[1]) * a1[2]
            y = .5/n * (r1[2] - r0[2]) * a1[0] + .5/n * - (r1[0] - r0[0]) * a1[2]
            z = .5/n * - (r1[1] - r0[1]) * a1[0] + .5/n * (r1[0] - r0[0]) * a1[1]
            result += np.array([x,y,z])/norm(a1)
        return result
    
    def clockwise_edges(self):
        # Returns the edges connected to a vertex in clockwise order
        # NB: currently not used in simulations
        if len(self.cells) > 2:
            e = [self.edges[0]]
            c = []

            while(len(e) < len(self.edges)):
                for cell in self.cells:
                    if e[-1] in cell.edges:
                        # if previous edge points towards vertex
                        if (cell.edge_direction(e[-1]) == -1 and e[-1].vertices[0] == self) or \
                        ((cell.edge_direction(e[-1]) == 1 and e[-1].vertices[1] == self)):
                            # find edge pointing away from vertex
                            for edge in cell.edges:
                                if edge not in e and self in edge.vertices:
                                    e.append(edge)
                                    c.append(cell)
            for cell in self.cells:
                if cell not in c:
                    c.append(cell)
        else:
            e = self.edges
            c = self.cells
        return (e, c)
    
    def is_apical(self):
        # Checks whether a vertex is on the apical side of the tissue
        return any([e.is_apical() for e in self.edges])
    
    def add_edge(self, edge):
        if edge not in self.edges:
            self.edges.append(edge)
            
    def remove_edge(self, edge):
        if edge in self.edges:
            self.edges.remove(edge)
           
    def add_cell(self, cell):
        if cell not in self.cells:
            self.cells.append(cell)
    
    def remove_cell(self, cell):
        if cell in self.cells:
            self.cells.remove(cell)
     
class Cell:
    
    def __init__(self, name = None, color = 'red', celltype = None, target_volume = None, \
                 apical_tension = None, basal_tension = None,\
                 apical_edges = None, lateral_edges = None, basal_edges = None,\
                 apical_direction = None, basal_direction = None):
        
        if target_volume is None:
            target_volume = p.init_volume
        if basal_tension is None:
            basal_tension = p.basal_tension
        if apical_tension is None:
            apical_tension = p.apical_tension
        if apical_edges is None:
            self.apical_edges = []
        else:
            self.apical_edges = apical_edges
        if lateral_edges is None:
            self.lateral_edges = []
        else:
            self.lateral_edges = lateral_edges
        if basal_edges is None:
            self.basal_edges = []
        else:
            self.basal_edges = basal_edges
        if apical_direction is None:
            self.apical_direction = []
        else:
            self.apical_direction = apical_direction
        if basal_direction is None:
            self.basal_direction = []
        else:
            self.basal_direction = basal_direction
        self.basal_tension = basal_tension
        self.apical_tension = apical_tension
        self.init_vertices = []
        self.name = name
        self.color = color
        self.celltype = celltype
        self.target_volume = target_volume
            
    def center(self, vertices):
        # Returns the center of a surface or the cell itself given a list of vertices
        loc = np.array([0. for dim in vertices[0].location])
        for vertex in vertices:
            loc += vertex.location
        loc = loc/len(vertices)
        return loc
    
    def surface(self, edges):
        # Returns the area of a surface. Input of edges determines which surface
        vertices = self.ordered_vertices(edges)
        c = self.center(vertices)
        locs = [v.location for v in vertices]
        area = 0
        n = len(vertices)
        for i in range(n):
            locs = np.roll(locs, 1, axis = 0)
            area += norm(triangle_area_vector(locs[0], locs[1], c))
        return area
    
    def surface_vector(self, edges):
        # Also returns the surface area, but with a direction
        vertices = self.ordered_vertices(edges)
        locs = [v.location for v in vertices]
        area = 0
        n = len(vertices)
        for i in range(n):
            locs = np.roll(locs, 1, axis = 0)
            area += np.cross(locs[0], locs[1])
        return area/2
    
    def perimeter(self, edges):
        # Returns the perimeter of a surface. Input of edges determines which surface
        length = 0
        for edge in edges:
            length += edge.length()
        return length
    
    def pyramid_volume(self, edges, centroid):
        # Returns the volume of a part of the cell, given two edges and the cell centroid
        vertices = self.ordered_vertices(edges)
        c = self.center(vertices)
        locs = [v.location for v in vertices]
        volume = 0
        n = len(vertices)
        for i in range(n):
            locs = np.roll(locs, 1, axis = 0)
            a = locs[0]; b = locs[1]
            volume += np.linalg.det(np.matrix([[a[0],b[0],c[0],centroid[0]],\
                                               [a[1],b[1],c[1],centroid[1]],\
                                        [a[2],b[2],c[2],centroid[2]],[1,1,1,1]]))
        return volume/6
    
    def volume(self):
        # Returns the volume of the cell
        c = self.center(self.all_vertices())
        volume = 0
        
        # Calculate dot product of every surface with centroid
        volume += self.pyramid_volume(self.apical_edges, c)
        volume += self.pyramid_volume(self.basal_edges, c)
        for lateral in self.lateral_surfaces():
            volume += self.pyramid_volume(lateral, c)
        return volume
    
    def adjacent_cells(self):
        # Returns all cells directly adjacent to this cell
        result = []
        for edge in self.apical_edges:
            for cell in edge.cells:
                if cell is not self:
                    result.append(cell)
        return result
    
    def dist_to_paneth(self, cells):
        # Returns distance between this cell and the closest paneth cell
        result = math.inf
        for c in cells:
            if cells[c].celltype == 'paneth':
                result = min(result, norm(cells[c].center(cells[c].all_vertices()) - self.center(self.all_vertices())))
        return result
    
    def lateral_surfaces(self):
        # Returns a list of lists of edges defining each lateral surface
        lateral_surfaces = []
        # For each lateral edge, find the apical edge that points away from one of the lateral edges' vertices
        for de in self.lateral_edges:
            surface_edges = [de]
            for ae in self.apical_edges:
                if self.directed_edge(ae)[0] in de.vertices:
                    surface_edges.append(ae)
                    # Find the lateral edge connected to the other side of the apical edge
                    for de2 in self.lateral_edges:
                        if self.directed_edge(ae)[1] in de2.vertices:
                            surface_edges.append(de2)
            # Also find the basal edge pointing towards one of the lateral edges' vertices
            for be in self.basal_edges:
                if self.directed_edge(be)[1] in de.vertices:
                    surface_edges.append(be)
            lateral_surfaces.append(surface_edges)
        return lateral_surfaces
    
    def ordered_vertices(self, edges):
        # Returns the vertices on a given surface, in (counter)clockwise order as seen from the outside of the cell
        # NB: whether this is clockwise or counterclockwise depends on initiation of the model
        verts = []
        e0 = edges[0]
        es = [e0]
        # Check if e0 has a direction. If not, use the direction of a connected edge.
        if self.edge_direction(e0):
            verts.append(self.directed_edge(e0)[0])
            verts.append(self.directed_edge(e0)[1])
        else:
            for edge in edges:
                if self.edge_direction(edge):
                    v = self.directed_edge(edge)[0]
                    if v in e0.vertices:
                        # Orientate e0 in such a way that if points towards v
                        verts.append(e0.other_vertex(v))
                        verts.append(v)
        while len(verts) < len(self.vertices(edges)):
            for edge in edges:
                if verts[-1] in edge.vertices and edge not in es and \
                len(verts) < len(self.vertices(edges)):
                    # again: check if edge has a direction
                    if self.edge_direction(edge):
                        verts.append(self.directed_edge(edge)[1])
                    else:
                        for edge2 in edges:
                            if self.edge_direction(edge2):
                                v = self.directed_edge(edge2)[0]
                                if v in edge.vertices:
                                    verts.append(v)
                    es.append(edge)
        # lateral surfaces need to be in the exact opposite direction than the apical and basal edges indicate
        if any([e in self.lateral_edges for e in edges]):
            verts = list(reversed(verts))
        return verts
    
    def vertices(self, edges):
        # Returns the vertices on a given surface in trivial order
        vertices = set()
        for edge in edges:
            vertices.add(edge.vertices[0])
            vertices.add(edge.vertices[1])
        return list(vertices)
      
    def all_vertices(self):
        # Returns all vertices connected to this cell
        all_verts = set()
        for site in [self.basal_edges, self.apical_edges, self.lateral_edges]:
            for edge in site:
                all_verts.add(edge.vertices[0])
                all_verts.add(edge.vertices[1])
        return list(all_verts)    
    
    def edge_direction(self, edge):
        # Given an edge of the cell, return it's direction with respect to this cell
        if edge in self.apical_edges:
            return self.apical_direction[self.apical_edges.index(edge)]
        elif edge in self.basal_edges:
            return self.basal_direction[self.basal_edges.index(edge)]
        else:
            return None
        
    def directed_edge(self, edge):
        # Returns a direction given an edge of this cell
        if self.edge_direction(edge) == 1:
            return [edge.vertices[0],edge.vertices[1]]
        elif self.edge_direction(edge) == -1:
            return [edge.vertices[1],edge.vertices[0]]
        else:
            return None
        
    def ordered_edges(self, es, e0 = None):
        # Returns edges on a surface in (counter)clockwise order, see ordered_vertices
        if not e0:
            e0 = es[0]
        edges = [e0]
        verts = []
        verts.append(self.directed_edge(e0)[0])
        verts.append(self.directed_edge(e0)[1])
        while len(edges) < len(es):
            for edge in es:
                if verts[-1] in edge.vertices and edge not in edges and len(edges) < len(es):
                    verts.append(self.directed_edge(edge)[1])
                    edges.append(edge)
        return edges
    
    def add_edge(self, edge, site):
        # Edge: the edge to add. Edges: the edge list to add it to.
        if site == 'basal':
            edges = self.basal_edges
            direction = self.basal_direction
        if site == 'apical':
            edges = self.apical_edges
            direction = self.apical_direction
        if site == 'lateral':
            self.lateral_edges.append(edge)
        if site != 'lateral':
            direct = None
            for e in edges:
                if e.vertices[0] == edge.vertices[0] \
                or e.vertices[1] == edge.vertices[1]:
                    direct = -1 * self.edge_direction(e)
                elif e.vertices[0] == edge.vertices[1] \
                or e.vertices[1] == edge.vertices[0]:
                    direct = self.edge_direction(e)
            edges.append(edge)
            direction.append(direct)
            
    def remove_edge(self, edge, site):
        if site == 'lateral':
            self.lateral_edges.remove(edge)
        if site == 'apical':
            i = self.apical_edges.index(edge)
            del self.apical_edges[i]
            del self.apical_direction[i]
        if site == 'basal':
            i = self.basal_edges.index(edge)
            del self.basal_edges[i]
            del self.basal_direction[i]

    def random_color(self):
        # Randomly sets a cell's color
        colors = ["crimson","limegreen", "navy", "darkorange", "turquoise","fuchsia"]
        self.color = colors[random.randint(0,len(colors)-1)]
            
class Edge:
    
    def __init__(self, vertices, cells = None, name = None):
        self.vertices = vertices
        if cells is None:
            self.cells = []
        else:
            self.cells = cells
        self.name = name
            
    def length(self):
        # Returns the length of this edge
        return norm(self.vertices[1].location - self.vertices[0].location)
    
    def center(self):
        # Returns the center of this edge
        return (self.vertices[0].location + self.vertices[1].location)/2
    
    def gradient(self, vertex = None):
        # Returns the gradient of a line given a reference vertex
        if vertex is None:
            vertex = self.vertices[0]
        if self.vertices[0] == vertex:
            line = self.vertices[1].location - self.vertices[0].location
        if self.vertices[1] == vertex:
            line = self.vertices[0].location - self.vertices[1].location
        return -1 * line / norm(line)
    
    def other_vertex(self, vertex):
        # Takes a vertex connected to this edge and returns the other
        i = self.vertices.index(vertex)
        return np.roll(self.vertices, 1)[i]
    
    def sister_edge(self):
        # Returns the edge on the other site that is connected to the same lateral edges
        refcell = self.cells[0]
        sides = refcell.lateral_surfaces()
        for s in sides:
            if self in s:
                for edge in s:
                    if self in refcell.apical_edges and edge in refcell.basal_edges:
                        return edge
                    if self in refcell.basal_edges and edge in refcell.apical_edges:
                        return edge
    
    def is_apical(self):
        # Checks whether this edge is on the apical surface
        refcell = self.cells[0]
        return self in refcell.apical_edges
    
    def is_basal(self):
        # Checks whether this edge is on the basal surface
        refcell = self.cells[0]
        return self in refcell.basal_edges
    
    def is_lateral(self):
        # Checks whether this edge is on the lateral surface
        # NB: this is formulated slightly different due to some errors during
        # transitions, but it functions the same as the other two
        return any([self in cell.lateral_edges for cell in self.cells])
    
    def add_vertex(self, vertex):
        if vertex not in self.vertices:
            self.vertices.append(vertex)
    
    def remove_vertex(self, vertex):
        if vertex in self.vertices:
            self.vertices.remove(vertex)
    
    def add_cell(self, cell):
        if cell not in  self.cells:
            self.cells.append(cell)
    
    def remove_cell(self, cell):
        if cell in self.cells:
            self.cells.remove(cell)
            
### Functions used in the classes ###

def triangle_area_vector(a, b, c):
    return .5 * np.cross(b - c, a - c)

def gs(v1, v2):
    # Gramm-Schmidt for 2 vectors
    result = v2 - np.dot(v1, v2) / np.dot(v1, v1) * v1
    return result/norm(result)