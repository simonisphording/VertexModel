# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 11:58:55 2019

@author: isphording
"""

import numpy as np
from numpy.linalg import norm
import parameters as p
import math
import random

### Classes describing vertices and cells ###

class Vertex:
    
    def __init__(self, location, edges = None, cells = None, bound = False, name = None):
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
    
    def force(self):
        # Calculates the force on a vertex
        surface = [0 for dimension in self.location]
        perimeter = [0 for dimension in self.location]
        curvature = [0 for dimension in self.location]
        line = [0 for dimension in self.location]
        for cell in self.cells:
            surface += -1 * p.elastic * (norm(cell.surface()) - cell.target_volume) * self.surface_gradient(cell)
            perimeter += -1 * p.contract * cell.perimeter() * self.perimeter_gradient(cell)
            if cell.celltype == 'paneth':
                curvature += -1 * p.gradient * self.curvature_gradient(cell)
            #else:
            #    curvature += p.gradient * self.curvature_gradient(cell)
        for edge in self.edges:
            line += -1 * p.tension * edge.gradient(vertex = self)
        return surface + line + perimeter + curvature
    
    def surface_gradient(self, cell):
        # Returns the gradient of the surface of a cell
        vertices = cell.ordered_vertices()
        result = [0 for dim in self.location]
        c = cell.center()
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
    
    def perimeter_gradient(self, cell):
        # Returns gradient of the perimeter of a cell
        result = [0 for dimension in self.location]
        for edge in self.edges:
            if edge in cell.edges:
                result += edge.gradient(vertex = self)
        return result
    
    def surface_vector_gradient(self, cell):
        vertices = cell.ordered_vertices()
        result = [0 for dim in self.location]
        i = vertices.index(self)
        rolled = np.roll(vertices, 1-i)
        r0, r1, r2 = [v.location for v in rolled[:3]]
        result[0] = r2[1] - r2[2] + r0[2] - r0[1]
        result[1] = r2[2] - r2[0] + r0[0] - r0[2]
        result[2] = r2[0] - r2[1] + r0[1] - r0[0]
        return result
    
    def clockwise_edges(self):
        # Returns all edges connected to this vertex in clockwise order
        # Currently not used
        if len(self.cells) > 2:
            e = [self.edges[0]]
            c = []

            while(len(e) < len(self.edges)):
                for cell in self.cells:
                    if e[-1] in cell.edges:
                        # if previous edge points away from vertex
                        if (cell.edge_direction(e[-1]) == 1 and e[-1].vertices[0] == self) or \
                        ((cell.edge_direction(e[-1]) == -1 and e[-1].vertices[1] == self)):
                            # find edge pointing towards vertex
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
    
    def __init__(self, edges = None, direction = None, target_volume = p.init_surface, name = None, color = 'red', celltype = None):
        if edges is None:
            self.edges = []
        else:
            self.edges = edges
        if direction is None:
            self.direction = []
        else:
            self.direction = direction
        self.name = name
        self.color = color
        self.celltype = celltype
        self.target_volume = target_volume
            
    def center(self):
        # Returns the center of the cell
        loc = np.array([0. for dim in self.vertices()[0].location])
        for vertex in self.vertices():
            loc += vertex.location
        loc = loc/len(self.vertices())
        return loc
    
    def surface_vector(self):
        # Returns a surface vector
        vertices = self.ordered_vertices()
        locs = [v.location for v in vertices]
        area = np.array([0. for dim in vertices[0].location])
        n = len(vertices)
        for i in range(n):
            locs = np.roll(locs, 1, axis = 0)
            area += np.cross(locs[0], locs[1])
        return area/2
    
    def surface(self):
        # Returns the surface area of this cell
        vertices = self.ordered_vertices()
        c = self.center()
        locs = [v.location for v in vertices]
        area = 0
        n = len(vertices)
        for i in range(n):
            locs = np.roll(locs, 1, axis = 0)
            area += norm(triangle_area_vector(locs[0], locs[1], c))
        return area
    
    def perimeter(self):
        # Returns the perimeter of this cell
        length = 0
        for edge in self.edges:
            length += edge.length()
        return length
    
    def adjacent_cells(self):
        # Returns a list of adjacent cells
        result = []
        for edge in self.edges:
            for cell in edge.cells:
                if cell is not self:
                    result.append(cell)
        return result
    
    def ordered_vertices(self):
        # Returns the vertices of this cell, in clockwise order
        e0 = self.edges[0]
        edges = [e0]
        verts = []
        verts.append(self.directed_edge(e0)[0])
        verts.append(self.directed_edge(e0)[1])
        while len(verts) < len(self.vertices()):
            for edge in self.edges:
                if verts[-1] in edge.vertices and edge not in edges and len(verts) < len(self.vertices()):
                    verts.append(self.directed_edge(edge)[1])
                    edges.append(edge)
        return verts
    
    def vertices(self):
        # Returns the vertices of this cell in a random order
        vertices = set()
        for edge in self.edges:
            vertices.add(edge.vertices[0])
            vertices.add(edge.vertices[1])
        return list(vertices)
        
    def edge_direction(self, edge):
        # Returns the direction of an edge with respect to this cell
        return self.direction[self.edges.index(edge)]
    
    def directed_edge(self, edge):
        # Returns an edge with the correct direction for this cell
        if self.edge_direction(edge) == 1:
            return [edge.vertices[0],edge.vertices[1]]
        else:
            return [edge.vertices[1],edge.vertices[0]]
    
    def ordered_edges(self, e0 = None):
        # Returns the edges of this cell, in clockwise order
        if not e0:
            e0 = self.edges[0]
        edges = [e0]
        verts = []
        verts.append(self.directed_edge(e0)[0])
        verts.append(self.directed_edge(e0)[1])
        while len(edges) < len(self.edges):
            for edge in self.edges:
                if verts[-1] in edge.vertices and edge not in edges and len(edges) < len(self.edges):
                    verts.append(self.directed_edge(edge)[1])
                    edges.append(edge)
        return edges
    
    def add_edge(self, edge):
        direction = None
        for e in self.edges:
            if e.vertices[0] == edge.vertices[0] \
            or e.vertices[1] == edge.vertices[1]:
                direction = -1 * self.edge_direction(e)
            elif e.vertices[0] == edge.vertices[1] \
            or e.vertices[1] == edge.vertices[0]:
                direction = self.edge_direction(e)
        if direction is None:
            print('Warning: edge does not connect to any other edge in the cell')
        self.edges.append(edge)
        self.direction.append(direction)
            
    def remove_edge(self, edge):
        i = self.edges.index(edge)
        del self.edges[i]
        del self.direction[i]
        
    def random_color(self):
        # Randomly assigns a color to this cell used in plotting
        colors = ["crimson","limegreen", "navy", "darkorange", "turquoise","fuchsia"]
        self.color = colors[random.randint(0,len(colors)-1)]
    
    def dist_to_paneth(self, cells):
        # Returns the euclidean distance to the closest Paneth cell
        result = math.inf
        for c in cells:
            if cells[c].celltype == 'paneth':
                result = min(result, norm(cells[c].center() - self.center()))
        return result

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
        # Returns the gradient along this edge with respect to one of its vertices
        if vertex is None:
            vertex = self.vertices[0]
        if self.vertices[0] == vertex:
            line = self.vertices[1].location - self.vertices[0].location
        if self.vertices[1] == vertex:
            line = self.vertices[0].location - self.vertices[1].location
        return -1 * line / norm(line)
    
    def other_vertex(self, vertex):
        # Takes a vertex in this edge and returns the other
        i = self.vertices.index(vertex)
        return np.roll(self.vertices, 1)[i]
    
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
    # Gramm-Schmidt
    result = v2 - np.dot(v1, v2) / np.dot(v1, v1) * v1
    return result/norm(result)