import argparse
import numpy as np
import os

import matplotlib.pyplot as plt
from math import pi, sqrt
from geometrypack.dataio import read_OFF, write_PLY, write_OFF, write_uv_PLY
from geometrypack.data_structures import MinPriorityQueue
from geometrypack.algorithms import distance_vertex_plane, normalize
from geometrypack import drawing
from geometrypack.meshes import compute_neighborhood, compute_boundary_indices, rotate


lowcoefficient = 4/5
highcoefficient = 4/3


def compute_edges(indices):
    edges = {}
    # separates all edges using a dictionary to avoid duplicates
    for triangle in indices:
        for i in range(3):
            edge = (triangle[i], triangle[(i+1)%3])
            edgehash = str(min(edge)) + str(max(edge))
            if edgehash not in edges:
                edges[edgehash] = edge
    return edges

def split_long_edges(vertices, indices, high):
    # make a copy of vertices and indices as a simple python list
    # because we need to append multiple times
    vertices = [v for v in vertices]
    indices = [i for i in indices]
    edges = compute_edges(indices)
    num_triangles = len(indices)
    for edge in edges.values():
        length = np.linalg.norm(edge[0] - edge[1])
        if length > high:
            remaining_vertices = []
            tindex = []
            for i in range(num_triangles):
                t = set(indices[i])
                if len(t.intersection(edge)) == 2:
                    remaining_vertices.extend(list(t.difference(edge)))
                    tindex.append(i)
            middle_index = len(vertices)
            middle = (vertices[edge[0]] + vertices[edge[1]])/2
            vertices.append(middle)
            for j in range(len(tindex)):
                oldt = indices[tindex[j]]
                newt1, newt2 = [], []
                for index in oldt:
                    if index == edge[1]:
                        newt1.append(middle_index)
                        newt2.append(index)
                    elif index == edge[0]:
                        newt1.append(index)
                        newt2.append(middle_index)
                    else:
                        newt1.append(remaining_vertices[j])
                        newt2.append(remaining_vertices[j])
                indices[tindex[j]] = newt1
                indices.append(newt2)
    return np.array(vertices), np.array(indices)

def collapse_short_edges(vertices, indices, low, high):
    pass
    # while exists edge e with length(e) < low:
    #     let e = (a, b) and let a[1]...a[n] be the one-ring of a
    #     should_collapse = true
    #     for i in range(n):
    #         if length(b, a[i]) > high:
    #             should_collapse = false
    #     if should_collapse:
    #         collapse a into b along e

def flip(indices, edge, tindex = None, remaining_vertices = None):
    if triangles_indices is None or remaining_vertices is None:
        remaining_vertices = []
        tindex = []
        for i in range(num_triangles):
            t = set(indices[i])
            if len(t.intersection(edge)) == 2:
                remaining_vertices.extend(list(t.difference(edge)))
                tindex.append(i)

    for j in range(len(tindex)):
        t = indices[tindex[j]]
        for i in range(len(t)):
            if t[i] == edge[1]:
                t[i] = remaining_vertices[j]


# number of neighbors a vertex vindex has
def valence(vindex, indices):
    return len([t for t in indices if vindex in t])

def equalize_valences(indices, targetval=6):
    edges = compute_edges(indices)
    num_triangles = len(indices)
    for edge in edges:
        remaining_vertices = []
        tindex = []
        for i in range(num_triangles):
            t = set(indices[i])
            if len(t.intersection(edge)) == 2:
                remaining_vertices.extend(list(t.difference(edge)))
                tindex.append(i)
        if len(tindex) == 2:
        # let a, b, c, d be the vertices of the two triangles adjacent to edge
            a, b, c, d = edge[0], edge[1], remaining_vertices[0], remaining_vertices[1]
            deviation_pre = (abs(valence(a) - targetval(a))
                            + abs(valence(b) - targetval(b))
                            + abs(valence(c) - targetval(c))
                            + abs(valence(d) - targetval(d)))
            flip(indices, edge, tindex, remaining_vertices)
            deviation_post = (abs(valence(a) - targetval(a))
                            + abs(valence(b) - targetval(b))
                            + abs(valence(c) - targetval(c))
                            + abs(valence(d) - targetval(d)))
            if deviation_pre <= deviation_post:
                flip(indices, (remaining_vertices[0], remaining_vertices[1]), edge)

def tangential_relaxation(vertices, indices):
    # normals = compute_normals(vertices, indices)
    # TODO: interpolate normals
    neighborhood = compute_neighborhood(vertices, indices)
    for i in range(len(vertices)):
        neighbors = neighborhood[i]
        q = sum(neighbors)
        p = q + np.dot(normals[i] - q)*normal
        vertices[i] = p

def project_to_surface(remeshed_vertices, remeshed_indices, vertices, indices):
    pass

#polygon mesh processing incremental remeshing
def remesh(vertices, indices, target_length, maxiterations=1):
    neighborhood = compute_neighborhood(vertices, indices)
    low = lowcoefficient*target_length
    high = highcoefficient*target_length
    for i in range(maxiterations):
        remeshed_vertices, remeshed_indices = split_long_edges(vertices, indices, high)
        # collapse_short_edges(vertices, indices, low, high)
        equalize_valences(remeshed_indices)
        # tangential_relaxation(remeshed_vertices, remeshed_indices)
        # project_to_surface(remeshed_vertices, remeshed_indices, vertices, indices)

    print(len(vertices), len(indices))
    return vertices, indices

def main():
    vertices, indices = read_OFF("data/meshes/manequin.off")
    vertices, indices = remesh(vertices, indices, 0.000001)
    print(len(vertices), len(indices))
    write_OFF("teste.off", vertices, indices)

if __name__ == '__main__':
    main()
