import argparse
import numpy as np
import os

import matplotlib.pyplot as plt
from math import pi, sqrt
from geometrypack.dataio import read_OFF, write_PLY, write_OFF, write_uv_PLY
from geometrypack.data_structures import MinPriorityQueue
from geometrypack.algorithms import distance_vertex_plane, normalize
from geometrypack import drawing
from geometrypack.meshes import compute_neighborhood, compute_neighborhood_triangles, compute_boundary_indices, rotate
from normal_interpolation import normals_mean_interpolation, normals_barycentric_interpolation
from copy import deepcopy


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

def triangles_for_shared_edge(indices, edge):
    remaining_vertices, tindex = [], []
    num_triangles = len(indices)
    for i in range(num_triangles):
        t = set(indices[i])
        if len(t.intersection(edge)) == 2:
            remaining_vertices.extend(list(t.difference(edge)))
            tindex.append(i)
    return remaining_vertices, tindex

def split_long_edges(vertices, indices, high):
    # make a copy of vertices and indices as a simple python list
    # because we need to append multiple times
    vertices = [v for v in vertices]
    indices = [i for i in indices]
    edges = compute_edges(indices)
    num_triangles = len(indices)
    print("num edges:", len(edges))
    print("high bar:", high)
    for edge in edges.values():
        length = np.linalg.norm(vertices[edge[0]] - vertices[edge[1]])
        if length > high:
            # print("high", edge, length)
            remaining_vertices, tindex = triangles_for_shared_edge(indices, edge)
            middle = (vertices[edge[0]] + vertices[edge[1]])/2
            middle_index = len(vertices)
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
        else:
            # print("not high", edge, length)
            pass
    return np.array(vertices), np.array(indices)

def removearray(L,arr):
    ind = 0
    size = len(L)
    while ind != size and not np.array_equal(L[ind],arr):
        ind += 1
    if ind != size:
        L.pop(ind)
    else:
        raise ValueError('array not found in list.')

def collapse_short_edges(vertices, indices, low, high):
    vertices = [v for v in vertices]
    indices = [i for i in indices]
    edges = compute_edges(indices)
    vertices_neighborhood = compute_neighborhood(vertices, indices)
    triangles_neighborhood = compute_neighborhood_triangles(vertices, indices)
    triangles_to_be_removed = []
    for edge in edges.values():
        length = np.linalg.norm(edge[0] - edge[1])
        if length < low:
            should_collapse = True
            middle = (vertices[edge[0]] + vertices[edge[1]])/2
            for v in vertices_neighborhood[edge[0]]:
                if np.linalg.norm(vertices[v] - middle) > high:
                    should_collapse = False
                    break
            if should_collapse:
                for tindex in triangles_neighborhood[edge[0]]:
                    triangle = indices[tindex]
                    if edge[1] in triangle:
                        # print(edge[1], triangle)
                        triangles_to_be_removed.append(tindex)
                        continue
                    for j in range(len(triangle)):
                        if triangle[j] == edge[0]:
                            triangle[j] = edge[1]
                vertices[edge[1]] = middle
    # triangles_to_be_removed.sort()
    # triangles_to_be_removed.reverse()
    # print(triangles_to_be_removed)
    # for tindex in triangles_to_be_removed:
    #     indices.__delitem__(tindex)

    # print(indices)
    return np.array(vertices), np.array(indices)

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
    neighborhood = compute_neighborhood(vertices, indices)
    #normals = normals_mean_interpolation(vertices, indices, neighborhood)
    normals = normals_barycentric_interpolation(vertices, indices, neighborhood)
    barycentric_pos = []
    for i in range(len(vertices)):
        neighbors = neighborhood[i]
        q = np.array([0, 0, 0])
        for index in neighbors:
            q = q + vertices[index]
        q = q/len(neighbors)
        barycentric_pos.append(q)
    for i in range(len(vertices)):
        q = barycentric_pos[i]
        vertices[i] = q + np.dot(normals[i], vertices[i] - q)*normals[i]

def project_to_surface(remeshed_vertices, remeshed_indices, vertices, indices):
    pass

#polygon mesh processing incremental remeshing
def remesh(vertices, indices, target_length, maxiterations=10):
    # neighborhood = compute_neighborhood(vertices, indices)
    edges = compute_edges(indices)
    edge_lengths = []
    for edge in edges.values():
        edge_lengths.append(np.linalg.norm(edge[0] - edge[1]))
    target_length = np.mean(np.array(edge_lengths))*0.1
    print("target_length: ", target_length)
    low = lowcoefficient*target_length
    high = highcoefficient*target_length
    remeshed_vertices = deepcopy(vertices)
    remeshed_indices = deepcopy(indices)
    for i in range(maxiterations):
        print("splitting long edges")
        remeshed_vertices, remeshed_indices = split_long_edges(remeshed_vertices, remeshed_indices, high)
        # print("collapsing short edges")
        # remeshed_vertices, remeshed_indices = collapse_short_edges(remeshed_vertices, remeshed_indices, low, high)
        print("equalizing valences")
        equalize_valences(remeshed_indices)

        # print("tangential relaxation...")
        # tangential_relaxation(remeshed_vertices, remeshed_indices)
    # project_to_surface(remeshed_vertices, remeshed_indices, vertices, indices)

    print(len(remeshed_vertices), len(remeshed_indices))
    return remeshed_vertices, remeshed_indices

def main():
    vertices, indices = read_OFF("data/meshes/cube.off")
    vertices, indices = remesh(vertices, indices, 0.000001)
    print(len(vertices), len(indices))
    write_OFF("teste.off", vertices, indices)

if __name__ == '__main__':
    main()
