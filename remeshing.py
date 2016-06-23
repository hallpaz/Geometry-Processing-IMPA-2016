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

def split_long_edges(vertices, indices, high):
    queue = MinPriorityQueue()
    for triangle in indices:
        for i in range(3):
            edge = (triangle[i], triangle[(i+1)%3])
            queue.add_point(edge_info, np.linalg.norm(edge[0] - edge[1]))
    candidate = queue.pop_point()
    while(candidate[0] > high):
        edge = candidate[-1]
        middle = (edge[0] + edge[1])/2
        indices_triangles = [i_t for i_t in enumerate(indices) if set(edge).intersection(i_t[1]) == 2]
        # indices[triangles[0][0]] = [edge[0], middle, set(triangles[1]).difference(edge)[0]]
        # indices[triangles[1][0]] = [edge[0], middle, set(triangles[1]).difference(edge)[1]]
        # indices.append([edge[0], middle, set(triangles[1]).difference(edge)[0]])
        # indices.append([edge[0], middle, set(triangles[1]).difference(edge)[0]])

    # I will need a priority queue here
    # while exists edge e with length(e) > high:
    #     split e at midpoint

def collapse_short_edges(low, high):
    pass
    # while exists edge e with length(e) < low:
    #     let e = (a, b) and let a[1]...a[n] be the one-ring of a
    #     should_collapse = true
    #     for i in range(n):
    #         if length(b, a[i]) > high:
    #             should_collapse = false
    #     if should_collapse:
    #         collapse a into b along e

def flip(indices, edge):
    pass

# number of neighbors a vertex vindex has
def valence(vindex, indices):
    return len([t for t in indices if vindex in t])

def equalize_valences(indices, targetval=6):
    # for edge in edges:
    #     le a, b, c, d be the vertices of the two triangles adjacent to edge
    #     deviation_pre = (abs(valence(a) - targetval(a))
    #                     + abs(valence(b) - targetval(b))
    #                     + abs(valence(c) - targetval(c))
    #                     + abs(valence(d) - targetval(d)))
    #     flip(indices, edge)
    #     deviation_post = (abs(valence(a) - targetval(a))
    #                     + abs(valence(b) - targetval(b))
    #                     + abs(valence(c) - targetval(c))
    #                     + abs(valence(d) - targetval(d)))
    #     if deviation_pre <= deviation_post:
    #         flip(indices, e)

def tangential_relaxation(vindex, normal, neighbors):
    q = sum(neighbors)
    p = q + np.dot(normal - q)*normal
    return p

def project_to_surface():
    pass

#polygon mesh processing incremental remeshing
def remesh(vertices, indices, target_length, maxiterations=10):
    neighborhood = compute_neighborhood(vertices, indices)
    low = lowcoefficient*target_length
    high = highcoefficient*target_length
    for i in range(maxiterations):
        split_long_edges(vertices, indices, high)
        collapse_short_edges(vertices, indices, low, high)
        equalize_valences(indices)
        tangential_relaxation()
        project_to_surface()
