import argparse
import numpy as np
import os

import matplotlib.pyplot as plt
from math import pi, sqrt
from geometrypack.dataio import read_OFF, write_PLY, write_OFF, write_uv_PLY
from geometrypack.algorithms import distance_vertex_plane, normalize
from geometrypack import drawing
from geometrypack.meshes import compute_neighborhood, compute_boundary_indices, rotate


lowcoefficient = 4/5
highcoefficient = 4/3

def split_long_edges(high):
    pass
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

def flip(edge):
    pass

def valence(vindex):
    pass

def equalize_valences():
    pass
    # for edge in edges:
    #     le a, b, c, d be the vertices of the two triangles adjacent to edge
    #     deviation_pre = (abs(valence(a) - targetval(a))
    #                     + abs(valence(b) - targetval(b))
    #                     + abs(valence(c) - targetval(c))
    #                     + abs(valence(d) - targetval(d)))
    #     flip(edge)
    #     deviation_post = (abs(valence(a) - targetval(a))
    #                     + abs(valence(b) - targetval(b))
    #                     + abs(valence(c) - targetval(c))
    #                     + abs(valence(d) - targetval(d)))
    #     if deviation_pre <= deviation_post:
    #         flip(e)

def tangential_relaxation(vindex, normal, neighbors):
    q = sum(neighbors)
    p = q + np.dot(normal - q)*normal
    return p

def project_to_surface():
    pass

#polygon mesh processing incremental remeshing
def remesh(target_length, maxiterations=10):
    low = lowcoefficient*target_length
    high = highcoefficient*target_length
    for i in range(maxiterations):
        split_long_edges(high)
        collapse_short_edges(low, high)
        equalize_valences()
        tangential_relaxation()
        project_to_surface()
