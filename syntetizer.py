import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import math
import random

from geometrypack.algorithms import normalize
from geometrypack.parametric_functions import cylinder
from geometrypack.dataio import read_OFF, write_OFF, write_PLY

data_folder = os.path.join("data", "meshes")
results_folder = os.path.join("results", "normal_interpolation")

FLIP_NONE = "regularflip"
FLIP_ALTERNATE = "alternateflip"
FLIP_RANDOM = "randomflip"

meridians = 30
parallels = 15
height = 4

def sample_cylinder_points(meridians, parallels, height=2, noise = False):
    inc_height = height/parallels
    inc_rad = 2*math.pi/meridians
    vertices = []
    h = 0
    while(h <= height):
        angle = 0
        for m in range(meridians):
            if noise:
                vertices.append(cylinder( angle + random.gauss(0, inc_rad/7), h + random.uniform(h-inc_height/3, h+inc_height/3) ))
            else:
                vertices.append(cylinder( angle, h))
            angle += inc_rad
        h+= inc_height

    return vertices

def cylinder_analytic_normal(vertex):
    return normalize(np.array([2*vertex[0], 2*vertex[1], 0]))

# def write_regular_cylinder_mesh(noise = False):
#     height = 2
#
#     vertices = sample_cylinder_points(meridians, parallels, height, noise)
#     indices = []
#     for p in range(parallels):
#         for m in range(meridians):
#             indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
#             indices.append([p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m ]) #inferior
#
#     indices = np.array(indices)
#     if noise:
#         write_OFF(os.path.join(data_folder,
#             "noisecylinder{}_{}_regular.off".format(meridians, parallels)), vertices, indices)
#     else:
#         write_OFF(os.path.join(data_folder,
#             "cylinder{}_{}_regular.off".format(meridians, parallels)), vertices, indices)
#     print("OFF escrito")
#
# def write_alternate_flip_cylinder_mesh(noise = False):
#     height = 2
#
#     vertices = sample_cylinder_points(meridians, parallels, height, noise)
#     indices = []
#     for p in range(parallels):
#         for m in range(meridians):
#             if p%2 == 0:
#                 indices.append([p*meridians + m, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
#                 indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians ]) #inferior
#             else:
#                 indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
#                 indices.append([p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m ]) #inferior
#
#     indices = np.array(indices)
#
#     if noise:
#         write_OFF(os.path.join(data_folder,
#             "noisecylinder{}_{}_alternateflip.off".format(meridians, parallels)), vertices, indices)
#     else:
#         write_OFF(os.path.join(data_folder,
#             "cylinder{}_{}_alternateflip.off".format(meridians, parallels)), vertices, indices)
#     print("OFF escrito")
#
# def write_random_flip_cylinder_mesh(noise = False):
#     height = 2
#
#     vertices = sample_cylinder_points(meridians, parallels, height, noise)
#     indices = []
#     for p in range(parallels):
#         for m in range(meridians):
#             k = random.randint(0, 1)
#             if k%2 == 0:
#                 indices.append([p*meridians + m, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
#                 indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians ]) #inferior
#             else:
#                 indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
#                 indices.append([p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m ]) #inferior
#
#     indices = np.array(indices)
#
#     if noise:
#         write_OFF(os.path.join(data_folder,
#             "noisecylinder{}_{}_randomflip.off".format(meridians, parallels)), vertices, indices)
#     else:
#         write_OFF(os.path.join(data_folder,
#             "cylinder{}_{}_randomflip.off".format(meridians, parallels)), vertices, indices)
#     print("OFF escrito")

def write_cylinder_mesh(flip_mode, noise = False, mer = meridians, par = parallels, hei = height):
    vertices = sample_cylinder_points(mer, par, hei, noise)

    k = 0
    indices = []
    for p in range(parallels):
        for m in range(meridians):
            if flip_mode == FLIP_RANDOM:
                k = random.randint(0, 1)
            elif flip_mode == FLIP_ALTERNATE:
                k = p
            if k%2 == 0:
                indices.append([p*meridians + m, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
                indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians ]) #inferior
            else:
                indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
                indices.append([p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m ]) #inferior

    indices = np.array(indices)

    if noise:
        write_OFF(os.path.join(data_folder,
            "noisecylinder{}_{}_{}.off".format(meridians, parallels)), vertices, indices)
    else:
        write_OFF(os.path.join(data_folder,
            "cylinder{}_{}_{}.off".format(meridians, parallels)), vertices, indices)
    print("Cylinder {} written".format(flip_mode))


def create_data():
    modes = [FLIP_NONE, FLIP_ALTERNATE, FLIP_RANDOM]
    for mode in modes:
        write_cylinder_mesh(mode)
        write_cylinder_mesh(mode, True)
    # write_alternate_flip_cylinder_mesh()
    # write_random_flip_cylinder_mesh()
    #
    # write_regular_cylinder_mesh(True)
    # write_alternate_flip_cylinder_mesh(True)
    # write_random_flip_cylinder_mesh(True)
