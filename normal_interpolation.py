import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import math
import random

from geometrypack.algorithms import normalize
from geometrypack.parametric_functions import cylinder
from geometrypack.dataio import read_OFF, write_OFF, write_PLY
from geometrypack.meshes import compute_neighborhood, compute_neighborhood_triangles
from geometrypack.drawing import colorize


data_folder = os.path.join("data", "meshes")
results_folder = os.path.join("results", "normal_interpolation")


def compute_normal(vertices, indices, tindex):
    t = indices[tindex]
    return triangle_normal([vertices[t[0]], vertices[t[1]], vertices[t[2]]])

# tvertices is a list of 3 numpy array
def triangle_normal(tvertices):
    normal = np.cross(tvertices[1] - tvertices[0], tvertices[2] - tvertices[1])
    return normalize(normal)

def compute_area(vertices, tindex):
    tvertices = vertices[tindex]
    return compute_area(tvertices)

def triangle_area(tvertices):
    a = np.linalg.norm(tvertices[0] - tvertices[1])
    b = np.linalg.norm(tvertices[1] - tvertices[2])
    c = np.linalg.norm(tvertices[2] - tvertices[0])
    # calculate the semi-perimeter
    s = (a + b + c) / 2
    # calculate the area
    area = (s*(s-a)*(s-b)*(s-c)) ** 0.5
    return area

def compute_angle(vertices, indices, triangle_index, pivot_index):
    local_index = 0
    t = indices[triangle_index]
    for i in t:
        if i == pivot_index:
            break
        local_index += 1
    # the indices of the other vertices on this triangle
    a, b = (local_index + 1)%3, (local_index + 2)%3
    pivot, v1, v2 = vertices[local_index], vertices[a], vertices[b]
    u, v = pivot - v1, pivot - v2

    c = np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v) # -> cosine of the angle
    angle = np.arccos(np.clip(c, -1, 1)) # if you really want the angle
    return angle


def vertex_mean_normal_interpolation(faces_normals, triangles_neighbors, vertex_index):
    normal = np.array([0, 0, 0])
    for triangle_index in triangles_neighbors:
        normal = normal + faces_normals[triangle_index]
    if len(neighbors) == 0:
        print("null neighborhood", vertex_index)
        return normal
    return normalize(normal)

def area_interpolation(vertices, indices, vneighboorhood, vertex_index):
    normal = np.array([0, 0, 0])
    for triangle_index in vneighborhood:
        tvertices = [vertices[i] for i in indices[triangle_index]]
        normal = normal + compute_area(tvertices) * triangle_normal(tvertices)
    normal = normalize(normal)
    return normal

def compute_barycentric_area(vertices, indices, triangle_index, pivot_index):
    t = indices[triangle_index]
    #index of the pivot on this triangle
    local_index = 0
    for i in t:
        if i == pivot_index:
            break
        local_index += 1
    # the indices of the other vertices on this triangle
    a, b = (local_index + 1)%3, (local_index + 2)%3
    pivot, v1, v2 = vertices[local_index], vertices[a], vertices[b]
    barycenter = (pivot + v1 + v2)/3
    area = 0
    #first triangle area:
    a = np.linalg.norm((pivot - v1)/2)
    b = np.linalg.norm(pivot - barycenter)
    c = np.linalg.norm((v1 + (pivot-v1)/2) - barycenter)
    s = (a + b + c) / 2
    # increment the area
    area += (s*(s-a)*(s-b)*(s-c)) ** 0.5

    #second triangle area:
    a = np.linalg.norm((pivot - v2)/2)
    b = np.linalg.norm(pivot - barycenter)
    c = np.linalg.norm((v2 + (pivot-v2)/2) - barycenter)
    s = (a + b + c) / 2
    # calculate the area
    area += (s*(s-a)*(s-b)*(s-c)) ** 0.5
    return area

def barycentric_area_interpolation(vertices, indices, vneighboorhood, vertex_index):
    normal = np.array([0, 0, 0])
    for triangle_index in vneighborhood:
        tvertices = [vertices[i] for i in indices[triangle_index]]
        normal = normal + compute_barycentric_area(triangle_index, vertex_index) * triangle_normal(triangle_index)
    normal = normalize(normal)
    return normal

def vertex_barycentric_normal_interpolation(vertices, indices, faces_normals, triangles_neighbors, vertex_index):
    normal = np.array([0, 0, 0])
    for triangle_index in triangles_neighbors:
        normal = normal + compute_barycentric_area(vertices, indices, triangle_index, vertex_index) * faces_normals[triangle_index]
    if len(triangles_neighbors) == 0:
        print("null neighborhood", vertex_index)
        return normal
    return normalize(normal)

def angle_interpolation(vertex_index):
    normal = np.array([0, 0, 0])
    for triangle_index in neighborhood[vertex_index]:
        normal = normal + compute_angle(triangle_index, vertex_index) * triangle_normal(triangle_index)
    normal = normalize(normal)
    return normal


def normals_mean_interpolation(vertices, indices, neighborhood = None):
    if neighborhood is None:
        neighborhood = compute_neighborhood_triangles(vertices, indices)
    faces_normals = [compute_normal(vertices, indices, tindex) for tindex in range(len(indices))]
    vertices_normals = [vertex_mean_normal_interpolation(faces_normals, neighborhood[i], i)
                        for i in range(len(vertices))]
    return vertices_normals

def normals_barycentric_interpolation(vertices, indices, neighborhood = None):
    if neighborhood is None:
        neighborhood = compute_neighborhood_triangles(vertices, indices)
    faces_normals = [compute_normal(vertices, indices, tindex) for tindex in range(len(indices))]
    vertices_normals = [vertex_barycentric_normal_interpolation(vertices, indices, faces_normals, neighborhood[i], i)
                        for i in range(len(vertices))]
    return vertices_normals

# calcula a normal para cada triângulo
# para cada vértice, identifique todos os triângulos que possuem esse vértice
# faça uma média ponderara dos valores da normal de cada triângulo que contém
# o vértice para enontrar a normal do vertice
# def interpolate_normals(filepath, methods, compute_analytic_normal, should_colorize = False):
#     vertices, indices = read_OFF(filepath)
#     filename = os.path.splitext(os.path.basename(filepath))[0]
#     print("BEGIN", filepath)
#     neighborhood = compute_neighborhood(vertices, indices)
#     print("computed neighborhood")
#
#     # triangle is a list of 3 numpy array
#     # def compute_normal(triangle_index):
#     #     t = indices[triangle_index]
#     #     a, b, c = t[0], t[1], t[2]
#     #
#     #     normal = np.cross(vertices[b] - vertices[a], vertices[c] - vertices[b])
#     #     return normalize(normal)
#
#     # def compute_area(triangle_index):
#     #     t = indices[triangle_index]
#     #     a = np.linalg.norm(vertices[t[0]] - vertices[t[1]])
#     #     b = np.linalg.norm(vertices[t[1]] - vertices[t[2]])
#     #     c = np.linalg.norm(vertices[t[2]] - vertices[t[0]])
#     #     # calculate the semi-perimeter
#     #     s = (a + b + c) / 2
#     #     # calculate the area
#     #     area = (s*(s-a)*(s-b)*(s-c)) ** 0.5
#     #     return area
#
#     def compute_barycentric_area(triangle_index, pivot_index):
#         t = indices[triangle_index]
#         #index of the pivot on this triangle
#         local_index = 0
#         for i in t:
#             if i == pivot_index:
#                 break
#             local_index += 1
#         # the indices of the other vertices on this triangle
#         a, b = (local_index + 1)%3, (local_index + 2)%3
#         pivot, v1, v2 = vertices[local_index], vertices[a], vertices[b]
#         barycenter = (pivot + v1 + v2)/3
#         area = 0
#         #first triangle area:
#         a = np.linalg.norm((pivot - v1)/2)
#         b = np.linalg.norm(pivot - barycenter)
#         c = np.linalg.norm((v1 + (pivot-v1)/2) - barycenter)
#         s = (a + b + c) / 2
#         # increment the area
#         area += (s*(s-a)*(s-b)*(s-c)) ** 0.5
#
#         #second triangle area:
#         a = np.linalg.norm((pivot - v2)/2)
#         b = np.linalg.norm(pivot - barycenter)
#         c = np.linalg.norm((v2 + (pivot-v2)/2) - barycenter)
#         s = (a + b + c) / 2
#         # calculate the area
#         area += (s*(s-a)*(s-b)*(s-c)) ** 0.5
#
#         return area
#
#     def compute_angle(triangle_index, pivot_index):
#         local_index = 0
#         t = indices[triangle_index]
#         for i in t:
#             if i == pivot_index:
#                 break
#             local_index += 1
#         # the indices of the other vertices on this triangle
#         a, b = (local_index + 1)%3, (local_index + 2)%3
#         pivot, v1, v2 = vertices[local_index], vertices[a], vertices[b]
#         u, v = pivot - v1, pivot - v2
#
#         c = np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v) # -> cosine of the angle
#         angle = np.arccos(np.clip(c, -1, 1)) # if you really want the angle
#         return angle
#
#     # def colorize(analyticnormal, estimatednormal, minimum, maximum):
#     #     error = np.linalg.norm(np.cross(np.array(analyticnormal), np.array(estimatednormal)))
#     #     error = (error - minimum)/(maximum-minimum)
#     #     return colorize(error)
#
#     #compute the normal at vertex 'vertex_index' by simple average
#     def mean_interpolation(vertex_index):
#         normal = np.array([0, 0, 0]) #will accumulate
#         for triangle_index in neighborhood[vertex_index]:
#             normal = normal + compute_normal(triangle_index)
#         if len(neighborhood[vertex_index]) == 0:
#             print(vertex_index, "null neighboorhood")
#         normal = normalize(normal/len(neighborhood[vertex_index]))
#         return normal
#
#     def area_interpolation(vertex_index):
#         normal = np.array([0, 0, 0])
#         for triangle_index in neighborhood[vertex_index]:
#             normal = normal + compute_area(triangle_index) * compute_normal(triangle_index)
#         normal = normalize(normal)
#         return normal
#
#     def barycentric_area_interpolation(vertex_index):
#         normal = np.array([0, 0, 0])
#         for triangle_index in neighborhood[vertex_index]:
#             normal = normal + compute_barycentric_area(triangle_index, vertex_index) * compute_normal(triangle_index)
#         normal = normalize(normal)
#         return normal
#
#     def angle_interpolation(vertex_index):
#         normal = np.array([0, 0, 0])
#         for triangle_index in neighborhood[vertex_index]:
#             normal = normal + compute_angle(triangle_index, vertex_index) * compute_normal(triangle_index)
#         normal = normalize(normal)
#         return normal
#
#     expected_normals = [compute_analytic_normal(v) for v in vertices]
#     for method in methods:
#
#     print("will begin {} method".format(method.__name__))
#         normals = [method(index) for index in range(len(vertices)) ]
#         print("Will write results")
#         folder = os.path.join(results_folder, filename)
#         os.makedirs(folder, 0o777, True)
#         if should_colorize:
#             raw_errors = [np.linalg.norm(np.cross(expected_normals[i], np.array(normals[i]))) for i in range(len(normals))]
#             #colors = [colorize(expected_normals[index], normals[index]) for index in range(len(vertices))]
#             minimum = min(raw_errors)
#             if minimum < 0.00000001:
#                 minimum = 0.0
#             maximum = max(raw_errors)
#             #print(maximum, minimum)
#             if maximum - minimum > 0.1:
#                 raw_errors = [(error-minimum)/(maximum-minimum) for error in raw_errors]
#             colors = [colorize(error) for error in raw_errors]
#             write_PLY(os.path.join(folder, '{}_color_{}.ply'.format(filename, method.__name__)), vertices, indices, normals, colors)
#         else:
#             write_PLY(os.path.join(folder, '{}_{}.ply'.format(filename, method.__name__)), vertices, indices, normals)
#
#
#     # print("will begin area method")
#     # normals = [area_interpolation(index) for index in range(len(vertices)) ]
#     # print("Will write results")
#     # if should_colorize:
#     #     raw_errors = [np.linalg.norm(np.cross(np.array(expected_normals[i]), np.array(normals[i]))) for i in range(len(normals))]
#     #     #colors = [colorize(expected_normals[index], normals[index]) for index in range(len(vertices))]
#     #     minimum = min(raw_errors)
#     #     if minimum < 0.00000001:
#     #         minimum = 0.0
#     #     maximum = max(raw_errors)
#     #     #print(maximum, minimum)
#     #     if maximum - minimum > 0.1:
#     #         raw_errors = [(error-minimum)/(maximum-minimum) for error in raw_errors]
#     #     colors = [colorize(error) for error in raw_errors]
#     #     write_PLY(os.path.join(folder, '{}_color_area.ply'.format(filename)), vertices, indices, normals, colors)
#     # else:
#     #     write_PLY(os.path.join(folder, '{}_area.ply'.format(filename)), vertices, indices, normals)
#     #
#     # print("will begin barycentric method")
#     # normals = [barycentric_area_interpolation(index) for index in range(len(vertices)) ]
#     # print("Will write results")
#     # if should_colorize:
#     #     raw_errors = [np.linalg.norm(np.cross(np.array(expected_normals[i]), np.array(normals[i]))) for i in range(len(normals))]
#     #     #colors = [colorize(expected_normals[index], normals[index]) for index in range(len(vertices))]
#     #     minimum = min(raw_errors)
#     #     if minimum < 0.00000001:
#     #         minimum = 0.0
#     #     maximum = max(raw_errors)
#     #     #print(maximum, minimum)
#     #     if maximum - minimum > 0.1:
#     #         raw_errors = [(error-minimum)/(maximum-minimum) for error in raw_errors]
#     #     colors = [colorize(error) for error in raw_errors]
#     #     write_PLY(os.path.join(folder, '{}_color_barycentric.ply'.format(filename)), vertices, indices, normals, colors)
#     # else:
#     #     write_PLY(os.path.join(folder, '{}_barycentric.ply'.format(filename)), vertices, indices, normals)
#     #
#     # print("will begin angle method")
#     # normals = [angle_interpolation(index) for index in range(len(vertices)) ]
#     # print("Will write results")
#     # if should_colorize:
#     #     raw_errors = [np.linalg.norm(np.cross(np.array(expected_normals[i]), np.array(normals[i]))) for i in range(len(normals))]
#     #     #colors = [colorize(expected_normals[index], normals[index]) for index in range(len(vertices))]
#     #     minimum = min(raw_errors)
#     #     if minimum < 0.00000001:
#     #         minimum = 0.0
#     #     maximum = max(raw_errors)
#     #     #print(maximum, minimum)
#     #     if maximum - minimum > 0.1:
#     #         raw_errors = [(error-minimum)/(maximum-minimum) for error in raw_errors]
#     #     colors = [colorize(error) for error in raw_errors]
#     #     write_PLY(os.path.join(folder, '{}_color_angle.ply'.format(filename)), vertices, indices, normals, colors)
#     # else:
#     #     write_PLY(os.path.join(folder, '{}_angle.ply'.format(filename)), vertices, indices, normals)
#
#     print("will begin ground truth")
#     normals = [compute_analytic_normal(v) for v in vertices ]
#     print("Will write results")
#     if should_colorize:
#         raw_errors = [np.linalg.norm(np.cross(np.array(expected_normals[i]), np.array(normals[i]))) for i in range(len(normals))]
#         colors = [colorize(error) for error in raw_errors]
#         write_PLY(os.path.join(results_folder, '{}_color_groundtruth.ply'.format(filename)), vertices, indices, normals, colors)
#     else:
#         write_PLY(os.path.join(results_folder, '{}_groundtruth.ply'.format(filename)), vertices, indices, normals)
#
# def rec_all_data(folder):
#     filenames = [fname for fname in os.listdir(folder) if fname.endswith('.off')]
#     for fname in filenames:
#         interpolate_normals(os.path.join(folder, fname), cylinder_analytic_normal, True)
#
# if __name__ == '__main__':
#     # main()
#     # write_regular_cylinder_mesh()
#     # write_alternate_flip_cylinder_mesh()
#     # write_random_flip_cylinder_mesh()
#
#     create_data()
#     rec_all_data(data_folder)
#     #interpolate_normals("data/meshes/cylinder150_100_regular.off")
#compute the normal at vertex 'vertex_index' by simple average
# def mean_interpolation(neighborhood, vertex_index):
#     normal = np.array([0, 0, 0]) #will accumulate
#     for triangle_index in neighborhood[vertex_index]:
#         normal = normal + triangle_normal(triangle_index)
#     if len(neighborhood[vertex_index]) == 0:
#         print(vertex_index, "null neighboorhood")
#     normal = normalize(normal/len(neighborhood[vertex_index]))
#     return normal
