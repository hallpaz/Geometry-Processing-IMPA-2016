import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import math
import random

from geometrypack.parametric_functions import cylinder
from geometrypack.dataio import read_OFF, write_OFF, write_PLY


def normalize(v):
    norm=np.linalg.norm(v)
    if norm < 0.00001:
       return v
    return v/norm

def sample_cylinder_points(meridians, parallels, height=1):
    inc_height = height/parallels
    inc_rad = 2*math.pi/meridians
    vertices = []
    h = 0
    while(h <= height):
        angle = 0
        for m in range(meridians):
            vertices.append(cylinder( angle, h ))
            angle += inc_rad
        h+= inc_height

    return vertices



def write_regular_cylinder_mesh():
    meridians = 150
    parallels = 100
    height = 2

    vertices = sample_cylinder_points(meridians, parallels, height)
    indices = []
    for p in range(parallels-1):
        for m in range(meridians):
            indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
            indices.append([p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m ]) #inferior

    indices = np.array(indices)

    write_OFF("cylinder{}_{}_regular.off".format(meridians, parallels), vertices, indices)
    print("OFF escrito")

def write_alternate_flip_cylinder_mesh():
    meridians = 150
    parallels = 100
    height = 2

    vertices = sample_cylinder_points(meridians, parallels, height)
    indices = []
    for p in range(parallels-1):
        for m in range(meridians):
            if p%2 == 0:
                indices.append([p*meridians + m, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
                indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians ]) #inferior
            else:
                indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
                indices.append([p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m ]) #inferior

    indices = np.array(indices)

    write_OFF("cylinder{}_{}_alternateflip.off".format(meridians, parallels), vertices, indices)
    print("OFF escrito")

def write_random_flip_cylinder_mesh():
    meridians = 150
    parallels = 100
    height = 2

    vertices = sample_cylinder_points(meridians, parallels, height)
    indices = []
    for p in range(parallels-1):
        for m in range(meridians):
            k = random.randint(0, 1)
            if k%2 == 0:
                indices.append([p*meridians + m, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
                indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians ]) #inferior
            else:
                indices.append([p*meridians + m, p*meridians + (m+1)%meridians, (p+1)*meridians + m]) #superior
                indices.append([p*meridians + (m+1)%meridians, (p+1)*meridians + (m+1)%meridians, (p+1)*meridians + m ]) #inferior

    indices = np.array(indices)

    write_OFF("cylinder{}_{}_randomflip.off".format(meridians, parallels), vertices, indices)
    print("OFF escrito")

def main():
    vertices = sample_cylinder_points(150, 100, 2)
    indices = []
    print("vai escrever")
    write_OFF("cylinder_150_100_delaunay.off", vertices, indices)


def interpolate_normals(filename, method):
    vertices, indices = read_OFF(filename)

    print("will compute neighbohood")
    neighborhood = [ [] for i in range(len(vertices)) ]
    for index in range(len(indices)):
        t = indices[index]
        neighborhood[t[0]].append(index)
        neighborhood[t[1]].append(index)
        neighborhood[t[2]].append(index)

    print("computed neighborhood")

    # triangle is a list of 3 numpy array
    def compute_normal(triangle_index):
        t = indices[triangle_index]
        normal = np.cross(vertices[t[1]] - vertices[t[0]], vertices[t[2]] - vertices[t[1]])
        return normalize(normal)

    def compute_area(triangle_index):
        t = indices[triangle_index]
        a, b, c = abs(vertices[t[0]] - vertices[t[1]]), abs(vertices[t[1]] - vertices[t[2]]), abs(vertices[t[2]] - vertices[t[0]])
        # calculate the semi-perimeter
        s = (a + b + c) / 2
        # calculate the area
        area = (s*(s-a)*(s-b)*(s-c)) ** 0.5
        return area

    # def compute_angle(triangle_index, vertex_index):
    #     t = indices[triangle_index]

    def mean_interpolation(vertex_index):
        normal = np.array([0, 0, 0])
        for triangle_index in neighborhood[vertex_index]:
            normal = normal + compute_normal(indices[triangle_index])
        normal = normalize(normal/len(neighborhood[vertex_index]))
        return normal

    def area_interpolation(vertex_index):
        normal = np.array([0, 0, 0])
        for triangle_index in neighborhood[vertex_index]:
            current_triangle = indices[triangle_index]
            normal = normal + compute_area(current_triangle) * compute_normal(current_triangle)
        normal = normalize(normal/len(neighborhood[vertex_index]))
        return normal

    # def barycentric_area_interpolation(vertex_index):
    #     normal = np.array([0, 0, 0])
    #     for triangle_index in neighborhood:
    #         current_triangle = indices[triangle_index]
    #         pivot = current_triangle.index(vertex_index)
    #         area = 0
    #
    #         barycenter = (current_triangle[0] + current_triangle[1] + current_triangle[2])/3
    #         a = vertices[current_triangle[pivot]]
    #         b = (0)/2
    #         c = 0
    #
    #         # calculate the semi-perimeter
    #         s = (a + b + c) / 2
    #         # calculate the area
    #         area += (s*(s-a)*(s-b)*(s-c)) ** 0.5
    #
    #         normal += compute_area(current_triangle) * compute_normal(current_triangle)
    #     normal = normalize(normal/len(neighborhood))
    #     return normal

    # def angle_interpolation(vertex_index):
    #     normal = np.array([0, 0, 0])
    #     for triangle_index in neighborhood:
    #         current_triangle = indices[triangle_index]
    #         normal += compute_angle(current_triangle, vertex_index) * compute_normal(current_triangle)
    #     normal = normalize(normal/len(neighborhood))
    #     return normal


    print("will begin mean method")
    normals = [mean_interpolation(index) for index in range(len(vertices)) ]
    print(normals[0])
    print("Will write results")
    write_PLY('regular_mean.ply', vertices, indices, normals)


    print("will begin area method")
    normals = [area_interpolation(index) for index in range(len(vertices)) ]
    print("Will write results")
    write_PLY('regular_area.ply', vertices, indices, normals)

    # print("will begin mean method")
    # normals = [mean_interpolation(index) for index in range(len(vertices)) ]
    # print("Will write results")
    # write_PLY('teste_mean.ply', vertices, indices, normals)


    # if method == "mean":
    #     print("will begin mean method")
    #     normals = [mean_interpolation(index) for index in range(len(vertices)) ]
    #     print("Will write results")
    #     write_PLY('teste_mean.ply', vertices, indices, normals)
    # elif method == "area":
    #     pass
    # else:
    #     pass

if __name__ == '__main__':
    # main()
    # write_regular_cylinder_mesh()
    # write_alternate_flip_cylinder_mesh()
    # write_random_flip_cylinder_mesh()

    interpolate_normals("data/meshes/cylinder150_100_regular.off", "mean")
