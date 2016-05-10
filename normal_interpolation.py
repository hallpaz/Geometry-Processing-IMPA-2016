import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import math
import random

from geometrypack.parametric_functions import cylinder
from geometrypack.dataio import read_OFF, write_OFF


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
            indices.append([p*meridians + m, (p+1)*meridians + m, p*meridians + (m+1)%meridians]) #superior
            indices.append([p*meridians + (m+1)%meridians, (p+1)*meridians + m, (p+1)*meridians + (m+1)%meridians ]) #inferior

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
                indices.append([p*meridians + m, (p+1)*meridians + m, p*meridians + (m+1)%meridians]) #superior
                indices.append([p*meridians + (m+1)%meridians, (p+1)*meridians + m, (p+1)*meridians + (m+1)%meridians ]) #inferior
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
                indices.append([p*meridians + m, (p+1)*meridians + m, p*meridians + (m+1)%meridians]) #superior
                indices.append([p*meridians + (m+1)%meridians, (p+1)*meridians + m, (p+1)*meridians + (m+1)%meridians ]) #inferior
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
    write_OFF("cilindro.off", vertices, indices)


def interpolate_normals(filename, method):
    vertices, indices = read_OFF(filename)

    neighborhood = [ np.array([t for t in range(len(indices) if index in indices[t])]) for index in range(len(vertices)) ]

    # triangle is a list of 3 numpy array
    def compute_normal(triangle_index):
        t = indices[triangle_index]
        normal = np.cross(vertices[t[1]] - vertices[t[0]], vertices[t[2]] - vertices[t[1]])
        return np.normalize(normal)

    def compute_area(triangle_index):
        t = indices[triangle_index]
        a, b, c = abs(vertices[t[0]] - vertices[t[1]]), abs(vertices[t[1]] - vertices[t[2]]), abs(vertices[t[2]] - vertices[t[0]])
        # calculate the semi-perimeter
        s = (a + b + c) / 2
        # calculate the area
        area = (s*(s-a)*(s-b)*(s-c)) ** 0.5
        return area

    def compute_angle(triangle_index, vertex_index):
        t = indices[triangle_index]

    def mean_interpolation(vertex_index):
        normal = np.array([0, 0, 0])
        for triangle_index in neighborhood:
            normal += compute_normal(indices[triangle_index])
        normal = np.normalize(normal/len(neighborhood))
        return normal

    def area_interpolation(vertex_index):
        normal = np.array([0, 0, 0])
        for triangle_index in neighborhood:
            current_triangle = indices[triangle_index]
            normal += compute_area(current_triangle) * compute_normal(current_triangle)
        normal = np.normalize(normal/len(neighborhood))
        return normal

    def barycentric_area_interpolation(vertex_index):
        normal = np.array([0, 0, 0])
        for triangle_index in neighborhood:
            current_triangle = indices[triangle_index]
            pivot = current_triangle.index(vertex_index)
            area = 0

            barycenter = (current_triangle[0] + current_triangle[1] + current_triangle[2])/3
            a = vertices[current_triangle[pivot]]
            b = ()/2
            c =

            # calculate the semi-perimeter
            s = (a + b + c) / 2
            # calculate the area
            area += (s*(s-a)*(s-b)*(s-c)) ** 0.5

            normal += compute_area(current_triangle) * compute_normal(current_triangle)
        normal = np.normalize(normal/len(neighborhood))
        return normal

    def angle_interpolation(vertex_index):
        normal = np.array([0, 0, 0])
        for triangle_index in neighborhood:
            current_triangle = indices[triangle_index]
            normal += compute_angle(current_triangle, vertex_index) * compute_normal(current_triangle)
        normal = np.normalize(normal/len(neighborhood))
        return normal)

    if method == "simple":
        pass
    elif method == ""

if __name__ == '__main__':
    main()
    write_regular_cylinder_mesh()
    write_alternate_flip_cylinder_mesh()
    write_random_flip_cylinder_mesh()
