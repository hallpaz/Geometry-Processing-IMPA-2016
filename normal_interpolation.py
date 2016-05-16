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



#calcula a normal para cada triângulo
# para cada vértice, identifique todos os triângulos que possuem esse vértice
# faça uma média ponderara dos valores da normal de cada triângulo que contém
# o vértice para enontrar a normal do vertice

def interpolate_normals(filename, method):
    vertices, indices = read_OFF(filename)
    print("BEGIN", vertices[0], vertices[-1])
    print("will compute neighbohood")
    # neighborhood[i] é uma lista que contém os índices de todos os triângulos
    # que contém o vértice i
    neighborhood = [ [] for i in range(len(vertices)) ]
    for index in range(len(indices)):
        # t recebe os 3 indices do triangulo da posição 'index'
        t = indices[index]
        neighborhood[t[0]].append(index)
        neighborhood[t[1]].append(index)
        neighborhood[t[2]].append(index)

    print("computed neighborhood")

    # triangle is a list of 3 numpy array
    def compute_normal(triangle_index):
        t = indices[triangle_index]
        a, b, c = t[0], t[1], t[2]

        normal = np.cross(vertices[b] - vertices[a], vertices[c] - vertices[b])
        return normalize(normal)

    def compute_area(triangle_index):
        t = indices[triangle_index]
        a = np.linalg.norm(vertices[t[0]] - vertices[t[1]])
        b = np.linalg.norm(vertices[t[1]] - vertices[t[2]])
        c = np.linalg.norm(vertices[t[2]] - vertices[t[0]])
        # calculate the semi-perimeter
        s = (a + b + c) / 2
        # calculate the area
        area = (s*(s-a)*(s-b)*(s-c)) ** 0.5
        return area

    def compute_barycentric_area(triangle_index, pivot_index):
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

    def compute_angle(triangle_index, pivot_index):
        local_index = 0
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

    #compute the normal at vertex 'vertex_index' by simple average
    def mean_interpolation(vertex_index):
        normal = np.array([0, 0, 0]) #will accumulate
        for triangle_index in neighborhood[vertex_index]:
            normal = normal + compute_normal(triangle_index)
        normal = normalize(normal/len(neighborhood[vertex_index]))
        return normal

    def area_interpolation(vertex_index):
        normal = np.array([0, 0, 0])
        for triangle_index in neighborhood[vertex_index]:
            normal = normal + compute_area(triangle_index) * compute_normal(triangle_index)
        normal = normalize(normal)
        return normal

    def barycentric_area_interpolation(vertex_index):
        normal = np.array([0, 0, 0])
        for triangle_index in neighborhood[vertex_index]:
            normal = normal + compute_barycentric_area(triangle_index, vertex_index) * compute_normal(triangle_index)
        normal = normalize(normal)
        return normal

    def angle_interpolation(vertex_index):
        normal = np.array([0, 0, 0])
        for triangle_index in neighborhood[vertex_index]:
            normal = normal + compute_angle(triangle_index, vertex_index) * compute_normal(triangle_index)
        normal = normalize(normal)
        return normal

    # print("will begin mean method")
    # normals = [mean_interpolation(index) for index in range(len(vertices)) ]
    # print(normals[0])
    # print("Will write results")
    # write_PLY('regular_mean.ply', vertices, indices, normals)
    #
    #
    # print("will begin area method")
    # normals = [area_interpolation(index) for index in range(len(vertices)) ]
    # print("Will write results")
    # write_PLY('regular_area.ply', vertices, indices, normals)
    #
    # print("will begin barycentric method")
    # normals = [barycentric_area_interpolation(index) for index in range(len(vertices)) ]
    # print("Will write results")
    # write_PLY('regular_barycentric.ply', vertices, indices, normals)

    print("will begin angle method")
    normals = [angle_interpolation(index) for index in range(len(vertices)) ]
    print("Will write results")
    write_PLY('regular_angle.ply', vertices, indices, normals)



if __name__ == '__main__':
    # main()
    # write_regular_cylinder_mesh()
    # write_alternate_flip_cylinder_mesh()
    # write_random_flip_cylinder_mesh()

    interpolate_normals("data/meshes/cylinder150_100_regular.off", "mean")
