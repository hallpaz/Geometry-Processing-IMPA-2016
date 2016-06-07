import argparse
import numpy as np
import os

import matplotlib.pyplot as plt
from math import pi
from geometrypack.dataio import read_OFF, write_PLY, write_OFF, write_uv_PLY
from geometrypack.algorithms import distance_vertex_plane, normalize
#from geometrypack.drawing import colorize
from geometrypack import drawing
from geometrypack.meshes import compute_neighborhood, compute_boundary_indices, rotate
#from mesh.soup_to_mesh import compute_boundary_indices

# JUST IMPORTED COMPUTE BOUNDARY INDICES
BIG_RADIUS = 10

def maxnormalize(v):
    norm = max(v)
    if norm < 0.0000001:
        return np.array([0, 0])
    return v/norm

def bounding_square(vertices):
    BIG = 9999999
    xmax = ymax = -BIG
    xmin = ymin = BIG
    for i in range(len(vertices)):
        if vertices[i][0] > xmax:
            xmax = vertices[i][0]
        if vertices[i][1] > ymax:
            ymax = vertices[i][1]
        if vertices[i][0] < xmin:
            xmin = vertices[i][0]
        if vertices[i][1] < ymin:
            ymin = vertices[i][1]
    return (xmin, ymin, xmax, ymax)

def uvcoordinates(vertices, box):
    width = box[2] - box[0]
    height = box[3] - box[1]
    uv = [((vertices[i][0] - box[0])/width, (vertices[i][1] - box[1])/height) for i in range(len(vertices))]
    return uv

def laplacian_adjust(vertices, neighborhood, border_indices, iterations=1):
    '''average vertices but keep borders fixed'''

    for k in range(iterations):
        for index in range(len(vertices)):
            if index in border_indices:
                #print("border", index)
                continue
            neighbors = neighborhood[index]
            average = np.array([0, 0])
            for n in neighbors:
                average = average + vertices[n]
            vertices[index] = (vertices[index] + average)/(len(neighbors)+1)
    return vertices

def stretch_border(vertices, indices):
    strectched_vertices = np.array([[v[0], v[1]] for v in vertices])
    border_indices = compute_boundary_indices(vertices, indices)
    neighborhood = compute_neighborhood(vertices, indices)
    radius = BIG_RADIUS

    seen_vertices_index = []
    old_ring = border_indices
    while(len(seen_vertices_index) < len(vertices)):
        current_ring = []
        for i in old_ring:
            current_ring.extend(neighborhood[i])
        # filter the vertices that were already seen
        current_ring = [index for index in current_ring if index not in seen_vertices_index]
        # remove duplicates
        current_ring = list(set(current_ring))
        for i in current_ring:
            strectched_vertices[i] = normalize(strectched_vertices[i]) * radius
        seen_vertices_index.extend(current_ring)
        old_ring = current_ring
        radius -= 1

    ax = plt.axes()
    drawing.plot_and_save("{}.png".format("rings"), False, ax, vertices=strectched_vertices, triangles=indices)
    drawing.plot_and_save("{}.png".format("rings"), True, ax, vertices=np.array([strectched_vertices[i] for i in border_indices]), vertices_color="red")
    # for index in border_indices:
    #     strectched_vertices[index] = maxnormalize(strectched_vertices[index])*BIG_RADIUS
    # drawing.plot_and_save("{}.png".format("square_rings"), False, ax, vertices=strectched_vertices, triangles=indices)
    # drawing.plot_and_save("{}.png".format("square_rings"), True, ax, vertices=np.array([strectched_vertices[i] for i in border_indices]), vertices_color="red")
    # Optional
    #reduce_hole(strectched_vertices, neighborhood, border_indices, old_ring)

    for k in range(5):
        ax = plt.axes()
        drawing.plot_and_save("manequin/{}.png".format(k), False, ax, vertices=strectched_vertices, triangles=indices)
        drawing.plot_and_save("manequin/{}.png".format(k), True, ax, vertices=np.array([strectched_vertices[i] for i in border_indices]), vertices_color="red")
        strectched_vertices = laplacian_adjust(strectched_vertices, neighborhood, border_indices, 2)
        print("Ok", k)
    return strectched_vertices

def projection(iterations):
    vertices, indices = read_OFF("data/meshes/manequin.off")
    border_indices = compute_boundary_indices(vertices, indices)

    neighborhood = compute_neighborhood(vertices, indices)
    projected_vertices = np.array([[v[0], v[1]] for v in vertices])
    for k in range(iterations):
        ax = plt.axes()
        # drawing.plot_and_save("manequin/{}.png".format(k), False, ax, vertices=projected_vertices, triangles=indices)
        # drawing.plot_and_save("manequin/{}.png".format(k), True, ax, vertices=np.array([projected_vertices[i] for i in border_indices]), vertices_color="red")
        projected_vertices = laplacian_adjust(projected_vertices, neighborhood, border_indices, 1000)
        print("Ok", k)

    bbox = bounding_square(projected_vertices)
    uv = uvcoordinates(projected_vertices, bbox)
    with open("uv.txt", "w") as uvfile:
        for coordinate in uv:
            uvfile.write("{} {}\n".format(coordinate[0], coordinate[1]))
    write_uv_PLY("uvmanequin.ply", vertices, indices, uv)

def square_projection():
    vertices, indices = read_OFF("data/meshes/manequin.off")
    border_indices = compute_boundary_indices(vertices, indices)

    neighborhood = compute_neighborhood(vertices, indices)
    projected_vertices = np.array([[v[0], v[1]] for v in vertices])
    for index in border_indices:
        projected_vertices[index] = maxnormalize(projected_vertices[index])
    for k in range(1000):
        ax = plt.axes()
        drawing.plot_and_save("manequin/{}.png".format(k), False, ax, vertices=projected_vertices, triangles=indices)
        drawing.plot_and_save("manequin/{}.png".format(k), True, ax, vertices=np.array([projected_vertices[i] for i in border_indices]), vertices_color="red")
        projected_vertices = laplacian_adjust(projected_vertices, neighborhood, border_indices)
        print("Ok", k)

def circle_projection():
    vertices, indices = read_OFF("data/meshes/manequin.off")
    border_indices = compute_boundary_indices(vertices, indices)

    neighborhood = compute_neighborhood(vertices, indices)
    projected_vertices = np.array([[v[0], v[1]] for v in vertices])
    for index in border_indices:
        projected_vertices[index] = normalize(projected_vertices[index])*10
    for k in range(1000):
        ax = plt.axes()
        drawing.plot_and_save("manequin/{}.png".format(k), False, ax, vertices=projected_vertices, triangles=indices)
        drawing.plot_and_save("manequin/{}.png".format(k), True, ax, vertices=np.array([projected_vertices[i] for i in border_indices]), vertices_color="red")
        projected_vertices = laplacian_adjust(projected_vertices, neighborhood, border_indices)
        print("Ok", k)


def main():
    vertices, indices = read_OFF("data/meshes/manequin90.off")
    border_indices = compute_boundary_indices(vertices, indices)

    #strectched_vertices = project_and_adjust(vertices, indices, border_indices)
    projected_vertices = project_and_adjust(vertices, indices, border_indices)
    #uv = compute_uvcoordinates_onsquare(projected_vertices, 1.0)
    uv = projected_vertices
    for coordinates in uv:
        if max(coordinates) > 1 or min(coordinates) < 0:
            print('range of values is inconsistent', coordinates)
            pass

def use_stretch():
    vertices, indices = read_OFF("data/meshes/manequin.off")
    strectched_vertices = stretch_border(vertices, indices)
    bbox = bounding_square(strectched_vertices)
    uv = uvcoordinates(strectched_vertices, bbox)
    with open("uv.txt", "w") as uvfile:
        for coordinate in uv:
            uvfile.write("{} {}\n".format(coordinate[0], coordinate[1]))
    write_uv_PLY("uvmanequin.ply", vertices, indices, uv)

def mark_border(filename):
    vertices, indices = read_OFF(filename)
    colors = [(0, 0, 255) for v in vertices]
    normals = [[0, 0, 0]]*len(vertices)
    border_indices = compute_boundary_indices(vertices, indices)
    print(len(border_indices))
    for i in border_indices:
        colors[i] = (255, 0, 0)
        print(vertices[i])
    name = os.path.splitext(os.path.basename(filename))[0]
    write_PLY("marked_" + name + ".ply", vertices, indices, normals, colors)

if __name__ == '__main__':
    # vertices, indices = read_OFF("data/meshes/manequin.off")
    # vertices = rotate(vertices, [1, 0, 0], -pi/2)
    # write_OFF("data/meshes/manequin90.off", vertices, indices)
    # mark_border('data/meshes/manequin90.off')
    #main()
    # square_projection()
    # circle_projection()
    # projection(1)
    use_stretch()
