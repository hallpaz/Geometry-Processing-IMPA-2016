import argparse
import numpy as np
import os

import matplotlib.pyplot as plt
from math import pi, sqrt
from geometrypack.dataio import read_OFF, write_PLY, write_OFF, write_uv_PLY
from geometrypack.algorithms import distance_vertex_plane, normalize
#from geometrypack.drawing import colorize
from geometrypack import drawing
from geometrypack.meshes import compute_neighborhood, compute_boundary_indices, rotate
#from mesh.soup_to_mesh import compute_boundary_indices

# JUST IMPORTED COMPUTE BOUNDARY INDICES
BIG_RADIUS = 10

def compute_area(tvertices):
    a = np.linalg.norm(tvertices[0] - tvertices[1])
    b = np.linalg.norm(tvertices[1] - tvertices[2])
    c = np.linalg.norm(tvertices[2] - tvertices[0])
    # calculate the semi-perimeter
    s = (a + b + c) / 2
    # calculate the area
    area = (s*(s-a)*(s-b)*(s-c)) ** 0.5
    return area

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

def edgesize_adjust(vertices3D, vertices2D, neighborhood, border_indices, iterations=1):
    '''average vertices but keep borders fixed'''

    for k in range(iterations):
        for index in range(len(vertices2D)):
            if index in border_indices:
                #print("border", index)
                continue
            neighbors = neighborhood[index]
            average = np.array([0, 0])
            weightsum = 0
            for n in neighbors:
                w = np.linalg.norm(vertices3D[index] - vertices3D[n])
                average = average + w*vertices2D[n]
                weightsum += w
            vertices2D[index] = average/weightsum
    return vertices2D

def edgeadjust3D(vertices3D, neighborhood, border_indices, iterations=1):
    for k in range(iterations):
        for index in range(len(vertices3D)):
            if index in border_indices:
                #print("border", index)
                continue
            neighbors = neighborhood[index]
            average = np.array([0, 0, 0])
            weightsum = 0
            for n in neighbors:
                w = np.linalg.norm(vertices3D[index] - vertices3D[n])
                average = average + w*vertices3D[n]
                weightsum += w
            vertices3D[index] = average/weightsum
    return vertices3D


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

def edgesizeprojection(iterations):
    vertices, indices = read_OFF("data/meshes/manequin.off")
    border_indices = compute_boundary_indices(vertices, indices)

    neighborhood = compute_neighborhood(vertices, indices)
    projected_vertices = np.array([[v[0], v[1]] for v in vertices])
    for k in range(iterations):
        ax = plt.axes()
        drawing.plot_and_save("manequin_es/{}.png".format(k), False, ax, vertices=projected_vertices, triangles=indices)
        drawing.plot_and_save("manequin_es/{}.png".format(k), True, ax, vertices=np.array([projected_vertices[i] for i in border_indices]), vertices_color="red")
        projected_vertices = edgesize_adjust(vertices, projected_vertices, neighborhood, border_indices, 1)
        print("Ok", k)

    bbox = bounding_square(projected_vertices)
    uv = uvcoordinates(projected_vertices, bbox)
    with open("uv_es1000.txt", "w") as uvfile:
        for coordinate in uv:
            uvfile.write("{} {}\n".format(coordinate[0], coordinate[1]))
    write_uv_PLY("uvmanequin_es.ply", vertices, indices, uv)

def onsurfacecomputation(iterations):
    vertices, indices = read_OFF("data/meshes/manequin.off")
    border_indices = compute_boundary_indices(vertices, indices)

    neighborhood = compute_neighborhood(vertices, indices)
    projected_vertices = np.array([[v[0], v[1], v[2]] for v in vertices])
    for k in range(iterations):
        ax = plt.axes()
        drawing.plot_and_save("manequin_sdt/{}.png".format(k), False, ax, vertices=projected_vertices, triangles=indices)
        drawing.plot_and_save("manequin_sdt/{}.png".format(k), True, ax, vertices=np.array([projected_vertices[i] for i in border_indices]), vertices_color="red")
        projected_vertices = edgeadjust3D(projected_vertices, neighborhood, border_indices, 1)
        print("Ok", k)

    # now, we do project on the plane
    # projected_vertices = np.array([[v[0], v[1]] for v in projected_vertices])
    bbox = bounding_square(projected_vertices)
    uv = uvcoordinates(projected_vertices, bbox)
    with open("uv_3Dpj1000.txt", "w") as uvfile:
        for coordinate in uv:
            uvfile.write("{} {}\n".format(coordinate[0], coordinate[1]))
    write_uv_PLY("uvmanequin_3Dpj.ply", vertices, indices, uv)
    try:
        write_uv_PLY("manequin3Ddistorted.ply", projected_vertices, indices, uv)
    except Exception as e:
        print(e)

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


def triangle_stretch_measure(p1, p2, p3, q1, q2, q3):
    '''p is the point in the plane domain. q is the vertex on the surface domain'''

    A = ((p2[0]-p1[0])*(p3[1]-p1[1]) - (p3[0]-p1[0])*(p2[1]-p1[1]))/2
    Ss = (q1*(p2[1]-p3[1]) + q2*(p3[1]-p1[1]) + q3*(p1[1]-p2[1]))/(2*A)
    St = (q1*(p3[0]-p2[0]) + q2*(p1[0]-p3[0]) + q3*(p2[0]-p1[0]))/(2*A)

    a = np.dot(Ss, Ss)
    b = np.dot(Ss, St)
    c = np.dot(St, St)

    delta = sqrt((a-c)**2 + 4*b*b)
    big_gamma_squared = ((a+c) + delta)/2
    small_gamma_squared = ((a+c) - delta)/2
    stretchvalue = sqrt((big_gamma_squared+small_gamma_squared)/2)
    return stretchvalue

def vertex_stretch_measure(vertices3D, vertices2D, surrounding_triangles):
    areasum = 0
    numsum = 0
    for t in surrounding_triangles:
        tvertices3D = [vertices3D[t[0]], vertices3D[t[1]], vertices3D[t[2]]]
        tvertices2D = [vertices2D[t[0]], vertices2D[t[1]], vertices2D[t[2]]]
        area = compute_area(tvertices3D)
        areasum += area
        # print(type(tvertices2D[0]), type(tvertices3D[0]))
        tsmeasure = triangle_stretch_measure(tvertices2D[0], tvertices2D[1], tvertices2D[2],
                                            tvertices3D[0], tvertices3D[1], tvertices3D[2])
        numsum += area*tsmeasure*tsmeasure
    stretchvalue = sqrt(numsum/areasum)
    return stretchvalue

def baseprojection(vertices, indices, iterations):

    border_indices = compute_boundary_indices(vertices, indices)

    neighborhood = compute_neighborhood(vertices, indices)
    projected_vertices = np.array([[v[0], v[1]] for v in vertices])

    projected_vertices = edgesize_adjust(vertices, projected_vertices, neighborhood, border_indices, iterations)
    ax = plt.axes()
    drawing.plot_and_save("manequin_base/{}.png".format(k), False, ax, vertices=projected_vertices, triangles=indices)
    drawing.plot_and_save("manequin_base/{}.png".format(k), True, ax, vertices=np.array([projected_vertices[i] for i in border_indices]), vertices_color="red")

    print("Ok base")
    return projected_vertices


def stretch_minimizer():
    iterations = 7
    vertices, indices = read_OFF("data/meshes/manequin.off")
    projected_vertices = np.array([[v[0], v[1]] for v in vertices])
    neighborhood = compute_neighborhood(vertices, indices)
    surrounding_triangles = [[t for t in indices if index in t] for index in range(len(vertices))]
    border_indices = compute_boundary_indices(vertices, indices)

    # base parametrization
    k = 777
    projected_vertices = edgesize_adjust(vertices, projected_vertices, neighborhood, border_indices, k)
    # hum...
    # bbox = bounding_square(projected_vertices)
    # projected_vertices = np.array(uvcoordinates(projected_vertices, bbox))

    ax = plt.axes()
    drawing.plot_and_save("manequin_base/{}.png".format(k), False, ax, vertices=projected_vertices, triangles=indices)
    drawing.plot_and_save("manequin_base/{}.png".format(k), True, ax, vertices=np.array([projected_vertices[i] for i in border_indices]), vertices_color="red")
    print("will start method")
    # uniform weights to begin
    weights = [1 for i in range(len(vertices))]
    for k in range(iterations):
        weights = [weights[j]/vertex_stretch_measure(vertices, projected_vertices, surrounding_triangles[j]) for j in range(len(weights))]
        for index in range(len(projected_vertices)):
            if index in border_indices:
                continue
            neighbors = neighborhood[index]
            average = np.array([0, 0])
            weightsum = 0
            # print(surrounding_triangles[index])
            for n in neighbors:
                # weights[n] = weights[n]/vertex_stretch_measure(vertices, projected_vertices, surrounding_triangles[n])
                average = average + weights[n]*projected_vertices[n]
                weightsum += weights[n]
            projected_vertices[index] = average/weightsum
        ax = plt.axes()
        drawing.plot_and_save("manequin_stretch/{}.png".format(k), False, ax, vertices=projected_vertices, triangles=indices)
        drawing.plot_and_save("manequin_stretch/{}.png".format(k), True, ax, vertices=np.array([projected_vertices[i] for i in border_indices]), vertices_color="red")

    bbox = bounding_square(projected_vertices)
    uv = uvcoordinates(projected_vertices, bbox)
    with open("uv_stretch1000.txt", "w") as uvfile:
        for coordinate in uv:
            uvfile.write("{} {}\n".format(coordinate[0], coordinate[1]))
    write_uv_PLY("uvmanequin_stretch.ply", vertices, indices, uv)

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
    # use_stretch()
    # edgesizeprojection(1000)
    # onsurfacecomputation(1000)
    stretch_minimizer()
