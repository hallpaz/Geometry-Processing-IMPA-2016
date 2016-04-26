import argparse
from collections import deque
from DataStructures import Triangle, Vertex
from model_operations import write_OFF, read_OFF
import os
import numpy as np


data_folder =  os.path.join("data", "subdivision")

def mesh_refinement_into_sphere(triangles, vertices, indices_map, num_iterations = 1):
    def index_of(vertex):
        key = str(vertex.x) + str(vertex.y) + str(vertex.z)
        if key in indices_map:
            return indices_map[key]
        else:
            index = len(indices_map)
            vertices.append(vertex)
            indices_map[key] = index
            return index

    refined_mesh = []
    for i in range(num_iterations):
        for t in triangles:
            # normalization to put over the unit sphere
            v1 = ((vertices[t.a] + vertices[t.b])/2).normalize()
            v2 = ((vertices[t.b] + vertices[t.c])/2).normalize()
            v3 = ((vertices[t.a] + vertices[t.c])/2).normalize()

            refined_mesh.append(Triangle( t.a, index_of(v1), index_of(v3)))
            refined_mesh.append(Triangle( t.b, index_of(v2), index_of(v1)))
            refined_mesh.append(Triangle( t.c, index_of(v3), index_of(v2)))
            refined_mesh.append(Triangle( index_of(v1), index_of(v2), index_of(v3)))
        triangles = refined_mesh

        write_OFF(os.path.join(data_folder, "esfera_iter" + str(5) + ".off"), vertices, triangles)
    return refined_mesh, vertices, indices_map

def surface(vertex):
    x = vertex.x
    y = vertex.y
    z = vertex.z
    return x*x*x*x - 5*x*x+ y*y*y*y - 5*y*y + z*z*z*z - 5*z*z + 11

def surface_gradient(vertex):
    x, y, z = vertex.x, vertex.y, vertex.z
    return Vertex((2*x*(2*x**2 - 5), 2*y*(2*y**2 - 5), 2*z*(2*z**2 - 5)))


def compute_normal(triangle_vertices):
    v1 = triangle_vertices[1] - triangle_vertices[0]
    v2 = triangle_vertices[2] - triangle_vertices[1]
    normal = Vertex.cross(v1, v2).normalize()
    return normal

def project_on_surface(vertex, normal, surface):
    t = 0.1
    a = surface(vertex)
    b = surface(vertex + t*normal)
    if a*b > 0 and abs(b) > abs(a):
        #wrong direction
        t = -t
        b = surface(vertex + t*normal)
    t = t/2
    middle = surface(vertex + t*normal)
    while(abs(middle) > 0.0001):
        


def mesh_refinement(triangles, vertices, indices_map, implicit_function, gradient, num_iterations = 1):
    def index_of(vertex, normal):
        key = str(vertex.x) + str(vertex.y) + str(vertex.z)
        if key in indices_map:
            return indices_map[key]
        else:
            index = len(indices_map)
            vertices.append(vertex)
            normals.append(normal)
            indices_map[key] = index
            return index

    refined_mesh = []
    for i in range(num_iterations):
        for t in triangles:
            # normalization to put over the unit sphere
            v1 = ((vertices[t.a] + vertices[t.b])/2)
            v2 = ((vertices[t.b] + vertices[t.c])/2)
            v3 = ((vertices[t.a] + vertices[t.c])/2)

            n1 = gradient(v1).normalize()
            n2 = gradient(v2).normalize()
            n3 = gradient(v3).normalize()

            #projection
            v1 = projection_on_surface(v1, n1, implicit_function)
            v2 = projection_on_surface(v2, n2, implicit_function)
            v3 = projection_on_surface(v3, n3, implicit_function)

            refined_mesh.append(Triangle( t.a, index_of(v1), index_of(v3)))
            refined_mesh.append(Triangle( t.b, index_of(v2), index_of(v1)))
            refined_mesh.append(Triangle( t.c, index_of(v3), index_of(v2)))
            refined_mesh.append(Triangle( index_of(v1), index_of(v2), index_of(v3)))
        triangles = refined_mesh

        write_OFF(os.path.join(data_folder, "projetada" + str(i) + ".off"), vertices, triangles)
    return refined_mesh, vertices, indices_map


def main(input_file, num_iterations=1):
    vertices, triangles = read_OFF(input_file)
    vertices = [v.normalize() for v in vertices]
    indices_map = { str(v.x) + str(v.y) + str(v.z) : i for i in range(len(vertices)) for v in vertices }
    mesh_refinement_into_sphere(triangles, vertices, indices_map, num_iterations)

    #mesh_refinement(triangles, vertices, indices_map, surface, num_iterations)


def compute_vertex_normals(normals):
    num = 0
    v = Vertex(0, 0, 0)
    for n in normals:
        num += 1
        v = v + n

    return v/n

def main(input_file, num_iterations=1):
    vertices, triangles = read_OFF(input_file)
    vertices = [project(v) for v in vertices]

    indices_map = { str(v.x) + str(v.y) + str(v.z) : i for i in range(len(vertices)) for v in vertices }
    #mesh_refinement_into_sphere(triangles, vertices, indices_map, num_iterations)
    # triangles_normals = [compute_normal([ vertices[t.a, vertices[t.b], vertices[t.c]) for t in triangles]
    # vertex_normals = []
    # for i in range(len(vertices)):
    #     triangle_numbers = [j for j in range(len(triangles)) if vertices[i] in [vertices[indices[j].a], vertices[indices[j].b], vertices[indices[j].c]]]
    #     vertex_normals.append(compute_vertex_normals([ normals[m] for m in triangle_numbers ] ))
    mesh_refinement(triangles, vertices, vertex_normals, indices_map, surface, num_iterations)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) ... The output is written to a txt file
    ''')
    parser.add_argument('input_file', help='input list of points with X Y per line (format: txt)')
    args = parser.parse_args()

    main(os.path.join(data_folder, args.input_file))
