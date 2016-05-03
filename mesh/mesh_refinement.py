import argparse
from collections import deque
from DataStructures import Triangle, Vertex
from model_operations import write_OFF, read_OFF
import os
import numpy as np


data_folder =  os.path.join("data", "subdivision")
projection_aproximation_tolerance = 0.001

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
    return Vertex( 4*(x**3) - 5*x, 4*(y**3) - 5*y, 4*(z**3) - 5*z )


def compute_normal(triangle_vertices):
    v1 = triangle_vertices[1] - triangle_vertices[0]
    v2 = triangle_vertices[2] - triangle_vertices[1]
    normal = Vertex.cross(v1, v2).normalize()
    return normal


def projection_on_surface(vertex, normal, surface):
    t = 0.1
    begin = vertex
    end = vertex + normal.scalar_mult(t)
    a = surface(begin)
    b = surface(end)
    i = 0
    if a*b > 0 and abs(b) > abs(a):
        print(a, b, "afastando")
        if abs(b) > abs(a):
            #wrong direction
            t = -t
            #print("inverteu: " + str(t))
            end = begin + normal.scalar_mult(t)
            b = surface(end)

    while a*b > 0:

        if abs(b) < abs(a):
            begin = end
            a = surface(begin)
            end = begin + normal.scalar_mult(t)
            b = surface(end)
        else:
            while a*b > 0:
                middle = (begin + end)/2
                c = surface(middle)
                if a*c < 0:
                    end = middle
                    b = c
                    break
                else:
                    if abs(c) > abs(a):
                        end = middle
                        b = c
                    else:
                        begin = middle
                        a = c
            print(a, b)

    middle = (begin + end)/2
    c = surface(middle)
    while(abs(c) > projection_aproximation_tolerance):
        #print(c)
        if c*a < 0:
            end = middle
        else:
            begin = middle
        a = surface(begin)
        middle = (end + begin)/2
        c = surface(middle)

    return middle


# def projection_on_surface(vertex, normal, surface):
#     t = 0.1
#     begin = vertex
#     end = vertex + normal.scalar_mult(t)
#     a = surface(begin)
#     b = surface(end)
#     i = 0
#     #print('iniciou: ' + str(t))
#     if a*b > 0 and abs(b) > abs(a):
#         #wrong direction
#         t = -t
#         print("inverteu: " + str(t))
#         end = begin + normal.scalar_mult(t)
#         b = surface(end)
#
#     # now we have the good direction, let's search for a cross
#     #t = t*10
#     # end = begin + normal.scalar_mult(t)
#     # b = surface(end)
#     while a*b > 0:
#         if abs(b) < abs(a):
#             begin = end
#             a = surface(begin)
#             end = begin + normal.scalar_mult(t)
#             b = surface(end)
#             print(a, b, t, "menor")
#         else:
#             #we are in an undesired location
#             if abs(t) < 0.001:
#                 print('saindo')
#                 exit()
#             end = (begin + end)/2
#             b = surface(end)
#             t = t/2
#             print(a, b, t, "maior")
#     print("VIVA!!!!!!!!!!!!!!!!!!")
#
#     middle = (begin + end)/2
#     c = surface(middle)
#     while(abs(c) > projection_aproximation_tolerance):
#         #print(c)
#         if c*a < 0:
#             end = middle
#         else:
#             begin = middle
#         a = surface(begin)
#         middle = (end + begin)/2
#         c = surface(middle)
#
#     return middle



def mesh_refinement(triangles, vertices, indices_map, implicit_function, gradient, threshold = 0.01):
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
    for t in triangles:
        stack = [t]
        while(stack):
            # normalization to put over the unit sphere
            current = stack.pop()
            centroid = (vertices[current.a] + vertices[current.b] + vertices[current.c])/3
            if surface(centroid) > threshold:
                v1 = ((vertices[current.a] + vertices[current.b])/2)
                v2 = ((vertices[current.b] + vertices[current.c])/2)
                v3 = ((vertices[current.a] + vertices[current.c])/2)

                n1 = gradient(v1).normalize()
                n2 = gradient(v2).normalize()
                n3 = gradient(v3).normalize()

                #projection
                v1 = projection_on_surface(v1, n1, implicit_function)
                v2 = projection_on_surface(v2, n2, implicit_function)
                v3 = projection_on_surface(v3, n3, implicit_function)

                stack.append(Triangle( current.a, index_of(v1), index_of(v3)))
                stack.append(Triangle( current.b, index_of(v2), index_of(v1)))
                stack.append(Triangle( current.c, index_of(v3), index_of(v2)))
                stack.append(Triangle( index_of(v1), index_of(v2), index_of(v3)))
            else:
                refined_mesh.append(current)

    triangles = refined_mesh

    write_OFF(os.path.join(data_folder, "surface_adaptative" + str(threshold) + ".off"), vertices, triangles)
    return refined_mesh, vertices, indices_map


# def main(input_file, num_iterations=1):
#     vertices, triangles = read_OFF(input_file)
#     vertices = [v.normalize() for v in vertices]
#     indices_map = { str(v.x) + str(v.y) + str(v.z) : i for i in range(len(vertices)) for v in vertices }
#     mesh_refinement_into_sphere(triangles, vertices, indices_map, num_iterations)

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
    vertices = [projection_on_surface(v, surface_gradient(v).normalize(), surface) for v in vertices]
    #print([surface(v) for v in vertices])

    base_name = input_file[:-4]
    write_OFF(os.path.join(base_name + "_0.off"), vertices, triangles)
    indices_map = { str(v.x) + str(v.y) + str(v.z) : i for i in range(len(vertices)) for v in vertices }
    #mesh_refinement_into_sphere(triangles, vertices, indices_map, num_iterations)
    # triangles_normals = [compute_normal([ vertices[t.a, vertices[t.b], vertices[t.c]) for t in triangles]
    # vertex_normals = []
    # for i in range(len(vertices)):
    #     triangle_numbers = [j for j in range(len(triangles)) if vertices[i] in [vertices[indices[j].a], vertices[indices[j].b], vertices[indices[j].c]]]
    #     vertex_normals.append(compute_vertex_normals([ normals[m] for m in triangle_numbers ] ))


    mesh_refinement(triangles, vertices, indices_map, surface, surface_gradient)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) ... The output is written to a txt file
    ''')
    parser.add_argument('input_file', help='input list of points with X Y per line (format: txt)')
    args = parser.parse_args()

    main(os.path.join(data_folder, args.input_file))
