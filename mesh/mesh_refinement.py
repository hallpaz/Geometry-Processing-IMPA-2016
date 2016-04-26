import argparse
from collections import deque
from DataStructures import Triangle, Vertex
from model_operations import write_OFF, read_OFF
import os

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

def main(input_file, num_iterations=1):
    vertices, triangles = read_OFF(input_file)
    vertices = [v.normalize() for v in vertices]
    indices_map = { str(v.x) + str(v.y) + str(v.z) : i for i in range(len(vertices)) for v in vertices }
    mesh_refinement_into_sphere(triangles, vertices, indices_map, num_iterations)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) ... The output is written to a txt file
    ''')
    parser.add_argument('input_file', help='input list of points with X Y per line (format: txt)')
    args = parser.parse_args()

    main(os.path.join(data_folder, args.input_file))
