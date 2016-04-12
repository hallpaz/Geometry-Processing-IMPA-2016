import argparse
import drawing
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import triangle

from crust import edges_from_triangle
from data_structures import Node
from graphs import BFS

data_folder = os.path.join("data", "region")
images_folder = os.path.join("images", "region")

def perimeter(points, triangle):
    p = 0
    p += np.linalg.norm(points[triangle[0]]-points[triangle[1]])
    p += np.linalg.norm(points[triangle[1]]-points[triangle[2]])
    p += np.linalg.norm(points[triangle[2]]-points[triangle[0]])
    return p

def heuristic_reconstruction(points: list, filename = None):
    points = np.array(points)

    ax = plt.axes()
    #drawing.plot_and_save(os.path.join(images_folder, filename+"_original.png"), False, ax, vertices=points, vertices_color="g")
    lim = ax.axis()

    deltriangles = triangle.delaunay(points)
    delaunay_edges = []
    for t in deltriangles:
        delaunay_edges.extend(edges_from_triangle(t))

    #try:
    drawing.plot(ax, vertices=points, segments=delaunay_edges)
    drawing.plot_and_save(os.path.join(images_folder, os.path.join("delaunay",filename+"delaunay.png")),
        True, ax, vertices=points, vertices_color="b")
    ax.axis(lim)


    triangles_length = [perimeter(points, t) for t in deltriangles]
    mean_length = sum(triangles_length)/len(triangles_length)
    std_length = np.std(triangles_length)
    filtered_triangles = [deltriangles[i] for i in range(len(deltriangles)) if triangles_length[i] < (mean_length + 0.2*std_length)]

    ax = plt.axes()
    drawing.plot_and_save(os.path.join(images_folder, os.path.join("boundary", filename+"_region.png")),
                True, ax, triangles = filtered_triangles, vertices=points, draw_vertices=False)

    #graph initialization
    num_of_triangles = len(filtered_triangles)
    graph = [Node(i) for i in range(num_of_triangles)]
    for node in graph:
        vertices = filtered_triangles[node.id]
        node.neighbors = [index for index in range(num_of_triangles) if len(set(vertices).intersection(filtered_triangles[index])) == 2]
        node.costs = [1.0 for i in node.neighbors]

    connected_one = []
    def append_to(node, connected_one:list):
        connected_one.append(node)

    BFS(graph, 1, append_to, connected_one)


    boundary_triangles = [filtered_triangles[node.id] for node in connected_one if len(node.neighbors) < 3]

    filtered_segments = []
    for t in filtered_triangles:
        filtered_segments.extend(edges_from_triangle(t))
    boundary_segments = [edge for edge in filtered_segments if filtered_segments.count(edge) == 1]
    # print(boundary_segments)

    ax = plt.axes()
    drawing.plot_and_save(os.path.join(images_folder, os.path.join("boundary", filename+"_boundary_triangles.png")),
                True, ax, triangles = boundary_triangles, vertices=points, segments=boundary_segments)


    #except Exception as e:
    #    print("Could not draw delaunay diagram. Please, check if you have the correct directory")


def main(input_file):
    dot_index = input_file.rfind(".")
    output_file = input_file[:dot_index] if dot_index > 0 else input_file.strip()
    input_file = os.path.join(data_folder, input_file)
    print("processing:", input_file)
    with open(input_file, "r") as myfile:
        points = []
        file_lines = myfile.readlines()
        for line in file_lines:
            content = line.split()
            points.append([float(content[0]), float(content[1])])
        points = np.array(points)

        heuristic_reconstruction(points, output_file)

def rec_all_data(folder):
    filenames = [fname for fname in os.listdir(folder) if fname.endswith('.txt')]
    for fname in filenames:
        main(fname)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) ... The output is written to a txt file
    ''')
    parser.add_argument('input_file', help='input list of points with X Y per line (format: txt)')
    args = parser.parse_args()
    if args.input_file == "all":
        rec_all_data(data_folder)
    else:
        main(args.input_file)
