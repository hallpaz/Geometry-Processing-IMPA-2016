import argparse
import drawing
import os
import triangle
import matplotlib.pyplot as plt
import numpy as np

from DataStructures import Vertex, Triangle, BoundingBox
from data_structures import Node, Color


data_folder = "data"
images_folder = "images"

BIG = 7777777
EMPTY = -1

def read_OFF(off_file):
    vertexBuffer = []
    indexBuffer = []
    with open(off_file, "r") as modelfile:
        first = modelfile.readline().strip()
        if first != "OFF":
            raise(Exception("not a valid OFF file ({})".format(first)))

        parameters = modelfile.readline().strip().split()
        min_distance = BIG
        minx = BIG
        miny = BIG
        minz = BIG
        maxx = -BIG
        maxy = -BIG
        maxz = -BIG

        if len(parameters) < 2:
            raise(Exception("OFF file has invalid number of parameters"))

        for i in range(int(parameters[0])):
            coordinates = modelfile.readline().split()
            X, Y, Z = float(coordinates[0]), float(coordinates[1]), float(coordinates[2])
            vertexBuffer.append(Vertex(X, Y, Z))
            if X < minx:
                minx = X
            if X > maxx:
                maxx = X
            if Y < miny:
                miny = Y
            if Y > maxy:
                maxy = Y
            if Z < minz:
                minz = Z
            if Z > maxz:
                maxz = Z

        for i in range(int(parameters[1])):
            indices = modelfile.readline().split()
            indexBuffer.append(Triangle(int(indices[1]), int(indices[2]), int(indices[3])))

        bounding_box = BoundingBox(minx, miny, minz-1, maxx, maxy, maxz+1)

    return vertexBuffer, indexBuffer#, bounding_box

def soup_to_mesh(filename: str, output_file: str, dimension = 2):
    vertices = []
    indices = []
    with open(filename, 'r') as soup_file:
        soup_lines = soup_file.readlines()
        for line in soup_lines:
            coordinates = line.split()
            t = []
            for i in range(3):
                v = Vertex(float(coordinates[0+dimension*i]), float(coordinates[1+dimension*i]), 0 if dimension == 2 else float(coordinates[2+dimension*i]))
                if v in vertices:
                    # repeated vertex, reetrieve index
                    index = vertices.index(v)
                    t.append(index)
                else:
                    # new vertex, new index
                    t.append(len(vertices))
                    vertices.append(v)
            indices.append(Triangle(t[0], t[1], t[2]))


    vertices = [str(v) for v in vertices]
    indices = [str(i) for i in indices]

    with open(output_file, 'w') as meshfile:
        meshfile.write(
        '''OFF
        %d %d 0
        %s%s
        '''%(len(vertices),len(indices), "".join(vertices), "".join(indices)))


### Complete this ###
# def share_edge(edge: set, triangle: list):
#
#     if edge.intersection(triangle) == edge

# def traverse_boundary(vertices, indices):
#     boundary_triangles = [t for t in indices if ]
#
#     #graph initialization
#     indices = [set([t.a, t.b, t.c]) for t in indices]
#     num_of_triangles = len(indices)
#     graph = [Node(i) for i in range(num_of_triangles)]
#     for node in graph:
#         current_vertices  = indices[node.id]
#         node.neighbors = [index for index in range(num_of_triangles) if len(set(current_vertices).intersection(indices[index])) == 2]
#         node.costs = [1.0 for i in node.neighbors]
#
#     connected_one = []
#     def append_to(node, connected_one:list):
#         connected_one.append(node)
#
#     BFS(graph, 1, append_to, connected_one)
#     boundary_triangles = [filtered_triangles[node.id] for node in connected_one if len(node.neighbors) < 3]

def setup_graph(vertices, indices):
    #graph initialization
    triangles = [[t.a, t.b, t.c] for t in indices]
    num_of_triangles = len(triangles)
    graph = [Node(i) for i in range(num_of_triangles)]
    for n in graph:
        n.data = triangles[n.id]
    for node in graph:
        # current_vertices  = triangles[node.id]
        node.neighbors = [graph[index] for index in range(num_of_triangles) if len(set(node.data).intersection(triangles[index])) == 2]
        node.costs = [1.0 for i in node.neighbors]

    print(len(graph))
    return graph


def isOnBoundary(node):
    if len(node.neighbors) < 3:
        return True
    return False

def traverse_boundary(graph:list, source, pivot = None):
    """Runs a BFS and compute graph size during the process"""

    if pivot is None:
        pivot = source.data[0]
        if pivot in source.neighbors[0].data and pivot in source.neighbors[1].data:
            pivot = source.data[1]

    for node in graph:
        node.color = Color.white

    boundary = []
    stack = []
    stack.append(source)
    boundary.append(source)
    while(stack):
        node = stack.pop()
        if node.color is not Color.black:
            for neighbor in node.neighbors:
                if neighbor.color is not Color.black:
                    if pivot in neighbor.data:
                        print("in", pivot, neighbor.data)
                        boundary.append(neighbor)
                        stack.append(neighbor)
                        if isOnBoundary(neighbor):
                            # neighbor is on the boundary
                            pivot = [newpivot for newpivot in neighbor.data if newpivot not in node.data][0]
                            #update pivot
                        break
                    else:
                        print(pivot, neighbor.data, node.data)

            node.color = Color.black
        else:
            print("black node")

    return boundary

def draw_boundary(filename: str):
    vertices, indices = read_OFF(os.path.join(data_folder, filename))
    graph = setup_graph(vertices, indices)
    points = np.array([(v.x, v.y) for v in vertices])

    pioneer = None
    for node in graph:
        #print(node.neighbors)
        if isOnBoundary(node):
            pioneer = node
            break
    boundary = traverse_boundary(graph, pioneer)
    print(len(boundary))
    incremental = []
    all_triangles = [b.data for b in boundary]
    for i in range(len(boundary)):
        t = boundary[i].data
        incremental.append((t[0], t[1]))
        incremental.append((t[1], t[2]))
        incremental.append((t[0], t[2]))
        ax = plt.axes()
        drawing.plot_and_save(os.path.join(images_folder, os.path.join("boundary", filename[:-4]+ str(i) +"_boundary_triangles.png")),
                True, ax, vertices=points, segments=incremental, triangles=all_triangles)

# def face_association_from_OFF(filename: str):
#     vertices, indices = read_OFF(filename)
#     triangle_connectivity = []
#
#     indices = [set([t.a, t.b, t.c]) for t in indices]
#     for t in indices:
#         neighbors = [j for j in indices if t.intersection(j) == 2]
#         connectivity.append( t.index() )
#         neighbors.append()
# #####################

def main(input_file, dimension = 2):
    dot_index = input_file.rfind(".")
    output_file = (input_file[:dot_index] if dot_index > 0 else input_file.strip()) + ".off"
    soup_to_mesh(input_file, output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) ... The output is written to a txt file
    ''')
    parser.add_argument('input_file', help='input list of points with X Y per line (format: txt)')
    args = parser.parse_args()
    if args.input_file == "all":
        rec_all_data(data_folder)
    elif args.input_file == "soup.off":
        draw_boundary(args.input_file)

    else:
        main(os.path.join(data_folder, args.input_file))
