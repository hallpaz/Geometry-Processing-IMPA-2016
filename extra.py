import argparse
import drawing
import os
import triangle

from DataStructures import Vertex, Triangle, BoundingBox
from data_structures import Node, Color


data_folder = "data"
images_folder = os.path.join("images", "boundary")

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

        neighbors = [index for index in range(num_of_triangles) if len(set(node.data).intersection(triangles[index])) == 2]

        # puts neighbors in correct order (opposite side of vertex)
        node.neighbors = [EMPTY, EMPTY, EMPTY]
        for i in neighbors:
            for v in range(3):
                if node.data[v] not in triangles[i]:
                    node.neighbors[v] = i
        # sort neighbors
        node.costs = [1.0 for i in node.neighbors]

    return graph


def traverse(graph:list, source, pivot, boundary):

    node = graph[source]
    if node.color is Color.white:
        neighbors = [graph[n] for n in node.neighbors if n is not EMPTY]
        for neighbor in neighbors:
            if pivot in neighbor.data:
                # new pivot
                if EMPTY in neighbor.neighbors:
                    i = neighbor.neighbors.index(EMPTY)
                    possible_indices = [j for j in range(3) if j != i]
                    pivot = possible_indices[0] if possible_indices[0] != pivot else possible_indices[1]
                    #update pivot
                boundary.append(graph[neighbor.id])
                traverse(graph, neighbor.id, pivot, boundary)
        node.color = Color.black

def traverse_boundary(graph:list, source, pivot = None):
    """Runs a BFS and compute graph size during the process"""

    if pivot is None:
        interior_vertex = graph[source].neighbors.index(EMPTY)
        possible_vertices = [i for i in range(3) if i != interior_vertex]
        pivot = possible_vertices[0]

    for node in graph:
        node.color = Color.white

    boundary = []
    traverse(graph, source, pivot, boundary)

    return boundary

def draw_boundary(filename: str):
    vertices, indices = read_OFF(os.path.join(data_folder, filename))
    graph = setup_graph(vertices, indices)

    pioneer = None
    for node in graph:
        print(node.neighbors)
        if EMPTY in node.neighbors:
            pioneer = node.id
            break
    boundary = traverse_boundary(graph, pioneer)

    incremental = []
    for i in range(len(boundary)):
        t = indices[boundary[i]]
        incremental.append([t.a, t.b], [t.b, t.c], [t.a, t.c])
        ax = plt.axes()
        drawing.plot_and_save(os.path.join(images_folder, os.path.join("boundary", filename[:-4]+ str(i) +"_boundary_triangles.eps")),
                True, ax, vertices=points, segments=incremental)

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





def mesh_refinement(triangles, vertices, indices_map, implicit_function, gradient, num_iterations = 1):
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

        write_OFF(os.path.join(data_folder, "surface_" + str(3) + ".off"), vertices, triangles)
    return refined_mesh, vertices, indices_map


def projection_on_surface(vertex, normal, surface):
    t = 0.01
    begin = vertex
    end = vertex + normal.scalar_mult(2)
    a = surface(begin)
    b = surface(end)
    i = 0
    print('iniciou: ' + str(t))
    while a*b > 0:
        print(a, b)
        if abs(b) > abs(a):
            #wrong direction
            t = -t
            #print("inverteu: " + str(t))
            end = begin + normal.scalar_mult(t)
            b = surface(end)
        else:
            if t > 0:
                t += 0.1

            #print("dobrou: " + str(t))
            begin = end
            end = vertex + normal.scalar_mult(t)
            a = surface(begin)
            b = surface(end)
        i += 1
        if i > 30:
            #exit()
            if abs(surface(begin)) < abs(surface(vertex)):
                return begin
            elif abs(surface(end)) < abs(surface(vertex)):
                return end
            else:
                return vertex

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
