import numpy as np
from collections import deque
from .data_structures import Color, Node
from scipy.linalg import expm3, norm


def rotation_matrix(axis, theta):
    return expm3(np.cross(np.eye(3), axis/norm(axis)*theta))
# angle in radians
def rotate(vertices, axis, angle):
    m = rotation_matrix(axis, angle)
    return [np.dot(m,v) for v in vertices]

def compute_neighborhood(vertices, indices):
    # neighborhood[i] é uma lista que contém os índices de todos os triângulos
    # que contém o vértice i
    neighborhood = [ set([]) for i in vertices]
    for index in range(len(indices)):
        # t recebe os 3 indices do triangulo da posição 'index'
        t = indices[index]
        for i in range(3):
            neighborhood[t[i]].add(t[(i+1)%3])
            neighborhood[t[i]].add(t[(i+2)%3])

    return neighborhood

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
    return Vertex(2*x*(2*x**2 - 5), 2*y*(2*y**2 - 5), 2*z*(2*z**2 - 5) )


def compute_normal(triangle_vertices):
    v1 = triangle_vertices[1] - triangle_vertices[0]
    v2 = triangle_vertices[2] - triangle_vertices[1]
    normal = Vertex.cross(v1, v2).normalize()
    return normal

def projection_on_surface(vertex, normal, surface):
    t = 1.0
    begin = vertex
    end = vertex + normal.scalar_mult(2)
    a = surface(begin)
    b = surface(end)
    while a*b > 0:
        print(a, b)
        if abs(b) > abs(a):
            #wrong direction
            t = -t
            end = vertex + t*normal
            b = surface(end)
        else:
            t = 2*t
            begin = end
            end = vertex + t*normal
            a = surface(begin)
            b = surface(end)

    middle = (begin + end)/2
    c = surface(middle)
    while(abs(c) > 0.0001):
        print(c)
        if c*a < 0:
            end = middle
        else:
            begin = middle
        a = surface(begin)
        c = surface(middle)
        middle = (end + begin)/2

    return middle



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

        write_OFF(os.path.join(data_folder, "impliciti_project" + str(i) + ".off"), vertices, triangles)
    return refined_mesh, vertices, indices_map


########################### Fuctions for Meshes #############################

def simple_mesh_smoothing(points, graph, anchor_indices = []):
    if len(points) > 3:
        # initialization
        for node in graph:
            if node.id in anchor_indices:
                node.color = Color.Black
            else:
                node.color = Color.White
        # initialize a list with the same length as points
        smoothed_points = [None for i in range(len(points))]
        # anchored
        for index in anchor_indices:
            smoothed_points[index] = points[index]

        queue = deque()
        for p in graph:
            if p.color is Color.White:
                queue.append(p)
                while(queue):
                    node = queue.popleft()
                    i = 0
                    newvalue = np.array([0 for i in range(len(node.value))])
                    for n in node.neighbors:
                        i += 1
                        newvalue = newvalue + points[n.id]#n.value
                        if n.color is Color.White:
                            # mark as 'in progress'
                            n.color = Color.Gray
                            queue.append(n)
                    # take the mean
                    newvalue /= i
                    smoothed_points[node.id] = newvalue
                    # mark as 'finished'
                    node.color = Color.Black

        return np.array(smoothed_points)

    return np.array(points)

def setup_graph(vertices, indices):
    #graph initialization
    triangles = [[t[0], t[1], t[2]] for t in indices]
    num_of_triangles = len(triangles)
    graph = [Node(i) for i in range(num_of_triangles)]
    for n in graph:
        n.data = triangles[n.id]
    for node in graph:
        # current_vertices  = triangles[node.id]
        node.neighbors = [graph[index] for index in range(num_of_triangles) if len(set(node.data).intersection(triangles[index])) == 2]
        node.costs = [1.0 for i in node.neighbors]

    #print(len(graph))
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

    boundary = []
    boundary_indices = [pivot]
    stack = []
    stack.append(source)
    boundary.append(source)
    while(stack):
        node = stack.pop()
        if node.color is not Color.Black:
            for neighbor in node.neighbors:
                if neighbor.color is not Color.Black:
                    if pivot in neighbor.data:
                        #print("in", pivot, neighbor.data)
                        boundary.append(neighbor)
                        stack.append(neighbor)
                        if isOnBoundary(neighbor):
                            # neighbor is on the boundary
                            pivot = [newpivot for newpivot in neighbor.data if newpivot not in node.data][0]
                            boundary_indices.append(pivot)
                            #update pivot
                        break
                    else:
                        #print(pivot, neighbor.data, node.data)
                        pass

            node.color = Color.Black
        else:
            print("Black node")

    return boundary, boundary_indices

def compute_boundary_indices(vertices, indices):
    graph = setup_graph(vertices, indices)
    points = np.array([(v[0], v[1]) for v in vertices])
    for node in graph:
        # node.color = Color.White
        node.color = Color.White

    should_continue = False
    boundary_indices = []
    k = 0
    while(True):
        pioneer = None
        for node in graph:
            #print(node.neighbors)
            if isOnBoundary(node) and node.color is Color.White:
                pioneer = node
                should_continue = True
                break
        if should_continue:
            boundary, bindex = traverse_boundary(graph, pioneer)
            boundary_indices.extend(bindex)
        else:
            break
        should_continue = False
        k += 1
    return boundary_indices
