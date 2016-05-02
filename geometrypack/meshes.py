

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
