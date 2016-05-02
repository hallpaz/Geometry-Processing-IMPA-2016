import numpy as np
from enum import IntEnum

class Direction(IntEnum):
    """docstring for Direction"""
    East = 0
    South = 1
    West = 2
    North = 3


def pavlids_boundary_tracing(data):
    start_col = 0
    start_row = 0
    found = False
    seen = 0
    for row in range(len(data)):
        for col in range(len(data[0])):
            seen += 1
            if data[row][col]:
                start_col = col
                start_row = row
                found = True
                print(row, col, data[row][col], len(data), len(data[0]))
                break
        if found:
            break

    if not found:
        return None
    current_direction = Direction.East
    boundary = [np.array([start_row, start_col])]

    def p1(direction, i, j):
        if direction is Direction.North:
            return i-1, j-1
        elif direction is Direction.West:
            return i+1, j-1
        elif direction is Direction.South:
            return i+1, j+1
        elif direction is Direction.East:
            return i-1, j+1
        else:
            raise "Invalid Direction"

    def p2(direction, i, j):
        if direction is Direction.North:
            return i-1, j
        elif direction is Direction.West:
            return i, j-1
        elif direction is Direction.South:
            return i+1, j
        elif direction is Direction.East:
            return i, j+1
        else:
            raise "Invalid Direction"

    def p3(direction, i, j):
        if direction is Direction.North:
            return i-1, j+1
        elif direction is Direction.West:
            return i-1, j-1
        elif direction is Direction.South:
            return i+1, j-1
        elif direction is Direction.East:
            return i+1, j+1
        else:
            raise "Invalid Direction"

    rotations = 0
    row, col = start_row, start_col
    while rotations <= 3:
        #print("it", row, col)
        i, j = p1(current_direction, row, col)
        try:
            next_point = data[i][j]
        except Exception as e:
            print("row: ", row, "col: ", col)
            print(i, j, "direction: ", current_direction, "error: ", e)
            raise
        if next_point:
            if i == start_row and j == start_col:
                break
            boundary.append(np.array([i*10, j*10]))
            row, col = i, j
            current_direction = Direction((current_direction - 1)%4)
            rotations = 0
        else:
            i, j = p2(current_direction, row, col)
            try:
                next_point = data[i][j]
            except Exception as e:
                print("row: ", row, "col: ", col)
                print(i, j, "direction: ", current_direction, "error: ", e)
                raise
            if next_point:
                if i == start_row and j == start_col:
                    break
                boundary.append(np.array([i*10, j*10]))
                row, col = i, j
                rotations = 0
                #current_direction = Direction((current_direction)%4)
            else:
                i, j = p3(current_direction, row, col)

                try:
                    next_point = data[i][j]
                except Exception as e:
                    print("row: ", row, "col: ", col)
                    print(i, j, "direction: ", current_direction, "error: ", e)
                    raise
                if next_point:
                    if i == start_row and j == start_col:
                        break
                    boundary.append(np.array([i*10, j*10]))
                    row, col = i, j
                    rotations = 0
                else:
                    current_direction = Direction((current_direction+1)%4)
                    rotations += 1
    print("rotations: ", rotations)
    print("row: ", row, "col: ", col)
    return boundary

def crust(points: list, filename = None)-> list:
    points = np.array(points)

    ax = plt.axes()
    drawing.plot_and_save(os.path.join(images_folder, filename+"_original.png"), False, ax, vertices=points, vertices_color="g")
    lim = ax.axis()

    vertices, edges, ray_origins, ray_directions = triangle.voronoi(points)
    try:
        drawing.plot_and_save(os.path.join(images_folder, os.path.join("voronoi", filename+"_voronoi.png")),
         False, ax, vertices=vertices, edges=edges, ray_origins=ray_origins, ray_directions=ray_directions)
        ax.axis(lim)
    except Exception as e:
        print("Could not draw voronoi diagram. Please, check if you have the correct directory")

    extended_points = np.concatenate((points, vertices))

    triangles = triangle.delaunay(extended_points)
    delaunay_edges = []
    for t in triangles:
        delaunay_edges.extend(edges_from_triangle(t))

    rec_edges = [e for e in delaunay_edges if ((e[0] < len(points)) and (e[1] < len(points)))]
    try:
        drawing.plot(ax, vertices=extended_points, edges=delaunay_edges, segments=rec_edges)
        drawing.plot_and_save(os.path.join(images_folder, os.path.join("delaunay",filename+"delaunay.png")),
            True, ax, vertices=points, vertices_color="b")
        ax.axis(lim)
    except Exception as e:
        print("Could not draw delaunay diagram. Please, check if you have the correct directory")


    #final reconstruction
    ax = plt.axes()
    drawing.plot_and_save(os.path.join(images_folder, filename+"_reconstruction.png"),
                            True, ax, vertices=points, segments=rec_edges, draw_vertices=False)
    # print(points)
    # print(rec_edges)
    ordered_points = []

    return ordered_points

def NNcrust(points: list, filename = None)-> list:
    points = np.array(points)

    ax = plt.axes()
    drawing.plot_and_save(os.path.join(images_folder, filename+"_original.png"), False, ax, vertices=points, vertices_color="g")
    lim = ax.axis()

    triangles = triangle.delaunay(points)
    delaunay_edges = []
    for t in triangles:
        delaunay_edges.extend(edges_from_triangle(t))

    rec_edges = set([])
    for i in range(len(points)):
        edges_with_p = [edge for edge in delaunay_edges if i in edge]
        if edges_with_p:
            costs = [distance(points, edge) for edge in edges_with_p]
            min_index = costs.index(min(costs))
            pq = edges_with_p[min_index]
            rec_edges.add(pq)

        ps_candidates = [edge for edge in edges_with_p if
            cos_angle(points[i], points[pq[0] if pq[0] != i else pq[1]], points[edge[0] if edge[0] != i else edge[1]]) < 0]
        if ps_candidates:
            costs = [distance(points, edge) for edge in ps_candidates]
            min_index = costs.index(min(costs))
            ps = ps_candidates[min_index]
            rec_edges.add(ps)

    # rec_edges = [e for e in delaunay_edges if ((e[0] < len(points)) and (e[1] < len(points)))]
    try:
        drawing.plot(ax, vertices=points, edges=delaunay_edges, segments=rec_edges)
        drawing.plot_and_save(os.path.join(images_folder, os.path.join("delaunay",filename+"delaunay.png")),
            True, ax, vertices=points, vertices_color="b")
        ax.axis(lim)
    except Exception as e:
        print("Could not draw delaunay diagram. Please, check if you have the correct directory")


    #final reconstruction
    ax = plt.axes()
    drawing.plot_and_save(os.path.join(images_folder, filename+"_reconstruction.png"),
                            True, ax, vertices=points, segments=rec_edges, draw_vertices=False)
    # print(points)
    # print(rec_edges)
    ordered_points = []

    return ordered_points

def peucker_reduction(points, error_tolerance):
    if not points:
        return None
    dmax = 0
    index = 0
    fitting_line = Segment(points[0], points[-1])
    i = 1
    for point in points[1:-1]:
        d = fitting_line.distance_to_point(point)
        if d > dmax:
            index = i
            dmax = d
        i += 1

    if dmax > error_tolerance:
        res1 = peucker_reduction(points[0: index+1], error_tolerance)
        res2 = peucker_reduction(points[index:], error_tolerance)
        result = res1[:-1] + res2
    else:
        result = [points[0], points[-1]]

    return result


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
    drawing.plot_and_save(os.path.join(images_folder, os.path.join("delaunay",filename+"delaunay.eps")),
        True, ax, vertices=points, vertices_color="b")
    ax.axis(lim)


    triangles_length = [perimeter(points, t) for t in deltriangles]
    mean_length = sum(triangles_length)/len(triangles_length)
    std_length = np.std(triangles_length)
    filtered_triangles = [deltriangles[i] for i in range(len(deltriangles)) if triangles_length[i] < (mean_length + 1*std_length)]

    ax = plt.axes()
    drawing.plot_and_save(os.path.join(images_folder, os.path.join("boundary", filename+"_region1std.eps")),
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
    drawing.plot_and_save(os.path.join(images_folder, os.path.join("boundary", filename+"_boundary_triangles1std.eps")),
                True, ax, triangles = boundary_triangles, vertices=points, segments=boundary_segments)

    ax = plt.axes()
    drawing.plot_and_save(os.path.join(images_folder, os.path.join("boundary", filename+"_boundary_triangles1stdstd.eps")),
                True, ax, vertices=points, segments=boundary_segments)


########################### Fuctions for Curves #############################
def simple_mean_curve_smoothing(points):
    if len(points) > 3:
        # python deals with negative indices as expected
        return np.array([points[i-1]/2 + points[(i+1)%len(points)]/2
                    for i in range(len(points))])

    return points

def simple_mean_smoothing_curve_with_anchor_points(points, anchor_indices):
    if len(points) > 3:
        # python deals with negative indices as expected
        return np.array([points[i-1]/2 + points[(i+1)%len(points)]/2
                for i in range(len(points))
                if i not in anchor_indices])

    return points


def laplacian_curve_smoothing(points):
    if len(points) > 3:
        # python deals with negative indices as expected
        return np.array([points[i-1]/4 + points[i]/2 + points[(i+1)%len(points)]/4
                for i in range(len(points))])

    return points


def laplacian_smoothing_curve_with_anchor_points(points, anchor_indices):
    if len(points) > 3:
        # python deals with negative indices as expected
        return np.array([points[i-1]/4 + points[i]/2 + points[(i+1)%len(points)]/4
                for i in range(len(points))
                if i not in anchor_indices])

    return points
