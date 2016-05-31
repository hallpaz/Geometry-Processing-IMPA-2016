import argparse
import numpy as np
import os
from geometrypack.data_structures import MinPriorityQueue
from geometrypack.dataio import read_OFF, write_PLY, write_OFF
from geometrypack.algorithms import distance_vertex_plane, planeFit, normalize
from geometrypack.drawing import colorize
from geometrypack.meshes import compute_neighborhood
from triangle import delaunay


results_folder = os.path.join("results", "simplification")
data_folder = os.path.join(os.path.join("data", "meshes"), "simplification")

BIG_ERROR = -1

def compute_error(vertex, vertex_neighborhood, max_error = BIG_ERROR):
    #print(vertex_neighborhood)
    if len(vertex_neighborhood) <= 2:
        return max_error
    planepoint, planenormal = planeFit(vertex_neighborhood)
    d = distance_vertex_plane(vertex, planepoint, planenormal)
    if d < 0.00001:
        d = 0.0
    return d


def simplification_heatcolor(vertices, indices):
    neighborhood = compute_neighborhood(vertices, indices)
    errors = [compute_error(vertices[vertex_index], [vertices[i] for i in neighborhood[vertex_index]])
            for vertex_index in range(len(vertices))]

    max_value = max(errors)
    print("erro maximo", max_value)
    for i in range(len(errors)):
        if errors[i] == BIG_ERROR:
            errors[i] = max_value
    min_value = min(errors)
    print("erro minimo", min_value)
    if max_value - min_value > 0.1:
        errors = [(e - min_value)/(max_value-min_value) for e in errors]
    colors = [colorize(e) for e in errors]
    return colors

def simplify_mesh(vertices, indices, threshold = None, reduction_ratio = None):
    # must use tuple because I need a hasheable type
    simplified_vertices = [(v[0], v[1], v[2]) for v in vertices ]
    simplified_indices = [[t[0], t[1], t[2]] for t in indices ]

    if threshold is None and reduction_ratio is None:
        print("You must provide either a threshold or a reduction ratio")
        return [], []

    queue = MinPriorityQueue()
    neighborhood = compute_neighborhood(vertices, indices)
    errors = [compute_error(vertices[vertex_index], [vertices[i] for i in neighborhood[vertex_index]])
            for vertex_index in range(len(vertices))]
    max_error = max(errors)

    for index in range(len(simplified_vertices)):
        queue.add_point(simplified_vertices[index], errors[index] if errors[index] != BIG_ERROR else max_error)
    errors = []
        # queue.add_point(simplified_vertices[index], compute_error(simplified_vertices[index],
        #     [simplified_vertices[i] for i in neighborhood[index]]))
    # this list keeps the indices of all vertices that must be removed at the end
    tobe_removed = []
    candidate = queue.pop_point()
    if reduction_ratio is not None:
        beginsize = len(simplified_vertices)
        desired_size = beginsize*reduction_ratio
        #desired_size = reduction_ratio
        while beginsize - len(tobe_removed) > desired_size:
            # find the index of the vertex in the original sequence
            index = simplified_vertices.index(candidate[-1])
            # marks the vertex to be removed from the sequence
            tobe_removed.append(index)
            # take a reference to the list of the current vertex neighbor's
            current_neighbors = list(neighborhood[index])
            # we always take the first neighbor as pivot
            pivot = current_neighbors[0]
            npoints = [vertices[n] for n in current_neighbors]
            if len(npoints) > 2:
                # finds the plane that best fits (LSq) the neighbors vertices
                planecentroid, planenormal = planeFit(npoints)
                # normalizes the plane normal vector
                planenormal = normalize(planenormal)
                # finds the projection of the pivot vertex onto the best fitting plane
                pivot_proj = vertices[pivot] - np.dot(vertices[pivot] - planecentroid, planenormal) * planenormal
                # we'll use the direction from the "plane centroid" to the pivot's projection
                # as an axis of our frame of reference in this plane. We need 2D coordinates
                # to use triangle's delaunay triangulation function
                axis1 = normalize(pivot_proj - planecentroid)
                # the other axis is computed by cross product
                axis2 = np.cross(axis1, planenormal)
                # we compute the 2D coordinates for all neighbors in this new frame
                projected_points = [[np.dot(axis1, p - planecentroid),
                    np.dot(axis2, p - planecentroid)] for p in npoints]
                # delaunay triangulation of the hole left
                tri = delaunay(projected_points)
                for t in tri:
                    # I still can't characterize when I must invert orientation :/
                    simplified_indices.extend([[current_neighbors[t[0]], current_neighbors[t[2]], current_neighbors[t[1]]]])
                    # we update the list of indices with new triangles (we didnt actually removed the old ones)
                    # we also update the list of neighbors with new neighbors
                    for i in range(3):
                        neighborhood[current_neighbors[t[i]]].add(current_neighbors[t[(i+1)%3]])
                        neighborhood[current_neighbors[t[i]]].add(current_neighbors[t[(i+2)%3]])
            # we do remove the point from the neighbor's list because it changes
            # the value of the error
            for n in current_neighbors:
                # removes vertex from neighbor's list
                neighborhood[n].remove(index)

            # now we update the error cost at each neighbor
            for n in current_neighbors:
                queue.add_point(simplified_vertices[n], compute_error(vertices[n],
                        [simplified_vertices[i] for i in neighborhood[n]], max_error))

            # pops the new candidate
            candidate = queue.pop_point()
    else: #won't bother with optimization now
        print("should't be here!")
        #while(candidate[0] < threshold):

    print("will clean vertices and indices")
    # set to compute intersection
    tobe_removed = set(tobe_removed)
    simplified_indices = [t for t in simplified_indices if len(tobe_removed.intersection(t)) == 0 ]
    simplified_vertices = [simplified_vertices[i] for i in range(len(simplified_vertices)) if i not in tobe_removed]
    print("will update indices")
    tobe_removed = list(tobe_removed)
    tobe_removed.sort()
    # removes from backwards to avoid conflict
    tobe_removed.reverse()
    for removed_index in tobe_removed:
        for t in simplified_indices:
            for i in range(3):
                if t[i] > removed_index:
                    t[i] = t[i] -1

    return simplified_vertices, simplified_indices

def mesh_simplification(filepath, threshold= 1):
    #extract the filename withou extension
    filename = os.path.splitext(os.path.basename(filepath))[0]
    print(filename)
    vertices, indices = read_OFF(filepath)
    normals = [[0, 0, 0] for i in vertices]
    colors = simplification_heatcolor(vertices, indices)
    # creates dir if it does not exist
    folder = os.path.join(results_folder, filename)
    os.makedirs(folder, exist_ok=True)
    os.chmod(folder, 0o777)

    heatfile = os.path.join(folder, "heat_" + filename + ".ply")
    write_PLY(heatfile, vertices, indices, normals, colors)
    print("Simplification heatmap written to {}".format(heatfile))
    for ratio in range(90, 0, -10):
        newvertices, newindices = simplify_mesh(vertices, indices, reduction_ratio=ratio/100)
        normals = [[0, 0, 0] for i in vertices]
        colors = simplification_heatcolor(newvertices, newindices)
        write_PLY(os.path.join(folder, filename + "_p{}.ply".format(ratio)), newvertices, newindices, normals, colors)

def rec_all_data(folder):
    filenames = [fname for fname in os.listdir(folder) if fname.endswith('.off')]
    for fname in filenames:
        mesh_simplification(os.path.join(data_folder, fname))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) ... The output is written to a txt file
    ''')
    parser.add_argument('input_file', help='input list of points with X Y per line (format: txt)')
    args = parser.parse_args()
    if args.input_file == "all":
        rec_all_data(data_folder)
    else:
        mesh_simplification(args.input_file)
