import argparse
import drawing
import os
import matplotlib.pyplot as plt
import numpy as np
import triangle


data_folder = "data"
images_folder = os.path.join("images", "nncrust")

def distance(points, edge):
    p1 = points[edge[0]]
    p2 = points[edge[1]]
    #print(p1, p2)
    return np.linalg.norm(p1-p2)

def cos_angle(pivot, q, s):
    pq = q - pivot
    ps = s - pivot
    return np.dot(pq, ps)/(np.linalg.norm(pq)*np.linalg.norm(ps))

def edges_from_triangle(triangle):
    return [ (triangle[0], triangle[1]), (triangle[1], triangle[2]), (triangle[2], triangle[0]) ]

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

        NNcrust(points, output_file)

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
