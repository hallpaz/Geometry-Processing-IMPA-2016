import argparse
import drawing
import os
import matplotlib.pyplot as plt
import numpy as np
import triangle


data_folder = "data"
images_folder = os.path.join("images", "crust")

def edges_from_triangle(triangle):
    return [ (triangle[0], triangle[1]), (triangle[1], triangle[2]), (triangle[2], triangle[0]) ]

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
            points.append([content[0], content[1]])
        points = np.array(points)

        crust(points, output_file)

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
