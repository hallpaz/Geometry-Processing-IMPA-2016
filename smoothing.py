import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
from geometrypack.curves import simple_mean_curve_smoothing, laplacian_curve_smoothing
from geometrypack.dataio import read_points, write_points
from geometrypack import drawing

from doideira import draw
from copy import deepcopy


data_folder = os.path.join("results", "smoothing")
images_folder = os.path.join(data_folder, "images")

SIMPLE = 1
LAPLACIAN = 2

def smooth_curve(filename, method, smooth_level = 1):
    points = read_points(os.path.join("data/points", filename))

    output_filename = filename[:-4]

    if method == SIMPLE:
        smooth_method = simple_mean_curve_smoothing
    else:
        smooth_method = laplacian_curve_smoothing

    smoothed_points = deepcopy(points)
    for level in range(smooth_level):
        smoothed_points = smooth_method(deepcopy(smoothed_points))

        # write_points(smoothed_points, os.path.join(data_folder,
        #             output_filename + "_smoothed{}.txt".format(level)))

        ax = plt.axes()
        drawing.plot_and_save(os.path.join(os.path.join(images_folder, output_filename),
            output_filename + "_{}{}.png".format(smooth_method.__name__, level)), False,
            ax, vertices=points, vertices_color="k", edges=[(i, (i+1)%len(points)) for i in range(len(points))])

        drawing.plot_and_save(os.path.join(os.path.join(images_folder, output_filename),
            output_filename + "_{}{}.png".format(smooth_method.__name__, level)), True,
            ax, vertices=smoothed_points, vertices_color="g",
            segments=[(i, (i+1)%len(smoothed_points)) for i in range(len(smoothed_points))])



def smooth_mesh(filename, method):
    print("Xabu na malha")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) or a 3D mesh (vertices and indices)
    and outputs a smoothed version of the input''')

    parser.add_argument('input_file', help='input list of points with X Y per line (format: txt) or a .off file')
    parser.add_argument('method', help='smoothing method: 0 for SIMPLE, 1 for LAPLACIAN')
    parser.add_argument('iterations', help='number of iterations')
    args = parser.parse_args()

    filename = args.input_file
    if filename.endswith(".txt"):
        smooth_curve(filename, int(args.method), int(args.iterations))
    elif filename.endswith(".off"):
        smooth_mesh(filename, int(args.method), int(args.iterations))
    else:
        raise Exception("File must be a .txt or .off")
