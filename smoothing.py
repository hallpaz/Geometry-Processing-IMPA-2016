import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
from geometrypack.curves import simple_mean_curve_smoothing
from geometrypack.dataio import read_points, write_points
from geometrypack import drawing

from doideira import draw


data_folder = "data"
images_folder = "images"

SIMPLE = 1
LAPLACIAN = 2

def smooth_curve(filename, method):
    points = read_points(os.path.join("data/points", filename))

    print(points)
    print(len(points))
    if method == SIMPLE:
        smoothed_points = simple_mean_curve_smoothing(points)
        write_points(smoothed_points, os.path.join("results/points", "smoothed_" + filename))
        teste = [np.array([2, 3]), np.array([8, 12]), np.array([ 14, 3])]
        draw(teste, "teste.eps")
        #drawing.draw_segments(smoothed_points, "smoothed0_" + filename + ".eps")
        # ax = plt.axes()
        # drawing.plot_and_save(os.path.join("results", "smoothed0_" + filename + ".png"), False,
        #  ax, vertices=smoothed_points, vertices_color="g",
        #  segments=[np.array([i, (i+1)%len(smoothed_points)]) for i in range(len(smoothed_points))])

    else:
        print("XABU no metodo")

def smooth_mesh(filename, method):
    print("Xabu na malha")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) or a 3D mesh (vertices and indices)
    and outputs a smoothed version of the input''')

    parser.add_argument('input_file', help='input list of points with X Y per line (format: txt) or a .off file')
    parser.add_argument('method', help='smoothing method: 0 for SIMPLE, 1 for LAPLACIAN')
    args = parser.parse_args()

    filename = args.input_file
    if filename.endswith(".txt"):
        smooth_curve(filename, int(args.method))
    elif filename.endswith(".off"):
        smooth_mesh(filename, int(args.method))
    else:
        raise Exception("File must be a .txt or .off")
