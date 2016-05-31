import argparse
import heapq
import matplotlib.pyplot as plt
import numpy as np
import os


from geometrypack.data_structures import MinPriorityQueue
from geometrypack.dataio import read_points, write_points
from geometrypack.algorithms import distance_point_line
from geometrypack import drawing


results_folder = os.path.join("results", "simplification")
data_folder = os.path.join("data", "points")


def compute_error(points, current, i = -1):
    if i == -1:
        i = points.index(current)
    d = distance_point_line(current, points[i-1], points[(i+1)%len(points)])
    #print(d)
    return d

def simplify_curve(points, threshold = None, reduction_ratio = None):
    simplified_points = [(p[0], p[1]) for p in points ]

    if threshold is None and reduction_ratio is None:
        print("You must provide either a threshold or a reduction ratio")
        return []

    queue = MinPriorityQueue()
    i = 0
    for p in simplified_points:
        queue.add_point(p, compute_error(simplified_points, p, i))
        i += 1

    candidate = queue.pop_point()
    #print("threshold", threshold, "first", candidate)
    if reduction_ratio is not None:
        desired_size = len(simplified_points)*reduction_ratio
        while(len(simplified_points) > desired_size):
            # find the index of the point in the original sequence
            index = simplified_points.index(candidate[-1])
            # removes the point from the sequence
            simplified_points.remove(candidate[-1])
            #updates the point before
            queue.add_point(simplified_points[index-1],
                    compute_error(simplified_points, simplified_points[index-1]))
            #updates the point after (which now assumes the same index)
            queue.add_point(simplified_points[index%len(simplified_points)],
                    compute_error(simplified_points, simplified_points[index%len(simplified_points)]))
            candidate = queue.pop_point()

    else: #won't bother with optimization now
        while(candidate[0] < threshold):
            index = simplified_points.index(candidate[-1])
            simplified_points.remove(candidate[-1])
            queue.add_point(simplified_points[index-1],
                    compute_error(simplified_points, simplified_points[index-1]))
            queue.add_point(simplified_points[index%len(simplified_points)],
                    compute_error(simplified_points, simplified_points[index%len(simplified_points)]))
            candidate = queue.pop_point()

    return simplified_points



def curve_simplification(filepath, threshold= 1):
    #extract the filename
    filename = os.path.basename(filepath)

    points = read_points(filepath)
    #drawing.draw_closed_curve(points, "original.eps")
    ratios = [i/10 for i in range(1, 10)]
    for r in ratios:
        ax = plt.axes()
        drawing.plot_and_save(os.path.join(os.path.join(results_folder, filename[:-4] + ".png")),
        False, ax, vertices=points, vertices_color="k",
        edges=[(i, (i+1)%len(points)) for i in range(len(points))],
        edges_color = "black")

        simplified_points = simplify_curve(points, reduction_ratio=r)
        #drawing.draw_closed_curve(points, "simplificado.eps")
        ax = plt.axes()
        drawing.plot_and_save(os.path.join(os.path.join(results_folder, "simplified" + str(r*10) +filename[:-4] + ".png")),
        True, ax, vertices=np.array(simplified_points), vertices_color="k",
        edges=[(i, (i+1)%len(simplified_points)) for i in range(len(simplified_points))],
        edges_color = "green")
        print("simplified", filepath, r)

    write_points(simplified_points, os.path.join(results_folder, "simplified_"+filename))

def rec_all_data(folder):
    filenames = [fname for fname in os.listdir(folder) if fname.endswith('.txt')]
    for fname in filenames:
        curve_simplification(os.path.join(data_folder, fname))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) ... The output is written to a txt file
    ''')
    parser.add_argument('input_file', help='input list of points with X Y per line (format: txt)')
    args = parser.parse_args()
    if args.input_file == "all":
        rec_all_data(data_folder)
    else:
        curve_simplification(args.input_file)
