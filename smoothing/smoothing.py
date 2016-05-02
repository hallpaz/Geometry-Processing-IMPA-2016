import argparse
import drawing
import matplotlib.pyplot as plt
import numpy as np
import os
import triangle
from collections import deque

data_folder = "data"
images_folder = "images"

SIMPLE = 1
LAPLACIAN = 2

########################### Fuctions for Curves #############################
def simple_mean_curve_smoothing(points):
    if len(points) > 3:
        # python deals with negative indices as expected
        return [points[i-1]/2 + points[i+1%len(points)]/2
                    for i in range(len(points))]

    return points

def simple_mean_smoothing_curve_with_anchor_points(points, anchor_indices):
    if len(points) > 3:
        # python deals with negative indices as expected
        return [points[i-1]/2 + points[i+1%len(points)]/2
                for i in range(len(points))
                if i not in anchor_indices]

    return points


def laplacian_curve_smoothing(points):
    if len(points) > 3:
        # python deals with negative indices as expected
        return [points[i-1]/4 + points[i]/2 + points[i+1%len(points)]/4
                for i in range(len(points))]

    return points


def laplacian_smoothing_curve_with_anchor_points(points, anchor_indices):
    if len(points) > 3:
        # python deals with negative indices as expected
        return [points[i-1]/4 + points[i]/2 + points[i+1%len(points)]/4
                for i in range(len(points))
                if i not in anchor_indices]

    return points
##########################################################

########################### Fuctions for Meshes #############################

# def simple_mean_mesh_smoothing(points, graph):
#     if len(points) > 3:
#         # initialization
#         for node in graph:
#             node.color = Color.White
#         # initialize a list with the same length as points
#         smoothed_points = [None for i in range(len(points))]
#         queue = deque()
#         for p in graph:
#             if p.color is Color.White:
#                 queue.append(p)
#                 while(queue):
#                     node = queue.popleft()
#                     i = 0
#                     newvalue = np.array([0 for i in range(len(node.value))])
#                     for n in node.neighbors:
#                         i += 1
#                         newvalue += n.value
#                         if n.color is Color.White:
#                             # mark as 'in progress'
#                             n.color = Color.Gray
#                             queue.append(n)
#                     # take the mean
#                     newvalue /= i
#                     smoothed_points[node.id] = newvalue
#                     # mark as 'finished'
#                     node.color = Color.Black
#
#         return smoothed_points = [points[i-1]/2 + points[i+1%len(points)]/2
#                                         for i in range(len(points))]
#
#     return points

def simple_mean_smoothing_mesh(points, anchor_indices = []):
    if len(points) > 3:
        # initialization
        for node in graph:
            if node.id in anchor_indices:
                node.color = Color.Black
            else:
                node.color = Color.White
        # initialize a list with the same length as points
        smoothed_points = [None for i in range(len(points))]
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
                        newvalue += n.value
                        if n.color is Color.White:
                            # mark as 'in progress'
                            n.color = Color.Gray
                            queue.append(n)
                    # take the mean
                    newvalue /= i
                    smoothed_points[node.id] = newvalue
                    # mark as 'finished'
                    node.color = Color.Black

        return smoothed_points

    return points


# def laplacian_mesh_smoothing(points):
#     if len(points) > 3:
#         # python deals with negative indices as expected
#         return smoothed_points = [points[i-1]/4 + points[i]/2 + points[i+1%len(points)]/4
#                                     for i in range(len(points))]
#
#     return points
#
#
# def laplacian_smoothing_mesh_with_anchor_points(points, anchor_indices):
#     if len(points) > 3:
#         # python deals with negative indices as expected
#         return smoothed_points = [points[i-1]/4 + points[i]/2 + points[i+1%len(points)]/4
#                                     for i in range(len(points))
#                                     if i not in anchor_indices]
#
#     return points

def read_points(filename: str)-> list:
    points = []
    with open(filename) as myfile:
        file_lines = myfile.readlines()
        for line in file_lines:
            content = line.split()
            content = [float(n) for n in content]
            # each element is a numpy array
            points.append(np.array(content))
    return points

def write_to_file(points:list, filename:str):
    if not points:
        return None
    with open(filename, "w") as myfile:
        for point in points:
            if len(point) == 2:
                myfile.write("{} {}\n".format(point[0], point[1]))
            elif len(point) == 3:
                myfile.write("{} {} {}\n".format(point[0], point[1], point[3]))
            else:
                raise Exception("Points should have dimension 2 or 3")



def smooth_mesh(filename, method, anchor_boundary = False):
    vertices, indices = read_OFF(filename)
    if anchor_boundary:
        boundary = []
    if method == SIMPLE:
        simple_mean_smoothing_mesh(vertices, graph, boundary)
    elif method == LAPLACIAN:
        laplacian_mean_smoothing_mesh(vertices, graph, boundary)
    else:
        raise Exception("Invalid Method")


def smooth_curve(filename, method= 1):
    vertices = read_points(filename)
    output_file = filename[:-4]
    # segs = np.array([[i, (i+1)%len(vertices)] for i in range(len(vertices))])
    # print(segs)
    # #drawing.plot_and_save(output_file + "original.png", True, ax=plt.axes(), vertices=vertices,
    #     segments= segs)
    print("before draw")
    drawing.draw_segments(vertices, output_file + "original.png")

    if anchor_boundary:
        boundary = []
    if method == SIMPLE:
        vertices = simple_mean_curve_smoothing(vertices)
        segs = np.array([[i, (i+1)%len(vertices)] for i in range(len(vertices))])

        drawing.plot_and_save(output_file + "_simple1.png", True, ax=plt.axes(), vertices= vertices,
            segments= segs)
    elif method == LAPLACIAN:
        laplacian_mean_curve_smoothing(vertices)
    else:
        raise Exception("Invalid Method")

def main():
    smooth_curve(os.path.join(data_folder, "fronteira0.txt"))

if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description='''
    # This script reads a list of 2D points (x, y) ... The output is written to a txt file
    # ''')
    # parser.add_argument('input_file', help='input list of points with X Y per line (format: txt)')
    # args = parser.parse_args()
    # if args.input_file == "all":
    #     #rec_all_data(data_folder)
    #     pass
    # elif args.input_file == "soup.off":
    #     draw_boundary(args.input_file)
    #
    # else:
    #     main(os.path.join(data_folder, args.input_file))
    main()
