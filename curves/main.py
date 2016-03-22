import argparse
import numpy as np
from boundary_tracing import pavlids_boundary_tracing
from data_reduction import peucker_reduction, draw_segments
from os import path
from PIL import Image
from Point2D import Point2D


data_folder = "data"
images_folder = "images"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) and Douglas Peucker data reduction algorithm. The output is written to a txt file
    ''')
    parser.add_argument('image_file', help='input list of points with X Y per line (format: txt)')
    #parser.add_argument('output_file', help='output list of points (format: txt)')
    parser.add_argument('tolerance', help='True if images from the two set of points should be drawn')

    args = parser.parse_args()

    # image = Image.open(path.join(images_folder, args.image_file))
    # image = image.convert('1')
    # image.save("binary.jpg")
    #
    # image_data = np.invert(np.asarray(image).transpose())
    # image_data = image_data.tolist()


    # image_data =    [[ False , False , False , False , False , False, False ] ,
    #                  [ False , False , True , True , True , True, False ] ,
    #                  [ False , False , True , True , True , True, False ] ,
    #                  [ False , True , True , True , True , True, False ] ,
    #                  [ True , True , True , True , True , False, False ] ,
    #                  [ True , True , True , True , True , False, False ] ,
    #                  [ False , False , True , True , True , False, False ] ,
    #                  [ False , False , True , True , True , False, False ] ,
    #                  [ False , False , True , True , True , False, False ] ,
    #                  [ False , False , False , False , False , False, False ] ]


    image_data = []
    false_line = [False]*400
    for i in range(400):
        image_data.append(false_line)

    flag = False
    for i in range(400):
        for j in range(400):
            if image_data[i][j]:
                print("before", i, j)
                flag = True
                break
        if flag:
            break
    print("-------------------------")
    print(image_data[0])

    for i in range(200, 300):
        for j in range(100, 300):
            image_data[i][j] = True

    flag = False
    for i in range(400):
        for j in range(400):
            if image_data[i][j]:
                print("after", i, j)
                flag = True
                break
        if flag:
            break
    print("-------------------------")

    print(image_data[0])


    boundary = pavlids_boundary_tracing(image_data)
    #image_data = None

    output_file = args.image_file
    dot_index = output_file.rfind(".")
    if  dot_index > 0:
        output_file = output_file[:dot_index]

    draw_segments(boundary, path.join(images_folder, "original_" + output_file + ".eps"))

    tolerance = float(args.tolerance)

    reduced_points = peucker_reduction(boundary, tolerance)


    Point2D.write_to_file(reduced_points, path.join(data_folder, output_file + ".txt"))

    draw_segments(reduced_points, path.join(images_folder, "reduced_" + output_file + ".eps"))
