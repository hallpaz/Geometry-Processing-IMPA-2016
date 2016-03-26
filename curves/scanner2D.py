from PIL import Image
from drawing import draw_segments, draw_points_set, draw_segment_set
from Point2D import Point2D, CCW_test
from Segment import Segment

import argparse
import math
import os

BIG = 9999999

def getsample(image, position, dst):

    direction = dst - position
    if direction.x > 0: #->
        start_col = int(position.x)
        end_col = image.width
    elif direction.x < 0: #<-
        start_col = 0
        end_col = int(position.x)
    else:
        start_col = int(position.x)
        end_col = start_col + 1

    if direction.y > 0: #->
        start_row = int(position.y)
        end_row = image.height
    elif direction.y < 0: #<-
        start_row = 0
        end_row = int(position.y)
    else:
        start_row= int(position.y)
        end_row = start_row + 1

    ray = Segment(position, dst)
    intersection = None
    #print(image.getpixel((398,398)))
    intersections = []
    for j in range(start_row, end_row-2):
        for i in range(start_col, end_col-2):
            #print("pixel: ", i, j)#image.getpixel((i,j)))
            t = ray.ray_intersection(Segment(Point2D(i, j), Point2D(i+1, j)))
            if t is not None and (image.getpixel((i,j)) == 0 and image.getpixel((i+1, j)) == 0):
                intersections.append(t)
            t = ray.ray_intersection(Segment(Point2D(i+1, j), Point2D(i+1, j+1)))
            if t is not None and (image.getpixel((i+1,j)) == 0 and image.getpixel((i+1, j+1)) == 0):
                intersections.append(t)
            t = ray.ray_intersection(Segment(Point2D(i+1, j+1), Point2D(i, j+1)))
            if t is not None and (image.getpixel((i+1,j+1)) == 0 and image.getpixel((i, j+1)) == 0):
                intersections.append(t)
            t = ray.ray_intersection(Segment(Point2D(i, j+1), Point2D(i, j)))
            if t is not None and (image.getpixel((i,j+1)) == 0 and image.getpixel((i, j)) == 0):
                intersections.append(t)
    if intersections:
        mind = 9999999
        for p in intersections:
            d = Point2D.distance(position, p)
            if d < mind:
                mind = d
                intersection = p
        return intersection

        # t = ray.intersects(Segment(Point2D(i, j), Point2D(i+1, j)))
        # if t :#and (image.getpixel((i,j)) == 0 and image.getpixel((i+1, j)) == 0):
        #     intersections.append(Point2D(i,j))
        # t = ray.intersects(Segment(Point2D(i+1, j), Point2D(i+1, j+1)))
        # if t :#and (image.getpixel((i+1,j)) == 0 and image.getpixel((i+1, j+1)) == 0):
        #     intersections.append(Point2D(i+1,j))
        # t = ray.intersects(Segment(Point2D(i+1, j+1), Point2D(i, j+1)))
        # if t :#and (image.getpixel((i+1,j+1)) == 0 and image.getpixel((i, j+1)) == 0):
        #     intersections.append(Point2D(i+1,j+1))
        # t = ray.intersects(Segment(Point2D(i, j+1), Point2D(i, j)))
        # if t :#and (image.getpixel((i,j+1)) == 0 and image.getpixel((i, j)) == 0):
        #     intersections.append(Point2D(i,j+1))
        # if intersections:
        #     mind = 9999999
        #     for p in intersections:
        #         d = Point2D.distance(position, p)
        #         if d < mind:
        #             mind = d
        #             intersection = p
        #     return intersection
    return None

def get_ortho_sample(image, position, dst):

    direction = dst - position
    incx = 1
    incy = 1
    if direction.x > 0: #->
        start_col = int(position.x)+1
        end_col = image.width-1
    elif direction.x < 0: #<-
        start_col = int(position.x)-1
        end_col = 0
        incx = -1
    else:
        start_col = int(position.x)
        end_col = start_col + 1
        #print("cols: ", start_col, end_col)

    if direction.y > 0: #->
        start_row = int(position.y)+1
        end_row = image.height-1
        #print("start row", start_row)
    elif direction.y < 0: #<-
        start_row = int(position.y)-1
        end_row = 0
        incy = -1
    else:
        start_row= int(position.y)
        end_row = start_row + 1

    ray = Segment(position, dst)
    intersection = None
    #print(image.getpixel((398,398)))
    intersections = []
    for j in range(start_row, end_row, incy):
        for i in range(start_col, end_col, incx):
            #print("pixel: ", i, j, image.getpixel((i,j)))
            t = ray.ray_intersection(Segment(Point2D(i, j), Point2D(i+1, j)))
            if t is not None and (image.getpixel((i,j)) == 0 and image.getpixel((i+1, j)) == 0):
                #intersections.append(t)
                return t
            t = ray.ray_intersection(Segment(Point2D(i+1, j), Point2D(i+1, j+1)))
            if t is not None and (image.getpixel((i+1,j)) == 0 and image.getpixel((i+1, j+1)) == 0):
                #intersections.append(t)
                return t
            t = ray.ray_intersection(Segment(Point2D(i+1, j+1), Point2D(i, j+1)))
            if t is not None and (image.getpixel((i+1,j+1)) == 0 and image.getpixel((i, j+1)) == 0):
                #intersections.append(t)
                return t
            t = ray.ray_intersection(Segment(Point2D(i, j+1), Point2D(i, j)))
            if t is not None and (image.getpixel((i,j+1)) == 0 and image.getpixel((i, j)) == 0):
                #intersections.append(t)
                return t
    # if intersections:
    #     mind = 9999999
    #     for p in intersections:
    #         d = Point2D.distance(position, p)
    #         if d < mind:
    #             mind = d
    #             intersection = p
    #     return intersection
    return None

def get_samples_from_image(image, num_samples = 100):
    minx = BIG
    maxx = -BIG
    miny = BIG
    maxy = -BIG
    samples = []
    for j in range(image.height):
        for i in range(image.width):
            if image.getpixel((i,j)) == 0: #is black
                if i < minx:
                    minx = i
                if i > maxx:
                    maxx = i
                if j < miny:
                    miny = j
                if j > maxy:
                    maxy = j

    points = [Point2D(minx, miny), Point2D(minx, maxy), Point2D(maxx, maxy), Point2D(maxx, miny)]
    #draw_segments(points, "original.eps")

    print("box: ", (minx, miny), (maxx, maxy))
    qsamples = num_samples//4
    center = Point2D(minx + (maxx - minx)/2, miny + (maxy - miny)/2)
    incy = (maxy - miny-1)//(qsamples-1)
    for j in range(qsamples):
        pos = Point2D(0, miny + j*incy)
        sample = get_ortho_sample(image, pos, Point2D(center.x, miny + j*incy))
        if sample is not None:
            samples.append(sample)
            print("sample: ", i)

    incx = (maxx - minx-1)//(qsamples-1)
    for i in range(qsamples):
        pos = Point2D(minx + i*incx, 0)
        sample = get_ortho_sample(image, pos, Point2D(minx + i*incx, center.y))
        if sample is not None:
            samples.append(sample)
            print("sample: ", i)

    xlimit = image.width-1
    for j in range(qsamples):
        pos = Point2D(xlimit, miny + j*incy)
        sample = get_ortho_sample(image, pos, Point2D(center.x, miny + j*incy))
        if sample is not None:
            samples.append(sample)
            print("sample: ", i)

    ylimit = image.height-1
    for i in range(qsamples):
        pos = Point2D(minx + i*incx, ylimit)
        sample = get_ortho_sample(image, pos, Point2D(minx + i*incx, center.y))
        if sample is not None:
            samples.append(sample)
            print("sample: ", i)

    #draw_segment_set(seg_list, "segmentos.eps")
    return samples


def get_circles_from_image(image, num_samples = 100):
    minx = BIG
    maxx = -BIG
    miny = BIG
    maxy = -BIG
    samples = []
    for j in range(image.height):
        for i in range(image.width):
            if image.getpixel((i,j)) == 0: #is black
                if i < minx:
                    minx = i
                if i > maxx:
                    maxx = i
                if j < miny:
                    miny = j
                if j > maxy:
                    maxy = j

    points = [Point2D(minx, miny), Point2D(minx, maxy), Point2D(maxx, maxy), Point2D(maxx, miny)]
    draw_segments(points, "original.eps")

    maxx += 1
    minx -= 1
    maxy += 1
    miny -= 1
    print((minx, miny), (maxx, maxy))
    center = Point2D(minx + (maxx - minx)/2, miny + (maxy - miny)/2)
    radius = math.sqrt(((maxx - minx)/2)**2 + ((maxy - miny)/2)**2)
    print(center, radius)
    inc = 2*math.pi/num_samples

    seg_list = []

    for i in range(num_samples):
        pos = center + radius*Point2D(math.cos(i*inc), math.sin(i*inc))
        sample = getsample(image, pos, center)
        seg_list.append(Segment(pos, center))
        if sample is not None:
            samples.append(sample)
            print("sample: ", i)

    draw_segment_set(seg_list, "circle.eps")
    return samples


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads an image and simulates a scanner 2D The output is written to a txt file
    ''')
    parser.add_argument('image_file', help='input list of points with X Y per line (format: txt)')
    #parser.add_argument('output_file', help='output list of points (format: txt)')
    parser.add_argument('num_samples', help='Number of sensors')

    args = parser.parse_args()

    basename = args.image_file
    bar_index = basename.rfind("/")
    if bar_index <= 0:
        image_file = os.path.join("images", basename)

    image = Image.open(image_file)
    image = image.convert('1')
    num_samples = int(args.num_samples)

    samples = get_samples_from_image(image, num_samples)
    #samples = get_circles_from_image(image, num_samples)

    output_file = basename
    dot_index = output_file.rfind(".")
    if  dot_index > 0:
        output_file = output_file[:dot_index]

    data_file = os.path.join("data", output_file + ".txt")
    Point2D.write_to_file(samples, data_file)
    output_file = os.path.join("sensor", output_file + str(num_samples) + "_sensors.eps")
    draw_points_set(samples, output_file, True)
