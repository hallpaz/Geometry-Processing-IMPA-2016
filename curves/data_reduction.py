from Point2D import Point2D
from Segment import Segment
import argparse
import turtle
import math

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


def draw_segments(points, filename, hold = False):
    if not points:
        return
    t = turtle.Turtle()

    t.penup()
    t.goto(points[0].x, points[0].y)
    t.pendown()
    for point in points[1:]:
        t.goto(point.x, point.y)
    t.goto(points[0].x, points[0].y)

    ts = t.getscreen()

    ts.getcanvas().postscript(file=filename)
    if hold:
        ts.mainloop()
    # try:
    #     ts.clear()
    # except Exception as e:
    #     pass



def circle(radius = 200, center = 0, samples = 100):
    points = []
    inc = 2*math.pi/samples
    for i in range(samples):
        points.append(Point2D(radius*math.cos(i*inc), radius*math.sin(i*inc)))

    return points

def reduce_circle():
    data_folder = "data/"
    images_folder = "images/"

    points = circle()

    print("Number of points ", len(points))

    reduced_points = peucker_reduction(points, 10)

    Point2D.write_to_file(reduced_points, data_folder + "output_points.txt")

    draw_segments(points, images_folder + "original_circle_points.eps")
    draw_segments(reduced_points, images_folder + "reduced_circle_points.eps")
    print("Number of reduced points ", len(reduced_points))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) and Douglas Peucker data reduction algorithm. The output is written to a txt file
    ''')
    parser.add_argument('input_file', help='input list of points with X Y per line (format: txt)')
    parser.add_argument('output_file', help='output list of points (format: txt)')
    parser.add_argument('tolerance', help='True if images from the two set of points should be drawn')
    #parser.add_argument('should_draw', help='True if images from the two set of points should be drawn')
    args = parser.parse_args()

    points = Point2D.read_from_file(args.input_file)

    print("Number of points ", len(points))
    tolerance = float(args.tolerance)
    reduced_points = peucker_reduction(points, tolerance)

    Point2D.write_to_file(reduced_points, args.output_file)

    draw_segments(points, "original_points.eps")
    draw_segments(reduced_points, "reduced_points.eps")
    print("Number of reduced points ", len(reduced_points))
