from Point2D import Point2D
from Segment import Segment
import turtle
import random
import math

def peucker_reduction(points, error_tolerance):
    dmax = 0
    index = 0
    fitting_line = Segment(points[0], points[-1])
    i = 1
    for point in points[1:len(points)-1]:
        d = fitting_line.distance_to_point(point)
        if d > dmax:
            print(d)
            index = i
            dmax = d
        i += 1

    if dmax > error_tolerance:
        print(error_tolerance)
        res1 = peucker_reduction(points[0: index], error_tolerance)
        res2 = peucker_reduction(points[index:], error_tolerance)
        result = res1 + res2
    else:
        result = [points[0], points[-1]]

    return result


def draw_segments(points, filename, hold = False):
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
    try:
        ts.clear()
    except Exception as e:
        pass



def circle(radius = 200, center = 0, samples = 100):
    points = []
    inc = 2*math.pi/samples
    for i in range(samples-1):
        points.append(Point2D(radius*math.cos(i*inc), radius*math.sin(i*inc)))

    return points


if __name__ == '__main__':
    #points = Point2D.read_from_file("input_points.txt")

    points = circle()

    print("Number of points ", len(points))

    reduced_points = peucker_reduction(points, 100)

    Point2D.write_to_file(reduced_points, "output_points.txt")

    draw_segments(points, "original_points.eps")
    draw_segments(reduced_points, "reduced_points.eps", True)
    print("Number of reduced points ", len(reduced_points))
