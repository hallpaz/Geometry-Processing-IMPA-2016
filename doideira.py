import turtle
import numpy as np

points = [np.array([2, 3]), np.array([8, 12]), np.array([ 14, 3])]

def draw(points, filename):
    if not points:
        return
    t = turtle.Turtle()

    t.penup()
    t.goto(points[0][0], points[0][1])
    t.pendown()
    for point in points[1:]:
        t.goto(point[0], point[1])
    t.goto(points[0][0], points[0][1])

    ts = t.getscreen()

    ts.getcanvas().postscript(file=filename)

    t.getscreen().mainloop()


if __name__ == '__main__':
    draw(points, "hello.eps")




# t = turtle.Turtle()
#
# t.penup()
# t.goto(points[0][0], points[0][1])
# t.pendown()
# for point in points[1:]:
#     t.goto(point[0], point[1])
# t.goto(points[0][0], points[0][1])
#
# ts = t.getscreen()
#
# #ts.getcanvas().postscript(file=filename)
#
# t.getscreen().mainloop()
