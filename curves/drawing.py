import turtle

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

def draw_segment_set(segments, filename, hold = False):
    colors = ["red", "green", "blue", "yellow"]
    i = 0
    if not segments:
        return
    t = turtle.Turtle()
    t.penup()
    for segment in segments:
        t.goto(segment.first.x, segment.first.y)
        t.pendown()
        t.goto(segment.second.x, segment.second.y)
        t.penup()
        t.color(colors[i])
        i = (i + 1)%4

    ts = t.getscreen()

    ts.getcanvas().postscript(file=filename)
    if hold:
        ts.mainloop()
    try:
        ts.clear()
    except Exception as e:
        pass

def draw_points_set(points, filename, hold = False):
    colors = ["red", "green", "blue", "yellow"]
    i = 0
    if not points:
        return
    t = turtle.Turtle()
    ts = t.getscreen()
    height = ts.canvheight
    t.penup()
    for point in points:
        t.goto(point.x, height - point.y)
        t.pendown()
        t.dot()
        t.penup()

    ts.getcanvas().postscript(file=filename)
    if hold:
        ts.mainloop()
    try:
        ts.clear()
    except Exception as e:
        pass

def draw_voronoi_diagram():
    pass

def draw_delaunay_diagram():
    pass
