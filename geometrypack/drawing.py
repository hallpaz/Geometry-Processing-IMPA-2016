import turtle
import matplotlib.pyplot as plt

# points is a lis of numpy arrays
def draw_segments(points, filename, hold = False):
    if not points:
        return

    t = turtle.Turtle()
    print("TURTLE")

    t.penup()
    t.goto(points[0][0], points[0][1])
    t.pendown()
    for point in points[1:]:
        t.goto(point[0], point[1])
    t.goto(points[0][0], points[0][1])

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


def plot(ax, **kw):

    #vertices(ax, **kw)
    ax.axes.set_aspect('equal')

    if 'segments' in kw: segments(ax, **kw)
    if 'triangles' in kw: triangles(ax, **kw)
    if 'holes' in kw: holes(ax, **kw)
    if 'edges' in kw: edges(ax, **kw)

    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

def plot_and_save(fname, should_close, ax, **kw):
    plot(ax, **kw)
    #plt.savefig(fname, format='png', transparent=True)
    plt.savefig(fname, format='png', transparent=True, bbox_inches='tight', pad_inches=0)
    if should_close:
        plt.close()


def vertices(ax, **kw):
    verts = kw['vertices']
    ax.scatter(*verts.T, color= kw['vertices_color'] if 'vertices_color' in kw else 'k')
    if 'labels' in kw:
        for i in range(verts.shape[0]):
            ax.text(verts[i,0], verts[i,1], str(i))
    if 'markers' in kw:
        vm = kw['vertex_markers']
        for i in range(verts.shape[0]):
            ax.text(verts[i,0], verts[i,1], str(vm[i]))

def segments(ax, **kw):
    verts = kw['vertices']
    segs = kw['segments']
    for beg, end in segs:
        x0, y0 = verts[beg,:]
        x1, y1 = verts[end,:]
        ax.fill([x0, x1], [y0, y1],
                facecolor='none', edgecolor= kw['segments_color'] if 'segments_color' in kw else 'r', linewidth=1,
                zorder=0)

def triangles(ax, **kw):
    verts = kw['vertices']
    ax.triplot(verts[:,0], verts[:,1], kw['triangles'], 'k-')

def holes(ax, **kw):
    ax.scatter(*kw['holes'].T, marker='x', color= kw['holes_color'] if 'holes_color' in kw else 'r')

def edges(ax, **kw):
    """
    Plot regular edges and rays (edges whose one endpoint is at infinity)
    """
    verts = kw['vertices']
    edges = kw['edges']
    for beg, end in edges:
        x0, y0 = verts[beg, :]
        x1, y1 = verts[end, :]
        ax.fill([x0, x1], [y0, y1], facecolor= kw['face_color'] if 'face_color' in kw else 'none',
        edgecolor= kw['edges_color'] if 'edges_color' in kw else 'k', linewidth=.5)

    if ('ray_origins' not in kw) or ('ray_directions' not in kw):
        return

    lim = ax.axis()
    ray_origin = kw['ray_origins']
    ray_direct = kw['ray_directions']
    for (beg, (vx, vy)) in zip(ray_origin.flatten(), ray_direct):
        x0, y0 = verts[beg, :]
        scale = 100.0 # some large number
        x1, y1 = x0 + scale*vx, y0 + scale*vy
        ax.fill([x0, x1], [y0, y1], facecolor='none', edgecolor='k', linewidth=.5)
    ax.axis(lim) # make sure figure is not rescaled by ifinite ray
