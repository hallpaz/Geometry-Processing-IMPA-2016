import numpy as np


def distance(points, edge):
    p1 = points[edge[0]]
    p2 = points[edge[1]]
    #print(p1, p2)
    return np.linalg.norm(p1-p2)

def cos_angle(pivot, q, s):
    pq = q - pivot
    ps = s - pivot
    return np.dot(pq, ps)/(np.linalg.norm(pq)*np.linalg.norm(ps))

def edges_from_triangle(triangle):
    return [ (triangle[0], triangle[1]), (triangle[1], triangle[2]), (triangle[2], triangle[0]) ]


def normalize(v):
    norm=np.linalg.norm(v)
    if norm < 0.00001:
       return v
    return v/norm

def distance_point_line(pnt, start, end):
    line_vec = np.array(start) - np.array(end)
    pnt_vec = np.array(start) - np.array(pnt)
    line_len = np.linalg.norm(line_vec)
    line_unitvec = normalize(line_vec)
    pnt_vec_scaled = pnt_vec * (1.0/line_len)
    t = np.dot(line_unitvec, pnt_vec_scaled)

    nearest = line_vec*t
    dist = np.linalg.norm(nearest - pnt_vec)

    return dist


def planeFit(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    from: http://stackoverflow.com/questions/12299540/plane-fitting-to-4-or-more-xyz-points
    """
    import numpy as np
    from numpy.linalg import svd
    
    points = [[p[0] for p in points], [p[1] for p in points], [p[2] for p in points]]
    points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trialing dimensions
    #print(points)
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:,np.newaxis]
    M = np.dot(x, x.T) # Could also use np.cov(x) here.
    return ctr, svd(M)[0][:,-1]

def distance_vertex_plane(vertex, planepoint, planenormal):
    v = np.array(vertex) - np.array(planepoint)
    d = abs(np.dot(v, np.array(normalize(planenormal))))
    return d
