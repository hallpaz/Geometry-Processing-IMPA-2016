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
