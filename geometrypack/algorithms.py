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
