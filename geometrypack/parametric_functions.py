import numpy as np
from math import cos, sin

def sphere(u, v, radius= 1):
    return np.array([radius*cos(u)*sin(v), radius*sin(u)*cos(v), radius*cos(v)])

def torus(u, v, radius):
    a = radius/4
    return np.array([(radius + a*cos(v))*cos(u), (radius + a*cos(v))*sin(u), a*sin(v)])

def cylinder(u, h, radius= 1):
    return np.array([radius*cos(u), radius*sin(u), h])

def heart(u, v):
    pass
