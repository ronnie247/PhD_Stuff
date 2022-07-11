from cmath import sqrt
import math
import os
import numpy as np
import matplotlib.pyplot as plt

#get equation of the plane
def getplane(filename, A1, A2, A3):
    R1 = np.zeros([3])
    R2 = np.zeros([3])
    R3 = np.zeros([3])

    f = open(filename,'r')
    lines = f.readlines()
    for i, line in enumerate(lines):
        vline = line.rstrip()
        words = vline.split()
        if (words[0] == 'ATOM'):
            xcoord = float(vline[30:38])
            ycoord = float(vline[38:46])
            zcoord = float(vline[46:54])
            if (int(words[1]) == A1):
                R1[0] = xcoord
                R1[1] = ycoord
                R1[2] = zcoord
            if (int(words[1]) == A2):
                R2[0] = xcoord
                R2[1] = ycoord
                R2[2] = zcoord
            if (int(words[1]) == A3):
                R3[0] = xcoord
                R3[1] = ycoord
                R3[2] = zcoord
    # These two vectors are in the plane
    v1 = R3 - R1
    v2 = R2 - R1
    # the cross product is a vector normal to the plane
    cp = np.cross(v1, v2)
    a, b, c = cp
    # This evaluates a * x3 + b * y3 + c * z3 which equals d
    d = np.dot(cp, R3)
    print('The equation is {0}x + {1}y + {2}z = {3}'.format(a, b, c, d))
    return a, b, c, d

#get distance from a plane
def planedist(a, b, c, d, x0, y0, z0):
    dist = (abs((a*x0)+(b*y0)+(c*z0)+d)) / (sqrt((a*a)+(b*b)+(c*c)))
    return dist

#get midpoint between two atoms
def midpt(p1, p2):
    return ((p1+p2)/2)
