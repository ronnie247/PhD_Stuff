from math import sqrt
import os
import numpy as np
import matplotlib.pyplot as plt

A1 = 13799
A2 = 13475
A3 = 13600

R1 = np.zeros([3])
R2 = np.zeros([3])
R3 = np.zeros([3])

f = open("ionchanneltest.pdb",'r')
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
# print(xcoord)
# print(ycoord)
# print(zcoord)
#Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# numpy.linspace(start, stop, num=number of points needed, endpoint=True, retstep=False, dtype=None, axis=0)
x = np.linspace(0, 160, 10)
y = np.linspace(0, 165, 10)
X, Y = np.meshgrid(x, y)
Z = (d - a * X - b * Y) / c

# plot the mesh. Each array is 2D, so we flatten them to 1D arrays
ax.plot(X.flatten(), Y.flatten(), Z.flatten(), 'bo ')
# plot the original points. We use zip to get 1D lists of x, y and z coordinates.
ax.plot(*zip(R1, R2, R3), color='r', linestyle=' ', marker='o')
# adjust the view so we can see the point/plane alignment
# ax.view_init(0, 22)
plt.tight_layout()
plt.savefig('plane.png')
plt.show()