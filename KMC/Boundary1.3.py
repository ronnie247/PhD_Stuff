from cmath import sqrt
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as m3d
import matplotlib.animation as ani
import dill as pickle

#get coords of point
def getcoord(filename, Atom):
    R = np.zeros([3])

    f = open(filename,'r')
    lines = f.readlines()
    for i, line in enumerate(lines):
        vline = line.rstrip()
        words = vline.split()
        if (words[0] == 'ATOM'):
            xcoord = float(vline[30:38])
            ycoord = float(vline[38:46])
            zcoord = float(vline[46:54])
            if (int(words[1]) == Atom):
                R[0] = xcoord
                R[1] = ycoord
                R[2] = zcoord
    f.close()
    return R

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
    # print('The equation is {0}x + {1}y + {2}z = {3}'.format(a, b, c, d))
    P = np.array([a, b, c, d])
    f.close()
    return P

#get distance from a plane
def planedist(a, b, c, d, x0, y0, z0):
    dist = (abs((a*x0)+(b*y0)+(c*z0)+d)) / (sqrt((a*a)+(b*b)+(c*c)))
    return dist

#get distance from a line
#Point A = i,j,k ; L = x0+tl, y0+tm, z0+tn
def linedist(i, x0, l, j, y0, m, k, z0=0, n=1):
    P1 = np.array([i, j, k])
    d = np.array([l, m, n])
    P2 = np.array([x0, y0, z0])
    P21 = P2 - P1
    P3 = np.cross(P21, d)
    dist = (np.linalg.norm(P3))/(np.linalg.norm(d))
    return dist

#get distance betweem two points
def pointdist(A, B):
    dist = sqrt(((A[0] - B[0])**2) +((A[1] - B[1])**2)+((A[2] - B[2])**2))
    return dist

#get midpoint between two atoms
def getmidpt(p1, p2):
    return ((p1+p2)/2)

#get line intersecting between two planes
# LineABCD = np.zeros([3])
# #(x,y,z) = (((dq-bs)+(br-cq)t)/(bp-aq), ((dp-as)-(cp-ar)t)/(bp-aq),t)
# LineABCD[0] = ((PlaneAD[3]*PlaneBC[1] - PlaneAD[1]*PlaneBC[3])+(((PlaneAD[1]*PlaneBC[2]) - (PlaneAD[2]*PlaneBC[1]))*i)) / ((PlaneAD[1]*PlaneBC[0])-(PlaneAD[0]*PlaneBC[1]))
# LineABCD[1] = ((PlaneAD[3]*PlaneBC[0] - PlaneAD[0]*PlaneBC[3])-(((PlaneAD[2]*PlaneBC[0]) - (PlaneAD[0]*PlaneBC[2]))*i)) / ((PlaneAD[1]*PlaneBC[0])-(PlaneAD[0]*PlaneBC[1]))
# LineABCD[2] = i
def lineplane(A, B):
    T = np.zeros([3])
    P0 = np.zeros([3])
    #(x,y,z) = (((dq-bs)+(br-cq)t)/(bp-aq), ((dp-as)-(cp-ar)t)/(bp-aq),t)
    P0[0] = (A[3]*B[1] - A[1]*B[3]) / ((A[1]*B[0])-(A[0]*B[1]))
    P0[1] = (A[3]*B[0] - A[0]*B[3]) / ((A[1]*B[0])-(A[0]*B[1]))
    P0[2] = 1 / ((A[1]*B[0])-(A[0]*B[1]))
    T[0] = ((A[1]*B[2]) - (A[2]*B[1])) / ((A[1]*B[0])-(A[0]*B[1]))
    T[1] = ((A[2]*B[0]) - (A[0]*B[2])) / ((A[1]*B[0])-(A[0]*B[1]))
    T[2] = 0 / ((A[1]*B[0])-(A[0]*B[1]))
    return P0, T

#pickle plot fn
def plotfun(plotdata):
    ax = plt.gca()
    ax.plot_surface(plotdata['xs'], plotdata['ys'], plotdata['zs'])
    # ax.scatter(plotdata['xs'], plotdata['ys'], plotdata['zs'])

file = "ionchanneltest.pdb"
A = np.array([4853, 4658, 4710, 4668, 4700, 4710, 4727, 4737, 4758, 4795])
B = np.array([18587, 18320, 18489, 18339, 18350, 18370, 18389, 18405, 18421, 18437])
C = np.array([8737, 8853, 8824, 8752, 8771, 8785, 8804, 8824, 8843, 8872])
D = np.array([13799, 13475, 13600, 13495, 13509, 13528, 13542, 13561, 13581, 13623])
l = len(A)

MidPts = np.zeros([2*l,3])
for i in range (l):
    RA = getcoord(file, A[i])
    RB = getcoord(file, B[i])
    RC = getcoord(file, C[i])
    RD = getcoord(file, D[i])
    AD = getmidpt(RA,RD)
    BC = getmidpt(RB, RC)
    # print (AD)
    for j in range (3):
        MidPts[i][j] = AD[j]
        MidPts[i+l][j] = BC[j]

# xx = MidPts[:,0]
# yy = MidPts[:,1]
# zz = np.zeros([2*l,2*l])
# for j, phase in enumerate(yy):
#     zz[j, :] = MidPts[j,2]
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(xx,yy,zz)
# plt.show()

# x = MidPts[:,0]
# y = MidPts[:,1]
# z = MidPts[:,2]
# # this will find the slope and x-intercept of a plane
# # parallel to the y-axis that best fits the data
# A_xz = np.vstack((x, np.ones(len(x)))).T
# m_xz, c_xz = np.linalg.lstsq(A_xz, z)[0]
# # again for a plane parallel to the x-axis
# A_yz = np.vstack((y, np.ones(len(y)))).T
# m_yz, c_yz = np.linalg.lstsq(A_yz, z)[0]
# # the intersection of those two planes and
# # the function for the line would be:
# # z = m_yz * y + c_yz
# # z = m_xz * x + c_xz
# # or:
# def lin(z):
#     x = (z - c_xz)/m_xz
#     y = (z - c_yz)/m_yz
#     return x,y
# fig = plt.figure()
# ax = m3d.Axes3D(fig)
# zz = np.linspace(20,80)
# xx,yy = lin(zz)
# ax.scatter(x, y, z)
# ax.plot(xx,yy,zz)
# plt.savefig('outputs/test1.png')
# plt.show()



# xs, ys = np.meshgrid(MidPts[:,0], MidPts[:,1])
# zs = np.zeros([2*l,2*l])
# for j, phase in enumerate(ys):
#     zs[j, :] = MidPts[j,2]
# fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
# plotdata = {'xs': xs, 'ys': ys, 'zs': zs, 'plotfun': plotfun}
# plotdata['plotfun'](plotdata)
# with open("MidPts.pickle", 'wb') as file:
#     pickle.dump(plotdata, file)
# plt.figure().add_subplot(111, projection='3d'),
# plotdata = pickle.load(open('MidPts.pickle', 'rb'))
# plotdata['plotfun'](plotdata)
# plt.show()