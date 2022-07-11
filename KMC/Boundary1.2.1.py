from cmath import sqrt
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as m3d
from mpl_toolkits.mplot3d import Axes3D
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
def midpt(p1, p2):
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

file = "ionchanneltest.pdb"
A = np.array([4853, 4658, 4710, 4668, 4700, 4710, 4727, 4737, 4758, 4795])
B = np.array([18587, 18320, 18489, 18339, 18350, 18370, 18389, 18405, 18421, 18437])
C = np.array([8737, 8853, 8824, 8752, 8771, 8785, 8804, 8824, 8843, 8872])
D = np.array([13799, 13475, 13600, 13495, 13509, 13528, 13542, 13561, 13581, 13623])

#Here the plane is from 3 points, need to find out the best possible plane fit
PlaneA = getplane(file, A[0], A[1], A[2]) 
PlaneB = getplane(file, B[0], B[1], B[2])
PlaneC = getplane(file, C[0], C[1], C[2])
PlaneD = getplane(file, D[0], D[1], D[2])

PlaneAD = (PlaneA + PlaneD)/2 #a,b,c,d
PlaneBC = (PlaneB + PlaneC)/2 #p,q,r,s
ConstPt, ParDir = lineplane(PlaneAD, PlaneBC)

Ri = np.zeros([3])
BoundaryPts = 1000
CrossSecPts = 1000
flagBounds = 0
flagAD = 0
flagBC = 0
Bounds = np.zeros([BoundaryPts,3]) #List of 1000 points that form the boundary
CrossSecAD = np.zeros([CrossSecPts,3]) #List of 500 points that form a cross section
CrossSecBC = np.zeros([CrossSecPts,3]) #List of 500 points that form a cross section
f = open(file,'r')
lines = f.readlines()
for i, line in enumerate(lines):
    vline = line.rstrip()
    words = vline.split()
    if (words[0] == 'ATOM'):
        Ri[0] = float(vline[30:38])
        Ri[1] = float(vline[38:46])
        Ri[2] = float(vline[46:54])
        temp_d = linedist(Ri[0], ConstPt[0], ParDir[0], Ri[1], ConstPt[1], ParDir[2], Ri[2], ConstPt[2], ParDir[2])
        if (temp_d < 100):
            Bounds[flagBounds][0] = Ri[0]
            Bounds[flagBounds][1] = Ri[1]
            Bounds[flagBounds][2] = Ri[2]
            if (flagBounds < BoundaryPts):
                flagBounds += 1
        temp_d_AD = planedist(PlaneAD[0], PlaneAD[1], PlaneAD[2], PlaneAD[3], Ri[0], Ri[1], Ri[2])
        if (temp_d_AD < 100):
            CrossSecAD[flagAD][0] = Ri[0]
            CrossSecAD[flagAD][1] = Ri[1]
            CrossSecAD[flagAD][2] = Ri[2]
            if (flagAD < CrossSecPts-1):
                flagAD += 1
        temp_d_BC = planedist(PlaneBC[0], PlaneBC[1], PlaneBC[2], PlaneBC[3], Ri[0], Ri[1], Ri[2])
        if (temp_d_BC < 100):
            CrossSecBC[flagBC][0] = Ri[0]
            CrossSecBC[flagBC][1] = Ri[1]
            CrossSecBC[flagBC][2] = Ri[2]
            if (flagBC < CrossSecPts-1):
                flagBC += 1
f.close()
o = open('Boundary.dat','w')
for i in range (BoundaryPts):
    o.write(str(Bounds[i][0]) + "   " + str(Bounds[i][1]) + "   " + str(Bounds[i][2]) + "\n")
o.close()
o = open('CrossSecAD.dat','w')
for i in range (CrossSecPts):
    o.write(str(CrossSecAD[i][0]) + "   " + str(CrossSecAD[i][1]) + "   " + str(CrossSecAD[i][2]) + "\n")
o.close()
o = open('CrossSecBC.dat','w')
for i in range (CrossSecPts):
    o.write(str(CrossSecBC[i][0]) + "   " + str(CrossSecBC[i][1]) + "   " + str(CrossSecBC[i][2]) + "\n")
o.close()

fig = plt.figure()
ax = plt.axes(projection='3d')
# # Data for a three-dimensional line
# zline = np.linspace(0, 15, 1000)
# xline = np.sin(zline)
# yline = np.cos(zline)
# ax.plot3D(xline, yline, zline, 'gray')

# Data for three-dimensional scattered points
temp = {(0, 0, 0)}
Bounds1 = []
for idx, row in enumerate(map(tuple, Bounds)):
    if row not in temp:
        Bounds1.append(row)
# print(Bounds1)
Bounds1 = np.array(Bounds1)
o = open('outputs/bounds1.dat','w')
for i in range (len(Bounds1)):
    o.write(str(Bounds1[i][0]) + "   " + str(Bounds1[i][1]) + "   " + str(Bounds1[i][2]) + "\n")
o.close()

# X = np.linspace(-20,20,100)
# Y = np.linspace(-20,20,100)
# X, Y = np.meshgrid(X,Y)
# Z = 4*xx**2 + yy**2
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X, Y, Z, cmap="plasma", linewidth=0, antialiased=False, alpha=0.5)
# plt.show()






xdataB = Bounds1[:,0]
ydataB = Bounds1[:,1]
zdataB = Bounds1[:,2]
# zdataAD = CrossSecAD[:,2]
# xdataAD = CrossSecAD[:,0]
# ydataAD = CrossSecAD[:,1]
# zdataBC = CrossSecBC[:,2]
# xdataBC = Bounds[:,0]
# ydataBC = CrossSecBC[:,1]
# ax.scatter3D(xdataB, ydataB, zdataB, c=zdataB, cmap='Greens')
# # ax.scatter3D(xdataAD, ydataAD, zdataAD, c=zdataAD, cmap='Reds')
# # ax.scatter3D(xdataBC, ydataBC, zdataBC, c=zdataBC, cmap='Blues')
# plt.show()

# xs, ys = np.meshgrid(xdataB, ydataB)
# Blen = len(zdataB)
# zs = np.zeros([Blen,Blen])
# for j, phase in enumerate(ys):
#     zs[j, :] = zdataB[j]
# fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
# plotdata = {'xs': xs, 'ys': ys, 'zs': zs, 'plotfun': plotfun}
# plotdata['plotfun'](plotdata)
# with open("outputs/Bounds.pickle", 'wb') as file:
#     pickle.dump(plotdata, file)
# plt.figure().add_subplot(111, projection='3d'),
# plotdata = pickle.load(open('outputs/Bounds.pickle', 'rb'))
# plotdata['plotfun'](plotdata)
# plt.show()

# lAD = len(A)
# lBC = len(B)
# MidPts = np.zeros([(lAD+lBC),3])
# for i in range (lAD):
#     R1 = getcoord(file, A[i])
#     R2 = getcoord(file, D[i])
#     MidPts[i] = midpt(R1,R2)
# for i in range (lAD, (lAD+lBC)):
#     R1 = getcoord(file, B[i-lAD])
#     R2 = getcoord(file, C[i-lAD])
#     MidPts[i] = midpt(R1,R2)

# # Generate some data that lies along a line
# X = MidPts[:, 0]
# Y = MidPts[:, 1]
# Z = MidPts[:, 2]
# data = np.concatenate((X[:, np.newaxis], Y[:, np.newaxis], Z[:, np.newaxis]), axis=1)
# # Perturb with some Gaussian noise
# data += np.random.normal(size=data.shape) * 0.4
# # Calculate the mean of the points, i.e. the 'center' of the cloud
# datamean = data.mean(axis=0)
# # Do an SVD on the mean-centered data.
# #Journal of Machine Learning Research 16 (2015) 3367-3402
# uu, dd, vv = np.linalg.svd(data - datamean)
# print(dd)
# # Now vv[0] contains the first principal component, i.e. the direction
# # vector of the 'best fit' line in the least squares sense.
# linepts = vv[0] * np.mgrid[160:165:122][:, np.newaxis]
# # shift by the mean to get the line in the right place
# linepts += datamean
# # Verify that everything looks right.
# ax = m3d.Axes3D(plt.figure())
# ax.scatter3D(*data.T)
# ax.plot3D(*linepts.T)
# plt.show()

