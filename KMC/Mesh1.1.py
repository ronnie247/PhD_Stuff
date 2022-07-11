from math import ceil, floor, sqrt
from random import randint, random
from telnetlib import XASCII
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

# def find_colour(_val):
#     # Colour value constants
#     _colours = {"blue": [0.0, 0.0, 1.0], "green": [0.0, 1.0, 0.00], "yellow": [1.0, 1.0, 0.0], "red": [1.0, 0.0, 0.0]}
#     # Map the value to a colour
#     _colour = [0, 0, 0]
#     if _val > 30:
#         _colour = _colours["red"]
#     elif _val > 20:
#         _colour = _colours["blue"]
#     elif _val > 10:
#         _colour = _colours["green"]
#     elif _val > 0:
#         _colour = _colours["yellow"]
#     return tuple(_colour)

array = np.loadtxt('bounds1.dat')
Blen = len(array)
# # print(Blen)
# data = np.zeros([])
# data = np.random.rand(80, 80) * 20
# # create discrete colormap
# cmap = colors.ListedColormap(['white', 'black'])
# bounds = [0,80,80]
# norm = colors.BoundaryNorm(bounds, cmap.N)
# fig, ax = plt.subplots()
# ax.imshow(data, cmap=cmap, norm=norm)
# # draw gridlines
# ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
# ax.set_xticks(np.arange(20, 100, 1));
# ax.set_yticks(np.arange(20, 100, 1));
# plt.show()

# fig = plt.figure()
# ax = Axes3D(fig)
# df = pd.DataFrame(array, columns = ['X','Y','Z'])
# max0 = df['X'].max()
# max1 = df['X'].max()
# max2 = df['X'].max()
# min0 = df['X'].min()
# min1 = df['X'].min()
# min2 = df['X'].min()
# # ax.plot_trisurf(df.X, df.Y, df.Z, cmap=cm.jet, linewidth=0.2)
# # plt.show()
# # dfback = pd.DataFrame(arrayback, columns = ['X','Y','Z'])
# x = df["X"].tolist()
# y = df["Y"].tolist()
# z = df["Z"].tolist()
# c = np.zeros([len(x)])
# for i in range (len(x)):
#   c[i] = 1+randint(1,10)
#   # c[i] = 1 + conditions(i)
# df['Value'] = c
# colo = df["Value"].tolist()
# # creating 3d figures
# # fig = plt.figure(figsize=(8, 5))
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# # configuring colorbar
# color_map = cm.ScalarMappable(cmap=cm.binary)
# color_map.set_array(colo)
# # creating the heatmap
# img = ax.scatter(x, y, z, marker='s', s=99, color='grey')
# plt.colorbar(color_map)
# # adding title and labels
# ax.set_title("3D Heatmap")
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('')
# # displaying plot
# plt.show()


# array_new = np.zeros([Blen,3])
# for i in range (Blen):
#   for j in range (3):
#     array_new[i][j] = round(array[i][j]) + 0.5
#   for k in range (2):
#     xy_array[i][j] = round() 


# ax = sns.heatmap(xz_array, linewidth=0.5)
# plt.show()

# L = len(array[:,0])
# xdata = np.zeros([L])
# ydata = np.zeros([L])
# zdata = np.zeros([L])
# x1 = np.zeros([L])
# y1 = np.zeros([L])
# z1 = np.zeros([L])
# for i in range (L):
#     xdata[i] = (round(array[i][0])/1.02)-63
#     ydata[i] = (round(array[i][1])/1.02)-9
#     zdata[i] = (round(array[i][2])/1.02)-32
#     x1[i] = xdata[i]
#     y1[i] = ydata[i]*sqrt(1/2) - zdata[i]*sqrt(1/2) + 11
#     z1[i] = ydata[i]*sqrt(1/2) + zdata[i]*sqrt(1/2) - 8
# fig, ax = plt.subplots()
# # plt.scatter(zdata,ydata, marker='s',c='white', linewidths=2, s=128, label= "Allowed region")
# plt.scatter(y1,z1, marker='s',c='white', linewidths=2, s=128, label= "Allowed region")
# ax.xaxis.set_minor_locator(MultipleLocator(0.5))
# ax.yaxis.set_minor_locator(MultipleLocator(0.5))
# ax.set_xticks(np.arange(-.5, 35, 0.5), minor=True)
# ax.set_yticks(np.arange(-.5, 35, 0.5), minor=True)
# # ax.grid(which='minor', alpha=0.2, linestyle = '-')
# ax.grid(which='major', alpha=0.5, linestyle = '-')
# # plt.xlabel("x-coordinate (scaled) in Angstroms")
# plt.ylabel("y-coordinate (scaled) in Angstroms")
# plt.xlabel("z-coordinate (scaled) in Angstroms")
# # plt.ylabel("z-coordinate (scaled) in Angstroms")
# ax.set_facecolor('tab:gray')
# plt.title("Cross Section of the ion channel showing the permeation path")
# fig.tight_layout()
# plt.grid()
# plt.show()
# xmin = round(xdata.min(),1)
# ymin = round(ydata.min(),1)
# zmin = round(zdata.min(),1)
# xmax = round(xdata.max(),1)
# ymax = round(ydata.max(),1)
# zmax = round(zdata.max(),1)


# print(xmin)
# print(xmax)
# print(ymin)
# print(ymax)
# xcol = 'white' if condition else b
# ycol = 
# zcol = 
# dataXY = np.zeros([100,100])
# for i in range (100):
#   for j in range (100):
#     if ((i in xdata) and (j in ydata)):
#       dataXY[i][j] = -0.5
#     else:
#       dataXY[i][j] = 0.5

# X, Y = np.meshgrid(np.linspace(60, 100, 40), np.linspace(9, 48, 40))
# Z = -0.5 if ((X in xdata) and (Y in ydata)) else 0.5
# fig, ax = plt.subplots()
# ax.imshow(Z)
# plt.show()

# # create discrete colormap
# cmap = colors.ListedColormap(['white', 'black'])
# bounds = [-1,0,1]
# norm = colors.BoundaryNorm(bounds, cmap.N)
# fig, ax = plt.subplots()
# ax.imshow(dataXY, cmap=cmap, norm=norm)
# # draw gridlines
# ax.grid(which='major', axis='both', linestyle='-', color='k', linewidth=2)
# ax.set_xticks(np.arange(20, 100, 1));
# ax.set_yticks(np.arange(20, 100, 1));
# plt.show()