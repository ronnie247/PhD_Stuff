#Kinetic Monte Carlo Project Part 3 Version 1
import os
import numpy as np
import math
import random
import pandas as pd
import matplotlib.pyplot as plt
# import matplotlib.animation as animate
# from matplotlib.animation import writers
# from mpl_toolkits import mplot3d
import plotly.express as px

# Make Lattice
NL = 1000
NA = 0.1*NL

Lattice = np.zeros([NL,NL])
energy_term = np.zeros([NL,NL])

# Set Temperature
kbT = 1.380649*(10**-23) * 273

# Occupy 10% of them
temp_NA = NA
for i in range (NL):
    for j in range (NL):
        # for k in range (NL):
            if (temp_NA > 0):
                r1 = random.random()
                if (r1 < (NA/NL)):
                    Lattice[i][j] = 1
                    temp_NA -= 1
                else :
                    Lattice[i][j] = 1
            energy_term[i][j] = kbT * random.uniform(-1.0,1.0)
# energy_term = np.random.normal(kbT, (kbT/2), (NL**3)) + np.random.normal(kbT, (kbT/2), (NL**3)) + np.random.normal(kbT, (kbT/2), (NL**3))
#assign energy term  at each point in normal distribution
#np.random.normal(mu, sigma, 1000)
#here mu = kbT

# only one particle moves, starting from origin
# x,y,z coordinates
pos_i = np.zeros((2,), dtype=int) #initial position = (1,1,1)
pos_i[0] = 1
pos_i[1] = 1
# pos_i[2] = 1
pos_j = np.zeros((2,), dtype=int) #position at any time t
pos_j[0] = NL/2
pos_j[1] = NL/2
# pos_j[2] = NL/2
pos_f = np.zeros((2,), dtype=int) #final position
total_x2 = 0
D = 0
t = 0

NSteps = 100
pos_traj = np.zeros([2,NSteps])
pos_traj[0][0] = pos_j[0]
pos_traj[1][0] = pos_j[1]
# pos_traj[2][0] = pos_j[2]
energy_direc = np.zeros(4)
rate_direc = np.zeros(4)
#6 possible directions - in order -x,+x,-y,+y,-z,+z

# Declare all flags
flag1 = 0 # The direction with minimum energy is stored in flag1 value [1,6]
flag2 = 0 # flag2 stores whether the movement will take place value 0 or 1
flag3 = flag2 # stores previous movement
flag_sum = 0 #Total number of movements

Lattice[pos_j[0]][pos_j[1]] = 2
# calculate energy change acting at that point in each direction
for i in range (1, NSteps):
    #if x or y or z is zero, then -x,-y,-z energy = infinite
    if (pos_j[0] <= 1):
        energy_direc[0] = math.inf
    if (pos_j[1] <= 1):
        energy_direc[2] = math.inf
    # if (pos_j[2] <= 1):
    #     energy_direc[4] = math.inf
    #if x or y or z is NL, then +x,+y,+z energy = infinite
    if (pos_j[0] >= NL-2):
        energy_direc[1] = math.inf
    if (pos_j[1] >= NL-2):
        energy_direc[3] = math.inf
    # if (pos_j[2] >= NL-2):
        # energy_direc[5] = math.inf
    # in place of this you can add periodic boundary conditions where
    # pos[j] = pos[j] % NL
    # Energy change in any direction is energy final - initial
    # not allowed to bump into atoms
    if (Lattice[pos_j[0]-1][pos_j[1]] != 1) :
        energy_direc[0] = energy_term[pos_j[0]-1][pos_j[1]] - energy_term[pos_j[0]][pos_j[1]]
    else :
        energy_direc[0] = math.inf
    if (Lattice[pos_j[0]+1][pos_j[1]] != 1) :
        energy_direc[1] = energy_term[pos_j[0]+1][pos_j[1]] - energy_term[pos_j[0]][pos_j[1]]
    else :
        energy_direc[1] = math.inf
    if (Lattice[pos_j[0]][pos_j[1]-1] != 1) :
        energy_direc[2] = energy_term[pos_j[0]][pos_j[1]-1] - energy_term[pos_j[0]][pos_j[1]]
    else :
        energy_direc[2] = math.inf
    if (Lattice[pos_j[0]][pos_j[1]+1] != 1) :
        energy_direc[3] = energy_term[pos_j[0]][pos_j[1]+1] - energy_term[pos_j[0]][pos_j[1]]
    else :
        energy_direc[3] = math.inf
    # if (Lattice[pos_j[0]][pos_j[1]][pos_j[2]-1] != 1) :
    #     energy_direc[4] = energy_term[pos_j[0]][pos_j[1]][pos_j[2]-1] - energy_term[pos_j[0]][pos_j[1]][pos_j[2]]
    # else :
    #     energy_direc[4] = math.inf
    # if (Lattice[pos_j[0]][pos_j[1]][pos_j[2]+1] != 1) :
    #     energy_direc[5] = energy_term[pos_j[0]][pos_j[1]][pos_j[2]+1] - energy_term[pos_j[0]][pos_j[1]][pos_j[2]]
    # else :
    #     energy_direc[5] = math.inf
    temp_min = energy_direc[0]
    
    # The direction with minimum energy is stored in flag1
    for j in range (4):
        if (energy_direc[j] < temp_min):
            flag1 = j
    # see if that direction is favoured, else other direction
    # flag2 stores whether the movement will take place
    #accepted move if energy decreases and lattice site = unocc
    if (energy_direc[flag1] < 0):
        flag2 = 1
        flag_sum +=1
    #accepted with probability if energy increases
    else :
        r2 = random.random()
        r3 = math.exp((-1*kbT))
        if (r2 < r3):
            flag2 = 1
            flag_sum +=1
        else : 
            flag2 = 0
    if (flag2 == 1): #if movement is allowed
        #and if the previous movement wasn't the same
        if ((flag1 == 0) and (flag2 !=flag3)):
            pos_j[0] -=1
        if ((flag1 == 1) and (flag2 !=flag3)):
            pos_j[0] +=1
        if ((flag1 == 2) and (flag2 !=flag3)):
            pos_j[1] -=1
        if ((flag1 == 3) and (flag2 !=flag3)):
            pos_j[1] +=1
        # if ((flag1 == 4) and (flag2 !=flag3)):
            # pos_j[2] -=1
        # if ((flag1 == 5) and (flag2 !=flag3)):
            # pos_j[2] +=1
        flag3 = -1*flag2
    Lattice[pos_j[0]][pos_j[1]] = 2
    for j in range (2):
        pos_traj[j][i] = pos_j[j]
    r4 = random.random()
    t -= (math.log(r4)/flag_sum)
# calculate time from rate formula, and distance
x2 = (pos_j[0]-pos_i[0])**2
y2 = (pos_j[1]-pos_i[1])**2
# z2 = (pos_j[2]-pos_i[2])**2
total_x2 = x2+y2
D = (total_x2/4)/t

# print (total_x2)
# print(t)
# print (D)
# for i in range (NSteps):
    # print(pos_traj[0][i])

f = open("KMC03.1.0.out","w")
f.write("\n KMC03.1.0.py Output")
f.write("\n Total x2:")
f.write(str(total_x2))
f.write("\n D:")
f.write(str(D))
f.write("\n T:")
f.write(str(t))
f.write("\n Trajectory:")
# str_traj = str(pos_traj)
for i in range (NSteps):
    f.write(str(pos_traj[0][i]))
    f.write("\t")
    f.write(str(pos_traj[0][i]))
    f.write("\t")
    # f.write(str(pos_traj[i][0]))
    f.write("\n")
f.close()

# np.savetxt('KMC03.1.1.out', pos_traj, delimiter='\t')

#Create Pandas Dataframe for Lattice and Trajectory
# Lattice_df = pd.DataFrame(data=Lattice)
# Trajectory_df = pd.DataFrame(data=pos_traj)

# plot trajectory by plotly
# L_fig = px.scatter_3d(Lattice_df, x='x', y='y', z='z', color='blue' if Lattice_df == 1 else 'green')
# T_fig = px.scatter_3d(Trajectory_df, x='x', y='y', z='z', color='red')
# T_fig = px.scatter_3d(Trajectory_df)
# L_fig.show()
# T_fig.show()

# # Data to plot:
# # Lattice as 0 = blue, 1 = green, 2 = Red
# fig1 = plt.figure()
# ax1 = plt.axes(projection='3d')
# # Data for a three-dimensional line
# xline = pos_traj[0]
# yline = pos_traj[1]
# zline = pos_traj[2]
# ax1.plot3D(xline, yline, zline, 'red')
# # Data for three-dimensional scattered points
# coord_0 = []
# coord_1 = []
# coord_2 = []
# for x in range (NSteps):
#     for y in range (NSteps):
#         for z in range (NSteps):
#             if (Lattice[x][y][z]==0):
#                 coord_0.append([x,y,z])
#             if (Lattice[x][y][z]==1):
#                 coord_1.append([x,y,z])
#             if (Lattice[x][y][z]==2):
#                 coord_2.append([x,y,z])
# fig2 = plt.figure()
# ax2 = plt.axes(projection='3d')
# for w in range (len(coord_0)):
#     ax2.plot3D(coord_0[w][0], coord_0[w][1], coord_0[w][2], color = 'green')
# for w in range (len(coord_1)):
#     ax2.plot3D(coord_1[w][0], coord_1[w][1], coord_1[w][2], color = 'blue')
# for w in range (len(coord_2)):
#     ax2.plot3D(coord_2[w][0], coord_2[w][1], coord_2[w][2], color = 'red')
# # plt.savefig('KMC03.1.png')
# plt.show()
# # save graph