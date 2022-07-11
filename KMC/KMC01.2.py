#Kinetic Monte Carlo Project Part 1 Version 2
import os
import numpy as np
import math

kB = 1.380649*(10**-23)

#Fill the lattice with N atoms
NAtoms = 1 #To begin with
#Input the coordinates for the Atoms
#This can also be used as a lattice point
Geometry = np.zeros([NAtoms,3])
GQ = np.zeros([NAtoms])
GEField = np.zeros([NAtoms,3])
GEEnergy = np.zeros([NAtoms])
GEForce = np.zeros([NAtoms,3])

#Read any file to get matrix
def file_read(file_name):
    input_file=open(file_name)
    file_content=input_file.readlines()
    input_file.close()
    R = 0
    C = 0
    temp_mat=[]
    for line in file_content:
        v_line=line.rstrip() #only remove the whitespaces after the line
        #f=v_line.split()
        if (len(v_line)>0):  #take in only non zero lines, 
            #in case any whitespace was there
            #it would have been removed by rstrip
                temp_mat.append(v_line.split())
    R=len(temp_mat)
    C=len(temp_mat[0])
    A = np.zeros([R,C])
    for i in range (R):
        for j in range (C):
            A[i][j]=np.double(temp_mat[i][j])
    return A

#Convert to compound indices
def CompIndices(a,b):
    ab = a*(a+1)/(2+b)
    return ab

#Define Lattice environment
# Inputs an existing geometry from MD and number of Atoms
# Outputs the filled lattice
# This step is redundant, can directly get a numpy lattice from file_read
def Init_MD_Lattice(L,q,N):
    for i in range (N):
        Geometry[i][0] = L[i][0]
        Geometry[i][1] = L[i][1]
        Geometry[i][2] = L[i][2]
        GQ[i] = q[i]
    return Geometry

#Define external field for ith atom
# Inputs an existing geometry from MD and number of Atoms
# Outputs the electric field values for the ith position
def EField_Geom(i):
    GEFieldX = (1/i+1)**2
    GEFieldY = (1/i+1)**2
    GEFieldZ = (1/i+1)**2
    return GEFieldX, GEFieldY, GEFieldZ

#Define Function for energy
# Inputs an existing geometry from MD and number of Atoms
# Outputs the electric field values for the ith position
def EField_Geom(gefi,i):
    GEEnergyI = GQ[i]*(gef[i][0]**2+gef[i][1]**2+gef[i][2]**2)
    return GEEnergyI

#Define Function for Force

#Define a rate/disorder function
#Define wait time parameters

#Define Direction of field at a given point from all possible directions
#Define Hop parameters

#Define iterations for any lattice point

#Run the modified Monte Carlo