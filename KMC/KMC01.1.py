#Kinetic Monte Carlo Project Part 1 Version 1
import os
import numpy as np
import math
import random

kB = 1.380649*(10**-23)
#Make a N dimensional Cube for the Lattice
N = 100 #default N
Lattice = np.zeros([N,N,N]) #This is either 0 or 1
Q = np.zeros([N,N,N]) #stores charges
LatticeEField = np.zeros([N,N,N])
LatticeEEnergy = np.zeros([N,N,N])
LatticeEForce = np.zeros([N,N,N])
# Occupancy = np.zeros([N,N,N]) #This is either 0 or 1
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
#Fill with Atoms
    
    #Method 1
    #Inputs the empty lattice and number of atoms
    #Outputs the filled Lattice
def Initialize_Lattice(L,Q,N):
    Ni = N
    while (Ni>0):
        r1 = int(N*random())
        r2 = int(N*random())
        r3 = int(N*random())
        if (L[r1][r2][r3] == 0):
            L[r1][r2][r3] = 1
            Q[r1][r2][r3] = random()
            Ni -=1
    return L
    
    #Method 2
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
    
    #Method 1
    #Inputs the empty lattice and number of atoms
    #Outputs the filled electric field lattice values
def LatticeElectricField(N):
#     n = len(L) #L is always a cube
    for i in range (N):
        for j in range (N):
            for k in range (N):
                    LatticeEField[i][j][k] = (1/(i+1)(j+1)(k+1))**2
                    #Assume 1/r**2 field for now
    return LatticeEField
    
    #Method 2
    # Inputs an existing geometry from MD and number of Atoms
    # Outputs the electric field values for the ith position
def EField_Geom(i):
    GEFieldX = (1/i+1)**2
    GEFieldY = (1/i+1)**2
    GEFieldZ = (1/i+1)**2
    return GEFieldX, GEFieldY, GEFieldZ
        
#Define Function for energy
    
    #Method 1
    #Inputs the lattice and charges
    #Outputs the net electric field at every lattice point
def Energy_Lattice(L,Q):
    n = len(L) #always a cube
    for i in range (n):
        for j in range (n):
            for k in range (n):
                    LatticeEEnergy[i][j][k] = Q[i][j][k]*LatticeEField[i][j][k]    
    return LatticeEEnergy
    
    #Method 2
    # Inputs an existing geometry from MD and number of Atoms
    # Outputs the electric field values for the ith position
def EField_Geom(gef,i):
    GEEnergyI = GQ[i]*(gef[i][0]**2+gef[i][1]**2+gef[i][2]**2)
    return GEEnergyI

#Define Function for Force

#Define a rate/disorder function
#Define wait time parameters

#Define Direction of field at a given point from all possible directions
#Define Hop parameters

#Define iterations for any lattice point

#Run the modified Monte Carlo
















