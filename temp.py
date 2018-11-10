# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np

# Initializing the number of levels in the hierarchy
N=5

#Initializing the number of initial stem cells
S_c=5

#Initializing the tissue
T=np.zeros((N,1))

#Inserting the stem cells
T[0,0]=S_c

#Initializing the division/mutation/differentiation rates array. The rist column is for the symmetric division, 
# the second is for symmetric differentiation and the third is cell death

R=np.zeros((N,3))

#Initializing each one of the rates

p_0=0.2
q_0=0.8
p_1=0.1
q_1=0.8
q_2=0.9

R[0,0]=p_0
R[0,1]=q_0
R[1,0]=p_1
R[1,1]=q_1
R[1,2]=1-(p_1+q_1)
R[2,1]=q_2
R[3,1]=q_2
R[4,1]=q_2
R[2,2]=1-q_2
R[3,2]=1-q_2
R[4,2]=1-q_2

N_T=np.sum(T)
print(N_T)

print(T)

A=np.ndarray.flatten(np.multiply(R,T)/N_T)

print(A)
print(np.trim_zeros(A))
print(R)