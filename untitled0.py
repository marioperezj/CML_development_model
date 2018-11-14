#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 13:45:20 2018

@author: mario
"""
import numpy as np
import math as math
#Creating the first class called tissue class, this class will storage different data as the division/differentiation/death ratios.

#Creating an array with the desired properties

n=5
events=3
t=np.zeros((n*events,events))

#Initializing the variables for the rates and the numbers of cells in each iteration

for i in np.arange(n):
    for k in np.arange(1,events+1):
        t[3*i+k-1:3*i+k,1]=k
        t[3*i+k-1:3*i+k,0]=i+1
        
c=np.zeros((n,1))

c[0]=5

t[0,2]=0.1
t[1,2]=0.9
t[3,2]=0.2
t[4,2]=0.7
t[5,2]=0.1
t[6,2]=0.3
t[7,2]=0.6
t[8,2]=0.1
t[9,2]=0.4
t[10,2]=0.5
t[11,2]=0.1
t[12,2]=0.9
t[14,2]=0.1

n_t=np.sum(c)
#print(n_t)

k_m=np.copy(t)

for i in np.arange(n):
    for k in np.arange(1,events+1):
        k_m[3*i+k-1:3*i+k,2]=c[i]*t[3*i+k-1:3*i+k,2]/n_t
        
k_m=k_m[k_m[:,2].argsort()]
k_m=k_m[k_m[:,2]!=0]

s = np.random.uniform(0,1,1)

def binary_search(A,T):
    L=0
    R=np.shape(k_m)[0]-1
    print(s)
    while L < R:
        m = math.floor((L + R) / 2)
        if A[m,2] < T:
            L = m + 1
        else:
            R = m
    return print(L)

binary_search(k_m,s)

#print(c)
#print(t)
#print(k_m)