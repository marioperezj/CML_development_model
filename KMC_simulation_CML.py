#Initializing the necessary Python packages for different features
#The numerical package is loaded for the use of arrays and math for functions, os is loaded 
#for the use of files to write on it, time is done to keep track of the performance
import numpy as np
import math as math
import os
import time
#Initializing the time to track the performance time of the code
start_time = time.time()
#Trying to remove the file out if it exists, otherwise pass an exception. The file out.tx will storage the output of the simulations
#The output of the simulations will be at this moment just the total number of cells in each level.
try:
    os.remove('out.txt')
except OSError:
    pass
cou=0 #This variable will be a counter of how many iteratios are in the simulation.
l=0 # This is a loop that allow us to run the same simulation many times to get the necessary statistics
#while l<1: #loop that creates many instances of the iteration
#    l=l+1 #updating the index of the loop
n=5  #NUmber of levels in the hierarchy
events=3 #Number of different events happening in the hierarchy for instance symetric cell division, death, mutation,etc.
t=np.zeros((n*events,events)) #Initializing the static array that will storage the levels and distinct rates of each one

for i in np.arange(n): #This procedure creates the array in first column is the number of level
    for k in np.arange(1,events+1): #In the second column it will be the number of rate for each level
        t[3*i+k-1:3*i+k,1]=k #The third column is the value of each one of the rates
        t[3*i+k-1:3*i+k,0]=i+1
        
c=np.zeros((n,2)) #Initializing the array that storage that will storage the number of cells in each level
c[0,0]=10 #Also initiaties the number of cells in each level, also the total differentiation rate per level

t[0,2]=0.5 #initializing the different rates for each one of the levels in the hierarchy.
t[1,2]=0.5
t[3,2]=0.5
t[4,2]=0.5
t[5,2]=0
t[6,2]=0.5
t[7,2]=0.5 
t[8,2]=0
t[9,2]=0.5
t[10,2]=0.5
t[11,2]=0
t[12,2]=0
t[14,2]=1

x=0 #Initializing the index for the Monte Carlo algorithm.
while x<100: #Stop the Kinetic Monte Carlo algorithm after a definite time. 
    
    if np.amax(c)==0: #Computing the sum of the total number of cells WRONG!!!
        break
    else:
        n_cells=np.amax(c)
         
    cou=cou+1 #Updating the counter for the number of total iterations in the KMC algorithm.
    k_m=np.copy(t) #Creating a dynamic array for the KMC simulation that will storage the different rates in each iteration
       
    for i in np.arange(n):
        if c[i,0]==0:
            c[i,1]=0
        else:
            c[i,1]=1
         
    for i in np.arange(n):#This procedure creates the array in first column is the number of level THIS IS FOR THE DYNAMIC ARRAY 
        for k in np.arange(1,events+1): #In the second column it will be the number of rate for each level CHANGING IN EACH ITERATION
            k_m[3*i+k-1:3*i+k,2]=(c[i,1]*t[3*i+k-1:3*i+k,2]) #The third column is the value of each one of the rates
    
    k_m=k_m[k_m[:,2]!=0] #Choosing only non zero transition rate for different states    
            
    cum_sum=np.cumsum(k_m[:,2],axis=0).reshape(np.shape(k_m)[0],1)
    
    k_m=np.append(k_m,cum_sum,axis=1)
    
    max_list=np.amax(k_m[:,3],axis=0)
               
    def binary_search(A,T): #Definition of the function performing a binary search
        L=0                #This function returns the leftmost item of the list ov events
        R=np.shape(k_m)[0]-1
        while L < R:
            m = math.floor((L + R) / 2)
            if A[m,3] < T:
                L = m + 1
            else:
                R = m
        return A[L,:]
    
    event=binary_search(k_m,max_list*np.random.uniform(0,1,1)) #Applying the binary search in the dynamic array k_m, passing a random number
       
    if event[1]==2:#Once we have the rigth state to advance we perform the actual change in the number of cells in each level
        c[event[0].astype(int),0]=c[event[0].astype(int),0]+2 #The rate number 2 is the symmetric differentiation
        c[event[0].astype(int)-1,0]=c[event[0].astype(int)-1,0]-1
    elif event[1]==1:
        c[event[0].astype(int)-1,0]=c[event[0].astype(int)-1,0]+1 #The rate number 1 is the symmetric cell division
    elif event[1]==3:
        c[event[0].astype(int)-1,0]=c[event[0].astype(int)-1,0]-1 #The rate number 3 is cell death
    else:
        print('Error event not found in list') #Error message for non defined rate of division
        
    ran=np.random.uniform(0,1,1) #Generating a new random number for updating the time in the KMC algorithm
    x=x+np.log(1/ran)/max_list #Computing the udpate in time
    with open('out.txt', 'a') as f: #Open the file to save the results
                np.savetxt(f,c[:,0].reshape(1,n),fmt='%-17d') #Saving the results as a table with 17 postions and numbers as integers
print(cou) #Print in the terminal the total number of iterations in the simulation
print("--- %s seconds ---" % (time.time() - start_time)) #Print the total time of the simulation
#print(c)
#print(t)
#print(k_m)