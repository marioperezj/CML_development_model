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
try:
    os.remove('mean.txt')
except OSError:
    pass
cou=0 #This variable will be a counter of how many iteratios are in the simulation.
l=0 # This is a loop that allow us to run the same simulation many times to get the necessary statistics
#while l<1: #loop that creates many instances of the iteration
#    l=l+1 #updating the index of the loop
gamma=4
n=3  #NUmber of levels in the hierarchy
events=6 #Number of different events happening in the hierarchy for instance symetric cell division, death, mutation,etc.
t=np.zeros((n*events,3)) #Initializing the static array that will storage the levels and distinct rates of each one

for i in np.arange(n): #This procedure creates the array in first column is the number of level
    for k in np.arange(1,events+1): #In the second column it will be the number of rate for each level
        t[events*i+k-1:events*i+k,1]=k #The third column is the value of each one of the rates
        t[events*i+k-1:events*i+k,0]=i+1
        
c=np.zeros((n,2)) #Initializing the array that storage that will storage the number of cells in each level
c[0,0]=10000 #Also initiaties the number of cells in each level, also the total differentiation rate per level
c[1,0]=1000
c[2,0]=0

#for i in np.arange(1,n):
#    c[i,0]=0

t[3,2]=1.0 #initializing the different rates for each one of the levels in the hierarchy
t[n*events-3,2]=1
t[n*events-2,2]=1
t[n*events-1,2]=2.0
t[(n-1)*events-1,2]=1.0

for i in np.arange(2,n):
    t[((n-i)*events)-1,2]=t[(n-(i-1))*events-1,2]/gamma
    
for i in np.arange(2,n):
    t[(i)*events-3,2]=1.0
    
for i in np.arange(1,n-1):
    t[(i)*events+4,2]=(2*t[(i-1)*events+5,2])/(t[(i)*events+5,2]*t[(i)*events+3,2])
#    t[(i+1)*events-3,2]=(2)*(t[(i)*events-1,2]/t[(i+1)*events-1,2])/t[(i)*events-2,2]
    
for i in np.arange(0,n):
    t[i*events+1,2]=0.5*t[(i)*events+5,2]*t[i*events+3,2]
    t[i*events,2]=0.5*t[(i)*events+5,2]*t[i*events+3,2]*(1-t[(i)*events+4,2])
    
x=1. #Initializing the index for the Monte Carlo algorithm.
mean=c[:,0].reshape(1,n)
delta_t=0
print(t[(t[:,1]==1)|(t[:,1]==2)|(t[:,1]==4)|(t[:,1]==5)])
while x<100000: #Stop the Kinetic Monte Carlo algorithm after a definite time. 
    pri_array=np.append(np.array([x]).reshape(1,1),c[:,0].reshape(1,n),axis=1)
    mean_final=np.append(np.array([x]).reshape(1,1),mean,axis=1)
    
    with open('out.txt', 'a') as f: #Open the file to save the results
        np.savetxt(f,pri_array.reshape(1,n+1),fmt='%5.10f') #Saving the results as a table with 17 postions and numbers as integers
    with open('mean.txt','a') as g:
        np.savetxt(g,mean_final,fmt='%5.10f')
    
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
            k_m[events*i+k-1:events*i+k,2]=c[i,1]*t[events*i+k-1:events*i+k,2]#The third column is the value of each one of the rates
    
    k_m=k_m[k_m[:,2]!=0] #Choosing only non zero transition rate for different states    
    k_m=k_m[(k_m[:,1]==1)|(k_m[:,1]==2)]
            
    cum_sum=np.cumsum(k_m[:,2],axis=0).reshape(np.shape(k_m)[0],1)
       
    k_m=np.append(k_m,cum_sum,axis=1)
    
#    print(k_m)
    
    max_list=np.amax(k_m[:,3],axis=0)
                
    def binary_search(A,T): #Definition of the function performing a binary search
        L=0                #This function returns the leftmost item of the list ov events
        R=np.shape(k_m)[0]-1
        while L < R:
            m = math.floor((L + R) / 2)
            if A[int(m),3] < T:
                L = m + 1
            else:
                R = m
        return A[int(L),:]
    
    step_state=max_list*np.random.uniform(0,1,1)    
#    print(step_state)
#    print(max_list)
    event=binary_search(k_m,step_state) #Applying the binary search in the dynamic array k_m, passing a random number
    
#    print(event)
       
    if (event[1]==2) and (event[0]!=n):#Once we have the rigth state to advance we perform the actual change in the number of cells in each level
        c[event[0].astype(int),0]=c[event[0].astype(int),0]+2 #The rate number 2 is the symmetric differentiation
        c[event[0].astype(int)-1,0]=c[event[0].astype(int)-1,0]-1
    elif event[1]==1:
        c[event[0].astype(int)-1,0]=c[event[0].astype(int)-1,0]+1 #The rate number 1 is the symmetric cell division
    elif event[1]==3:
        c[event[0].astype(int)-1,0]=c[event[0].astype(int)-1,0]-1 #The rate number 3 is cell death
    elif event[1]==2 and event[0]==n:
        c[n-1,0]=c[n-1,0]-1
    else:
        print('Error event not found in list') #Error message for non defined rate of division
        
    ran=np.random.uniform(0,1,1) #Generating a new random number for updating the time in the KMC algorithm
    delta_t=np.log(1/ran)/max_list
    mean=(np.add(delta_t*c[:,0].reshape(1,n),x*mean))/(x+delta_t)
    x=x+delta_t #Computing the udpate in time

print(cou) #Print in the terminal the total number of iterations in the simulation
print("--- %s seconds ---" % (time.time() - start_time)) #Print the total time of the simulation
#print(c)
#print(t)
#print(k_m)
