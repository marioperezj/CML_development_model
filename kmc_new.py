#This code was created by Mario Perez at ELTE University, it creates a Kinetic Monte Carlo simulation of a tissue. 
#This tissue is based on the following article https://www.nature.com/articles/ncomms14545 when the parameters definitions can be found
#The sumlation outputs a file with the following data in columns
#Time /t delta_time /t Number of cells in each level /t Level to affect /t Event type to affect
import numpy as np
import os
import time
import argparse
#Import different parameters for the simulation from the command line
parser = argparse.ArgumentParser()
parser.add_argument("-t", help="Time running simulation", type=int)
parser.add_argument("-p", help="p vaue for all the levels", type=float)
parser.add_argument("-gamma",help="Gamma values for the simulation",type=float)
parser.add_argument("-n",help="Number of levels in the simualtion",type=int)
parser.add_argument("-idfile",help="Name of the file or path to identift the file",type=str)
args = parser.parse_args()
start_time = time.time() #Measuring the initial time to test the performance of the simulation.
number_levels=args.n
rate_number=3 #Number of available events in the simulation, in this version it contains scd, scdif and adif
gamma=args.gamma
p=args.p
p_stem_cell=0.
sim_time=args.t
counter=0 #counter of the steps in the KMC
id_file=args.idfile+'.txt' #id file to print the output
try:
    os.remove(id_file) #Check if the file exists, remove it and create a new one
except OSError:
    pass
#Creating the matrix that will contains most of the data, it's organized by row (each row is a level)
#The columns are delta,p,q,rscd, rscdif, acdif
rates_matrix=np.zeros((number_levels,rate_number+3)) 
#Creating the arrays that will storage the dynamic number of cells and the defined number of cells
c_set=np.array([10,20,30,40,50])
c=np.array([10,20,30,40,50])
#Function that creates the matrix with the rates to be considered in the KMC simulation
#The function returns two items, the first one is the matrix with the rates, the second is a flag to assure that p and q are fixed properly
#The flag is zero if p and q are correct, different of zero otherwise
def construct_rates(rates_matrix_par,number_levels_par,gamma_par,p_par,p_stem_cell_par,c_set_par):
    sanity_p_q=1
    #Set number of cells in the stem cell level separately
    def set_p_value_stem_cell():
        rates_matrix_par[0,1]=p_stem_cell_par
    #Set values for delta, p and q for the boundary levels with particular values
    def set_initial_rates():
        rates_matrix_par[number_levels_par-1,0]=2.0
        rates_matrix_par[number_levels_par-2,0]=1.0
        rates_matrix_par[number_levels_par-1,1]=1.0
        rates_matrix_par[number_levels_par-1,2]=1.0
    #calculate gamma for all levels, starting from TDF level down to the stem cell level    
    def calculate_delta():
        for i in range(number_levels_par-3,-1,-1):
            rates_matrix_par[i,0]=rates_matrix_par[i+1,0]/gamma_par
    #set p value for progenitors levels        
    def set_p_values():
        for i in range(1,number_levels_par-1,1):
            rates_matrix_par[i,1]=p_par
    #Set q values for progenitors levels using delta and p
    def set_q_values():
        for i in range(1,number_levels_par-1,1):
            rates_matrix_par[i,2]=2.0*(rates_matrix_par[i-1,0]/rates_matrix_par[i,0])/rates_matrix_par[i,1]
    #Set the rate of csd accordingly to the model        
    def set_rscd():
        for i in range(0,number_levels_par,1):
            rates_matrix_par[i,3]=0.5*rates_matrix_par[i,0]*(1-rates_matrix_par[i,2])*rates_matrix_par[i,1]
    #Set the rate of csdif accordingly to the model              
    def set_rscdif():
        for i in range(0,number_levels_par,1):
            rates_matrix_par[i,4]=0.5*rates_matrix_par[i,0]*rates_matrix_par[i,1]
    #Set the rate of acdif accordingly to the model              
    def set_radif():
        for i in range(0,number_levels_par,1):
            rates_matrix_par[i,5]=rates_matrix_par[i,0]*(1-rates_matrix_par[i,1])
    #performing all the cahnges                 
    set_initial_rates()
    set_p_value_stem_cell()
    calculate_delta()
    set_p_values()
    set_q_values()
    set_rscd()
    set_rscdif()
    set_radif()
    #performing a sanity check that p and q are not greater than 1, and changing the flag to appropiate value
    if np.amax(rates_matrix_par[:,1:3])>1.0:
        sanity_p_q=0
    else:
        1
    reduced_rates=rates_matrix_par[:,3:6] #Creating the list of final reduced rates to use in the monte carlo iteration.
    #dividing by the number of cells per level, resulting in rates per cell
    for i in range(0,number_levels):
        reduced_rates[i,:]=reduced_rates[i,:]/c_set_par[i]
    return reduced_rates, sanity_p_q
#function that perform the cellular event in the level, with the outcomes of the KMC step, this function updates the global array c
def perform_event(event_to_perform_par,level_to_perform_par,number_levels_par):
    global c
    if event_to_perform_par==0:#perform scd
        c[level_to_perform_par]=c[level_to_perform_par]+1
    elif event_to_perform_par==1 and level_to_perform_par!=number_levels_par-1:#perform scdif in a progenitor level
        c[level_to_perform_par]=c[level_to_perform_par]-1
        c[level_to_perform_par+1]=c[level_to_perform_par+1]+2
    elif event_to_perform_par==1 and level_to_perform_par==number_levels_par-1:#perform scdif in the TDF level
        c[number_levels_par-1]=c[number_levels_par-1]-1
    elif event_to_perform_par==2: #Perform acdif 
        c[level_to_perform_par+1]=c[level_to_perform_par+1]+1
    else:
        print('Event not found') #Perform warning message in case of undefined event
#Function that performs one step of the KMC it receives as arguments c, and initial state the number of levels and time
#The function update the global variable c using perform_event and the time also choose one of the events using the proper algorithm
def kinetic_montecarlo_step(c_par,initial_rates_par,number_levels_par,t_par):
    initial_state=np.multiply(c_par.reshape(number_levels_par,1),initial_rates_par) #Modifiyng rates according to actual number
    rates_indices=np.where(initial_state>0.)#Selecting only non zero rates
    states_list=initial_state[rates_indices]# Getting the indices for the non zero states
    final_rates=np.cumsum(states_list)#create the cumulative sum of the rates
    max_rate=np.amax(final_rates)#getting the maximum of the sum to multiply the random number
    u=np.random.uniform(0,1,1)
    binary_result=np.searchsorted(final_rates,u*max_rate) #Performing binary search of the correct event
    level_to_perform=rates_indices[0][binary_result] #Getting the desired level to affect
    event_to_perform=rates_indices[1][binary_result] #Getting the rigth event to perform
    perform_event(event_to_perform,level_to_perform,number_levels_par) #Performing the event
    u_new=np.random.uniform(0,1,1)
    delta_t=np.log(1/u_new)/max_rate #get delta_t
    data=np.concatenate((t_par,delta_t,c,level_to_perform,event_to_perform),axis=0).reshape(1,number_levels_par+4) #CReating the proper array to print
    with open(id_file, 'a') as f:
        np.savetxt(f,data,fmt='%5.10f')#Saving the array to the file
    global t
    t=t_par+delta_t #Update the global varibale t

initial_rates,sanity_p_q=construct_rates(rates_matrix,number_levels,gamma,p,p_stem_cell,c_set) #Calling function contstruct rates
t=np.array([0.]) #Initializing t
if sanity_p_q!=0:#Creating sanity check otherwise stop simulation
    while t<sim_time: #Iterating many steps of the KMC
        kinetic_montecarlo_step(c,initial_rates,number_levels,t)
        counter+=1
else:
    print('The definitions of gamma and p are inconsistent (p or q > 1) please select new ones carefully')
print('The number of the KMC steps is:', counter)
print("--- Total time of the simulation in seconds: %s ---" % (time.time() - start_time)) #Print the total time of the simulation