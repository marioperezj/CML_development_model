#Creating the plot of the cumulative average of the uneven time spaced series of the KInetic montecarlo simulation
import numpy as np
import matplotlib.pyplot as plt
raw_data=np.genfromtxt('miaus.txt',delimiter=' ')
leng=len(raw_data[:,0])
number_levels=5 #NUmber of levels in the simulation should be consistent with the output file
delta_cell=np.cumsum(np.multiply(raw_data[:,2:number_levels+2],raw_data[:,1].reshape(leng,1)),axis=0)
for i in range(1,leng):
    delta_cell[i,:]=delta_cell[i,:]/raw_data[i,0]
for i in range(0,number_levels):
    plt.scatter(raw_data[:,0],delta_cell[:,i],marker='o',s=1,c='blue')
plt.grid(True)
#plt.ylim(top=70)
#Uncomment to adjust the range values due to unwanted long fluctuations at the beggining of the plot
#plt.ylim(bottom=0)
plt.xlabel("Time")
plt.ylabel("<N>")
plt.ticklabel_format(axis='x',style='sci',scilimits=(5,5))
plt.show()
print(delta_cell[leng-1,:])
