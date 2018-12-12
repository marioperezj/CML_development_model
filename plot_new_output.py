#Creating the plot of the cumulative average of the uneven time spaced series of the KInetic montecarlo simulation
import numpy as np
import matplotlib.pyplot as plt
#raw_data=np.genfromtxt('out_file_kmc.txt',delimiter=' ')
with open('out_file_kmc.txt',"r") as f:
    data=[x.split() for x in f.readlines()]
data=np.array(data)
length=len(sorted(data,key=len, reverse=True)[0])
data_final=np.array([i+[None]*(length-len(i)) for i in data],dtype=float)
leng=len(data_final[:,0])
#print(data_final)
print(data_final[:,4:length])
print(data_final[:,1].reshape(leng,1))
delta_cell=np.nancumsum(np.multiply(data_final[:,4:length],data_final[:,1].reshape(leng,1)),axis=0)
print(delta_cell)
for i in range(1,leng):
    delta_cell[i,:]=delta_cell[i,:]/data_final[i,0]
for i in range(6,length-4):
    plt.scatter(data_final[:,0],delta_cell[:,i],marker='o',s=3,label='column'+str(i))
#plt.grid(True)
##plt.ylim(top=70)
##Uncomment to adjust the range values due to unwanted long fluctuations at the beggining of the plot
##plt.ylim(bottom=0)
#plt.xlabel("Time")
#plt.ylabel("<N>")
#plt.ticklabel_format(axis='x',style='sci',scilimits=(5,5))
plt.legend(bbox_to_anchor=(1.04,1), borderaxespad=0)
plt.show()
#print(delta_cell[leng-1,:])
