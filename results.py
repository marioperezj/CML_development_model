import numpy as np
import matplotlib.pyplot as plt
data1=np.loadtxt('total_31.txt')
data2=np.loadtxt('total_33.txt')
data3=np.loadtxt('total_35.txt')
data1[:,1]=data1[:,1]/83274066388.18518
data2[:,1]=data2[:,1]/83274066388.18518
data3[:,1]=data3[:,1]/83274066388.18518
#data=np.loadtxt('results_3_e-5.txt')
#data=np.loadtxt('results.txt')
#data1=np.loadtxt('results2.txt')
##print(data1)
#x=data1[:,2]
#xx=data[:,3]
#xxx=np.concatenate([x,xx])
#xxx=xxx[(xxx>0.)]
#xxx=xxx[xxx<20.0]
#print(xxx)
#print(np.size(xxx))
#bins=np.arange(1,20)
#plt.axis([1,20,0,25])
#ticks=np.arange(1,20,1)
#plt.xticks(ticks)
#plt.hist((xxx), bins=bins)
#data=np.loadtxt('prog_29.txt')
#print(data)
#plt.ylim(1,1e8)
plt.yscale('log')
plt.grid(True)
#plt.xlabel("Time in chronic phase",fontsize=20)
#plt.ylabel("Number of simulations ",fontsize=20)
plt.tick_params(axis='both',labelsize='xx-large')
#for i in range(1,4):
#    print(i)
#    plt.scatter(data[:,0]/365.0,data[:,i],marker='o',s=10,label=str(i-1)+' mutations')
#plt.legend(bbox_to_anchor=(0.4,1), borderaxespad=0, fontsize=17)
##plt.scatter(data[:,0],data[:,2],marker='o',s=150,label=str(i))
plt.xlabel("Time (years)",fontsize=20)
plt.ylabel("Ratio BCR-ABL/ABL",fontsize=20)
plt.scatter(data1[:,0]/365.0,data1[:,1],marker='o',s=10,label='patient 1')
plt.scatter(data2[:,0]/365.0,data2[:,1],marker='o',s=10,label='patient 2')
plt.scatter(data3[:,0]/365.0,data3[:,1],marker='o',s=10,label='patient 3')
plt.legend(bbox_to_anchor=(0.1,0.8), borderaxespad=0, fontsize=15)