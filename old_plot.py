#Creating the plot of the cumulative average of the uneven time spaced series of the KInetic montecarlo simulation
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
with open('test2.txt',"r") as f:
    data=[x.split() for x in f.readlines()]
length_raw=len(sorted(data,key=len, reverse=True)[0])
num_lev=21
d={}
d=['level'+str(i) for i in range(0,length_raw)]
df = pd.DataFrame.from_records(data,columns=d)
df_float_raw=df.apply(pd.to_numeric)
df_float=df_float_raw.iloc[:,0:2+9*num_lev]
long=df_float.shape[0]
length=df_float.shape[1]
print(length)
print(long)
wild_type=df_float.iloc[:,2:num_lev+2]
dic=[[] for x in range(num_lev)]
for i in range(0,num_lev):
    dic[i]=['level'+str(k) for k in range(i+num_lev+2,length,num_lev)]
final_sum=np.zeros([long,num_lev])
for i in range(0,len(dic)):
    final_sum[:,i]=df_float[dic[i]].sum(axis=1)
final_data=pd.DataFrame(final_sum)
ratio_data=np.zeros([long,num_lev])
for i in range(0,num_lev-1):
    for k in range(0,long):
        ratio_data[k,i]=final_sum[k,i]/wild_type.iloc[k,i]
ratio_data=pd.DataFrame(ratio_data)
#print(df_float.iloc[:,2+1*num_lev:19+2*num_lev])
#total_wild=wild_type.iloc[:,1:9].sum(axis=1)
#total_bcr_abl=final_data.sum(axis=1)
#print(df_float.iloc[101362,2:13])
#print(df_float.iloc[101362,2+num_lev:13+num_lev])
#print(df_float.iloc[101362,2+2*num_lev:13+2*num_lev])
total_progenitors_bcr=final_data.iloc[:,1:9].sum(axis=1)
#print(df_float.loc[df_float.iloc[:,32+21+21+21]>30])
total_progenitors_wild=wild_type.iloc[:,1:9].sum(axis=1)
prog_ratio=np.zeros((long))
#tot_term=np.zeros((long))
#tot_term=wild_type.iloc[:,25]+final_data.iloc[:,25]
for k in range(0,long):
    prog_ratio[k]=total_progenitors_bcr[k]/total_progenitors_wild[k]
#plt.scatter(df_float.iloc[:,0]/365,prog_ratio[:],marker='o',s=1,label='column'+str(i))
show=2
for i in range(2+(show-1)*num_lev,2+(show)*num_lev):
    plt.scatter(df_float.iloc[:,0]/365.0,df_float.iloc[:,i],marker='o',s=8,label='Level '+str((i-2)%num_lev))
#plt.figure()
#for i in range(2+2*num_lev,15+2*num_lev):
#    plt.scatter(df_float.iloc[:,0]/365.0,df_float.iloc[:,i],marker='o',s=10,label='Level '+str((i-2)%num_lev))
#plt.grid(True)
#plt.yscale('log')   
##plt.ylim(0,1)
#plt.xlabel("Time (years)")
#plt.ylabel("BCR-ABL/BCR for terminal level")
#plt.legend(bbox_to_anchor=(1.04,1), borderaxespad=0)
#plt.figure()
#for i in range(0,num_lev):
#    plt.scatter(df_float.iloc[:,0]/365.0,ratio_data.iloc[:,i],marker='o',s=10,label='Level '+str(i))
#for i in range(0,num_lev):
#    plt.scatter(df_float.iloc[:,0]/365,ratio_data.iloc[:,i],marker='o',s=5,label='column'+str(i))
plt.yscale('log')  
plt.ylim(1,1e14)
#plt.xlim(0,8)
plt.tick_params(axis='both',labelsize='x-large')
#plt.xlim(2,8)
#plt.plot(df_float.iloc[:,0]/365,prog_ratio,marker='.')
#plt.axvline(x=0.38,label='first mutation',color='r',linewidth=3.0)
#plt.axvline(x=0.99,label='second mutation',color='g',linewidth=3.0)
#plt.axvline(x=7.57,label='third mutation',color='y',linewidth=3.0)
#plt.scatter(df_float.iloc[:,0]/365,total_progenitors_bcr,marker='o',s=1,label='column'+str(i))
plt.grid(True)
#Uncomment to adjust the range values due to unwanted long fluctuations at the beggining of the plot 
plt.xlabel("Time (years)",fontsize=15)
#plt.ylabel("Ratio of CML over wild type cells \n in progenitor levels",fontsize=15)
plt.ylabel("Number of CML cells with three mutations",fontsize=15)
#plt.legend(loc=9,prop={'size': 15})
plt.legend(bbox_to_anchor=(1.04,1), borderaxespad=0, fontsize=11)
plt.show()
