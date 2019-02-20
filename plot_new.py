import numpy as np
import matplotlib.pyplot as plt
import io
def readrow(row, cols):
    data = np.fromstring(row, sep=' ')
    data.resize((cols,))
    return data
#final_data=np.zeros((200,3))
#i=0
#while i<200:
#    file='simu'+str(i)+'.txt'
#    with open(file, 'rb') as f:
#        data = np.array([readrow(row, 85) for row in f])
#    final_data[i,0]=i
#    if np.size(data[data[:,72]>1e3])<1:
#        print('bad')
#    else:
#        final_data[i,2]=data[data[:,72]>1e3][0,0]/365.0
#    if np.size(data[data[:,21]>1e7])<1:
#        print('very_bad')
#    else:
#        final_data[i,1]=data[data[:,21]>1e7][0,0]/365.0
#    #print(final_data[i,:])
#    i=i+1
#with open('results.txt', 'a') as f:
#    np.savetxt(f,final_data,fmt='%5.10f')#Saving the array to the file
with open('simu54.txt', 'rb') as f:
    data = np.array([readrow(row, 85) for row in f])
with open('simu15_3.txt', 'rb') as f:
    data1 = np.array([readrow(row, 85) for row in f])
with open('simu56_3.txt', 'rb') as f:
    data2 = np.array([readrow(row, 85) for row in f])
final_data=np.sum(data[:,1:86],axis=1)/83274066388.18518
final_data1=np.sum(data1[:,1:86],axis=1)/83274066388.18518
final_data2=np.sum(data2[:,1:86],axis=1)/83274066388.18518
plt.yscale('log')
plt.grid(True)
plt.xlabel("Time (years)",fontsize=20)
plt.ylabel("Ratio BCR-ABL/ABL",fontsize=20)
plt.scatter(data1[:,0]/365.0,final_data1,marker='o',s=10,label='patient 1')
plt.scatter(data2[:,0]/365.0,final_data2,marker='o',s=10,label='patient 2')
plt.scatter(data[:,0]/365.0,final_data,marker='o',s=10,label='patient 3')
plt.legend(bbox_to_anchor=(0.3,0.8), borderaxespad=0, fontsize=17)
plt.tick_params(axis='both',labelsize='xx-large')
#with open('simu54.txt', 'rb') as f:
#    data = np.array([readrow(row, 85) for row in f])
#first_mutant=np.sum(data[:,1:10],axis=1)
#two_mutant=np.sum(data[:,22:31],axis=1)
#third_mutant=np.sum(data[:,43:52],axis=1)
#four_mutant=np.sum(data[:,64:73],axis=1)
#plt.yscale('log')
#plt.ylim(1,1e6)
#plt.grid(True)
#plt.scatter(data[:,0]/365.0,first_mutant,marker='o',s=10,label=str(0)+' mutations')
#plt.scatter(data[:,0]/365.0,two_mutant,marker='o',s=10,label=str(1)+' mutations')
#plt.scatter(data[:,0]/365.0,third_mutant,marker='o',s=10,label=str(2)+' mutations')
#plt.scatter(data[:,0]/365.0,four_mutant,marker='o',s=10,label=str(3)+' mutations')
#plt.legend(bbox_to_anchor=(0.3,0.8), borderaxespad=0, fontsize=17)
#plt.tick_params(axis='both',labelsize='xx-large')
#plt.xlabel("Time (years)",fontsize=20)
#plt.ylabel("Number of BCR-ABL cells",fontsize=20)
