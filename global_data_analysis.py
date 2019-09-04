import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(precision=12)
np.set_printoptions(suppress=True)
def readrow(row, cols): 
    data = np.fromstring(row, sep=' ')
    data.resize((cols,))
    return data
num_lev=21
num_mut=20

directory='correct_alpha_0.020_beta_0.70_mu_1.0e-05_exp9_1.20/'
cp_time=np.zeros(1)
bc_per=np.zeros(1)
total_per=np.zeros(1)
total_time=np.zeros(1)
for i in range(1,2000):
    file=directory+'blood_'+str(i)+'.txt'
    try:
        with open(file) as f:
            data=np.array([readrow(row, (num_mut+1)*num_lev+1) for row in f])
    except:
        continue
    time_array=data[:,0]/365.0
    dict_leuke=np.empty((num_lev,),dtype=object)    
    for i in range(num_lev):
        dict_leuke[i]=[k+i for k in range(1,(num_mut+1)*num_lev+1,num_lev)]
    total_leuke_blood=np.empty((num_lev,),dtype=object)
    for i in range(0,num_lev):
        total_leuke_blood[i]=np.sum(data[:,dict_leuke[i]],axis=1)
    total_wbc=np.sum(data[:,1:],axis=1)
    long=np.shape(total_leuke_blood[0])[0]
    bcr_abl_cells=np.zeros((long,num_lev))
    for k in range(0,long):
        for i in range(0,num_lev):
            bcr_abl_cells[k,i]=total_leuke_blood[i][k]/total_wbc[k]
    total_per=np.append(total_per,np.amax(bcr_abl_cells[-1,17]*100))
    bp_init=np.argmax(-np.diff(bcr_abl_cells[:,17]*100)>0)
    bc_per=np.append(bc_per,bcr_abl_cells[bp_init-10,17]*100)
    bc_per_var=bcr_abl_cells[bp_init-10,17]*100
    total_time=np.append(total_time,time_array[-1])
    if bc_per_var>50:
        print(bc_per_var,file)
        print(time_array[-1],file)
    if time_array[-1]>19.9:
        print(time_array[-1],file)
    if np.amax(bcr_abl_cells[:,17])*100>10.0:
        cp_time_var=time_array[bcr_abl_cells[:,17]*100>10.0][0]-time_array[bcr_abl_cells[:,17]*100>0.0][0]
        cp_time=np.append(cp_time,cp_time_var)
bc_per=np.delete(bc_per,0)
total_time=np.delete(total_time,0)
bins=20
plt.hist(total_time,bins)
plt.title('a')
plt.figure(1)
plt.hist(bc_per,bins,color='blue')
plt.title('Blast percentage when critical mutation appears',fontsize=20)
plt.xlabel('Time (years)',fontsize=20)
plt.ylabel('Number of patients',fontsize=20)
plt.tick_params(axis='both',labelsize='xx-large')
plt.figure(2)
cp_time=np.delete(cp_time,0)
total_per=np.delete(total_per,0)
plt.hist(cp_time,bins,color='blue')
plt.title('Duration of chronic phase',fontsize=20)
plt.xlabel('Time (years)',fontsize=20)
plt.ylabel('Number of patients',fontsize=20)
plt.tick_params(axis='both',labelsize='xx-large')
plt.figure(3)
plt.hist(total_per,bins)
plt.title('total_per')
plt.hist(cp_time,bins,color='blue')
plt.figure(4)
