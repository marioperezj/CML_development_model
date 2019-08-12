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

directory='Histogram_alpha_0.00_beta_1.00_mu_1.0e-05_new_0.50/'
cp_time=np.zeros(1)
bc_time=np.zeros(1)
for j in range(0,4):
    for i in range(0,50):
        file=directory+'blood_'+str(i)+'_'+str(j)+'_.txt'
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
        if np.amax(data[:,0])/365.0>9.5:
            print(file)
        if np.amax(data[:,0])/365.0<9.5:
            if np.amax(bcr_abl_cells[:,17])*100>10.0:
                cp_time_var=time_array[bcr_abl_cells[:,17]*100>10.0][0]-time_array[bcr_abl_cells[:,17]*100>0.0][0]
                bc_time_var=time_array[long-1]-time_array[bcr_abl_cells[:,17]*100>10.0][0]
                cp_time=np.append(cp_time,cp_time_var)
                bc_time=np.append(bc_time,bc_time_var)
cp_time=np.delete(cp_time,0)
print(cp_time)
print(bc_time)
plt.hist(cp_time,20)
plt.hist(bc_time,20)