#Creating the plot of the cumulative average of the uneven time spaced series of the KInetic montecarlo simulation
import numpy as np
import matplotlib.pyplot as plt
import math

num_lev=21
num_mut=31
def readrow(row, cols):
    data = np.fromstring(row, sep=' ')
    data.resize((cols,))
    return data
dummy='0_0_'
id_mu='1.0e-05'
alpha='0.00'
beta='1.00'
#id_mu='mu_8.1e-06_'
id_epsilon='0.50'
#='exp_1.50'
id_rho='.'
directory_name='Histogram_'+'alpha_'+alpha+'_beta_'+beta+'_mu_'+id_mu+'_new_'+id_epsilon+'/'
file_bm=directory_name+'bone_marrow_'+dummy
file_blood=directory_name+'blood_'+dummy
file_delta=directory_name+'delta_'+dummy
file_ratios=directory_name+'ratios_'+dummy
file=file_delta
num_col=4
beta=2
leak_par=5e-9
bm_original_cellularity=10485750000
mod_cellularity=bm_original_cellularity/0.5
x_axis=np.logspace(8,12,num=10000)
leaking_function=np.zeros_like(x_axis)
for i in np.arange(0,np.size(x_axis)):
    leaking_function[i]=1/(1+math.exp(-leak_par*(x_axis[i]-mod_cellularity)))
rho=np.zeros_like(x_axis)
for i in np.arange(0,np.size(x_axis)):
    if x_axis[i]-bm_original_cellularity<=0:
        rho[i]=1
    else:
        rho[i]=np.exp(-0.5*beta*((x_axis[i]-bm_original_cellularity)/bm_original_cellularity))
with open(file+'.txt', 'rb') as f:
    data = np.array([readrow(row, num_col) for row in f])
time=data[:,0]/365.0
#print(data[:,2])
def delta():
    plt.figure(1)
#    plt.ylim(1e-1,1e1)
#    plt.yscale('log')
#    plt.ylim(2000,3000)
#    plt.scatter(time,data[:,1],label='delta')
    plt.scatter(time,data[:,1],label='leaking')
#    plt.scatter(time,data[:,2],label='cells')
    plt.legend(loc='best', borderaxespad=0, fontsize=10)
    plt.show()
    
def cells():
    plt.figure(3)
    plt.scatter(time,data[:,2],label='cells')
    plt.legend(loc='best', borderaxespad=0, fontsize=10)
    plt.show()
    
def analytic():
#    plt.ylim(0,1)
    plt.grid('True')
#    plt.xscale('log')
#    plt.yscale('log')
    plt.scatter(x_axis,leaking_function)
#    plt.scatter(x_axis,rho)
    plt.show()
    
#cells()
#delta()
#analytic()
##################################################################################################
with open(file_ratios+'.txt', 'rb') as f:
    data_ratios = np.array([readrow(row, (num_mut+1)*num_lev+1) for row in f])

def ratios(lev_initiate,lev_final,num_mutation):
    plt.figure(2)
    for i in range(lev_initiate,lev_final):
        plt.scatter(time,data_ratios[:,num_mutation*num_lev+i+1],label=str(i%num_lev))
    plt.grid('True',which='both')
    plt.show()
    plt.ylim(-5,2)
    plt.legend(bbox_to_anchor=(1.0,1.0), borderaxespad=0, fontsize=10)
ratios(16,num_lev,5)