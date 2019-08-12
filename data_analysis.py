import numpy as np
import matplotlib.pyplot as plt

def readrow(row, cols):
    data = np.fromstring(row, sep=' ')
    data.resize((cols,))
    return data

num_lev=21
num_mut=30
num_col=3

#='result_lam_0.7_1_'
dummy='16'
id_mu='8.0e-06'
alpha='0.00'
beta='1.00'
#id_mu='mu_8.1e-06_'
id_epsilon='0.80'
#='exp_1.50'
id_rho='.'
directory_name='final_'+'alpha_'+alpha+'_beta_'+beta+'_mu_'+id_mu+'_modified_'+id_epsilon+'/'
file_bm=directory_name+'bone_marrow_'+dummy
file_blood=directory_name+'blood_'+dummy
file_delta=directory_name+'delta_'+dummy
file_ratios=directory_name+'ratios_'+dummy
blast_level=num_lev-4

with open(file_bm+'.txt', 'rb') as f:
    data = np.array([readrow(row, (num_mut+1)*num_lev+1) for row in f])
with open(file_blood+'.txt', 'rb') as f:
    data1 = np.array([readrow(row, (num_mut+1)*num_lev+1) for row in f])
with open(file_delta+'.txt', 'rb') as f:
    data_delta = np.array([readrow(row, num_col) for row in f])
time=data[:,0]/365.0

dict_leuke=np.empty((num_lev,),dtype=object)    
for i in range(num_lev):
    dict_leuke[i]=[k+i for k in range(num_lev+1,(num_mut+1)*num_lev+1,num_lev)]
total_leuke_blood=np.empty((num_lev,),dtype=object)
for i in range(0,num_lev):
    total_leuke_blood[i]=np.sum(data1[:,dict_leuke[i]],axis=1)
total_leuke_bone=np.empty((num_lev,),dtype=object)
for i in range(0,num_lev):
    total_leuke_bone[i]=np.sum(data[:,dict_leuke[i]],axis=1)
long=np.shape(total_leuke_blood[0])[0]
total_wbc_bone=np.sum(data[:,1:],axis=1)
total_wbc=np.sum(data1[:,1:],axis=1)
bcr_abl_cells=np.zeros((long,num_lev))
bcr_abl_cells_bone=np.zeros((long,num_lev))
bcr_abl_percentage=np.sum(data1[:,num_lev+1:],axis=1)/total_wbc
bcr_abl_percentage_bone=np.sum(data[:,num_lev+1:],axis=1)/total_wbc_bone
for k in range(0,long):
    for i in range(0,num_lev):
        bcr_abl_cells[k,i]=(total_leuke_blood[i][k]+data1[k,i+1])/total_wbc[k]
for k in range(0,long):
    for i in range(0,num_lev):
        bcr_abl_cells_bone[k,i]=(total_leuke_bone[i][k]+data[k,i+1])/total_wbc_bone[k]
tot_percentage_bcr_abl=np.sum(bcr_abl_cells[:,:],axis=1)
tot_percentage_bone=np.sum(bcr_abl_cells_bone[:,:],axis=1)

def graph_cells(lev_initiate,lev_final,num_mutation):
#    plt.figure(figsize=(16,8))
    plt.figure(1)
    plt.yscale('log')
#    plt.yscale('linear')
#    plt.xlim(2.5,4.0)
    plt.ylim(1e-1,1e11)
    for i in range(lev_initiate,lev_final):
        plt.scatter(time,data[:,num_mutation*num_lev+i+1],marker=',',s=10.0,label='Level'+str(i%num_lev))
    plt.legend(bbox_to_anchor=(1.0,1.0), borderaxespad=0, fontsize=20)
    plt.tick_params(axis='both',labelsize='xx-large')
    plt.grid(axis='both')
    plt.title('Cells in the bone marrow with bcr-abl and '+str(num_mutation)+' extra mutations',fontsize=20)
    plt.xlabel("Time (years)",fontsize=20)
    plt.ylabel("Number of cells",fontsize=20)
    plt.show()
#    plt.savefig(file_blood+name+'.png')
    
#     plt.figure(figsize=(16,8))
    plt.figure(2)
    plt.yscale('log')
#    plt.yscale('linear')
#    plt.xlim(2.5,4.0)
    plt.ylim(1e-1,1e11)
    for i in range(lev_initiate,lev_final):
        plt.scatter(time,data1[:,num_mutation*num_lev+i+1],marker=',',s=10.0,label='Level'+str(i%num_lev))
    plt.legend(bbox_to_anchor=(1.0,1.0), borderaxespad=0, fontsize=20)
    plt.tick_params(axis='both',labelsize='xx-large')
    plt.grid(axis='both')
    plt.title('Cells in the blood with bcr-abl and '+str(num_mutation)+' extra mutations',fontsize=20)
    plt.xlabel("Time (years)",fontsize=20)
    plt.ylabel("Number of cells",fontsize=20)
    plt.show()
#    plt.savefig(file_blood+name+'.png')
    
def plot_wbc():
#    plt.figure(1)
#     normal_value_low=[4e9 for i in range(np.size(data[:,0]))]
    normal_value_up=[1.0e+10 for i in range(long)]
#    zero_value=[0 for i in range(np.size(data[:,0]))]
#    chronic_phase=[10 for i in range(long)]
#    blast_crisis=[20 for i in range(long)]
    common_value_diagnosis=[1e11 for i in range(long)]
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.set_xlabel('time (years)')
    ax1.set_yscale('log')
#    ax1.set_ylim(1e9,1e14)
    ax1.set_ylabel('Number of WBC', color=color,fontsize=20)
    ax1.tick_params(axis='y', labelcolor=color,labelsize='xx-large')
    ax1.tick_params(axis='x',labelsize='xx-large')
#    ax1.fill_between(time,normal_value_low,normal_value_up,label='Normal region',alpha=0.3,color='salmon')
    ax1.plot(time,normal_value_up,label='Maximum healthy',color=color,linestyle=':')
    ax1.plot(time,common_value_diagnosis,label='Common value at diagnosis',color=color)
    ax1.scatter(time,total_wbc, color=color,label='WBC')
    ax1.legend(loc='best', borderaxespad=0, fontsize=20)
    
    color = 'tab:blue'
    ax2 = ax1.twinx()
    ticks=np.arange(0,101,10)
    ax2.set_yticks(ticks)
    plt.grid(axis='both')
    ax2.set_ylabel('Cell percentage', color=color,fontsize=20) 
    ax2.set_ylim(1e0,1e2)
#    ax2.set_yscale('log')
    ax2.tick_params(axis='y', labelcolor=color,labelsize='xx-large')
#    plt.fill_between(time,zero_value,chronic_phase,alpha=0.3,label='Chronic Phase')
#    plt.plot(time,chronic_phase,label='Chronic Phase',color=color,linestyle=':')
#    plt.fill_between(time,chronic_phase,blast_crisis,alpha=0.3,label='Acelerated phase')
#    plt.plot(time,blast_crisis,label='Acelerated phase',color=color,linestyle='-.')
#    plt.scatter(time,(cell_type_percentage[:,blast_level])*100,marker='D',s=1.0,label='Blasts', color=color)
#    col = ['violet','yellow','seagreen','orchid','cyan','goldenrod','darkseagreen','darkorange']
    ax2.scatter(time,(bcr_abl_percentage[:])*100,label='BCR-ABL percentage', color='blue')
    ax2.scatter(time,bcr_abl_cells[:,blast_level]*100,label='Percentage of blasts')
    ax2.tick_params(axis='x',labelsize='xx-large')
    ax2.legend(loc=9, borderaxespad=0, fontsize=20)
    fig.tight_layout()
#    plt.show()
    
def ratio_plots():
    plt.figure(7)
    global time
#    ticks=np.arange(0,101,5)
    for i in range(13,num_lev):
        plt.scatter(time,bcr_abl_cells[:,i]*100,label='k= '+str(i))
#    plt.yscale('log')
#    plt.ylim(1e-1,10e1)
    plt.grid('True',which='both')
    plt.legend(loc=9, borderaxespad=0, fontsize=20)
#    print(to_plot[long-1,:])
    plt.xlabel('time (years)',fontsize=20)
    plt.ylabel('Percentage of cells by level',fontsize=20)
    plt.title('Percentage of cells by level in blood',fontsize=20)
    plt.tick_params(axis='both',labelsize='xx-large')
#    plt.yticks(ticks)
    plt.show()
    
def ratio_plots_bone():
    plt.figure(8)
    global time
#    ticks=np.arange(0,56,5)
    for i in range(13,num_lev):
        plt.scatter(time,bcr_abl_cells_bone[:,i]*100,label=str(i))
#    plt.yscale('log')
    plt.ylim(1e-1,10e1)
    plt.grid('True',which='both')
    plt.legend(loc=9, borderaxespad=0, fontsize=20)
#    print(to_plot[long-1,:])
    plt.xlabel('time (years)',fontsize=20)
    plt.ylabel('Percentage of cells by level',fontsize=20)
    plt.title('Percentage of cells by level in bone marrow',fontsize=20)
    plt.tick_params(axis='both',labelsize='xx-large')
#    plt.yticks(ticks)
    plt.show()
    
def delta():
    plt.figure(5)
#    plt.ylim(1e-1,1e1)
#    plt.yscale('log')
#    plt.ylim(2000,3000)
    plt.scatter(data_delta[:,2],data_delta[:,1],label='leaking function')
    plt.xlabel('N(t)-N_0',fontsize=20)
    plt.ylabel('Leaking strength',fontsize=20)
    plt.title('Leaking function profile',fontsize=20)
#    plt.scatter(time,data[:,2],label='cells')
    plt.legend(loc='best', borderaxespad=0, fontsize=20)
    plt.tick_params(axis='both',labelsize='xx-large')
    plt.grid('True',which='both')
    plt.show()
    
def delta_time():
    plt.figure(6)
#    plt.ylim(1e-1,1e1)
#    plt.yscale('log')
#    plt.ylim(2000,3000)
    plt.scatter(time,data_delta[:,1],label='leaking function')
    plt.xlabel('time (years)',fontsize=20)
    plt.ylabel('Leaking strength',fontsize=20)
    plt.title('Leaking function profile',fontsize=20)
#    plt.scatter(time,data[:,2],label='cells')
    plt.legend(loc='best', borderaxespad=0, fontsize=20)
    plt.tick_params(axis='both',labelsize='xx-large')
    plt.grid('True',which='both')
    plt.show()
    
def delta_poster():
    plt.figure(7)
#    plt.ylim(1e-1,1e1)
#    plt.yscale('log')
#    plt.ylim(2000,3000)
#    for i in range(0,num_lev):
    plt.scatter(time,(data[:,0]+tot_percentage_bone[0])*data_delta[:,1],label='k= '+str(0))
    plt.scatter(time,(data[:,13]+tot_percentage_bone[13])*data_delta[:,1],label='k= '+str(13))
    plt.scatter(time,(data[:,19]+tot_percentage_bone[19])*data_delta[:,1],label='k= '+str(19))
    plt.xlabel('time (years)',fontsize=20)
    plt.ylabel('Number of cells',fontsize=20)
    plt.title('Migration of cells from bone marrow',fontsize=20)
#    plt.scatter(time,data[:,2],label='cells')
    plt.legend(loc='best', borderaxespad=0, fontsize=20)
    plt.tick_params(axis='both',labelsize='xx-large')
    plt.grid('True',which='both')
    plt.show()

delta()   
ratio_plots()
delta_time()
ratio_plots_bone()
graph_cells(0,num_lev,5)
#plot_wbc()
#delta_poster()