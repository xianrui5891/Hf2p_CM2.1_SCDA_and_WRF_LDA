"""
主程序位置
"""
import os
import sys

from cm2 import *
import numpy as np

import mpi4py
mpi4py.rc.initialize = False  
mpi4py.rc.finalize   = False 
from mpi4py import MPI


#from latent import latent_cda


#declare the variables of the bottom of atmosphere
U_bot = np.zeros((144,90),dtype=np.float64, order='F')
V_bot = np.zeros((144,90),dtype=np.float64, order='F')
T_bot = np.zeros((144,90),dtype=np.float64, order='F') 
P_bot = np.zeros((144,90),dtype=np.float64, order='F') 

#declare the variables of the surface of sea
SSU = np.zeros((360,200),dtype=np.float64, order='F')
SSV = np.zeros((360,200),dtype=np.float64, order='F')
SST = np.zeros((360,200),dtype=np.float64, order='F')
SSZ = np.zeros((360,200),dtype=np.float64, order='F')

#define the parameters of the model
num_cpld_calls = None
num_atms_calls = None
no_latent_scda = 100


#init the model
#print(cm2_cda_plugs.cm2_cda_maininit.__doc__)
num_cpld_calls, num_atms_calls = cm2_cda_plugs.cm2_cda_maininit()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
da_intervals = 12
da_aim_tag = 6


#operators
is_figure = False
is_save_nc = True
is_save_all_in_one = True
is_scda = True
sv = None
cda_model = None
plotter = None

if rank==0: 
    if is_scda:
        import torch
        from vae_cda import *
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        cda_model = latent_space_da(da_interval=da_intervals, dev=device)
    
    from data_output_utils import FieldPlotter
    plotter = FieldPlotter("./python_plots")

def save_data(time_str, p_bot, sst, data_type):
    plotter.plot_and_save(time_str=time_str, var_name='p_bot', data2d=p_bot, data_type=data_type)
    plotter.plot_and_save(time_str=time_str, var_name='sst', data2d=sst, data_type=data_type)



nu_scda = 0

print("num_cpld_calls,num_atms_calls in main",num_cpld_calls,num_atms_calls)

#start the integral process
for nc in range(1,num_cpld_calls+1):
    cm2_cda_plugs.atmos_step_pre(nc)

    for na in range(1,num_atms_calls + 1):
        if (nc%da_intervals == (da_aim_tag + 1)%da_intervals) and (na==1) and (nu_scda==1):
            atm_nu_scda = 1
        else:
            atm_nu_scda = 0

        cm2_cda_plugs.atmos_step(U_bot, V_bot, T_bot, P_bot, nc, na, atm_nu_scda)
    
        if rank==0: 
            print("PE_id=",rank,"for nc,na,nu,U_bot[72,:]=", nc, na, atm_nu_scda, U_bot[72,:])

        if rank==0 and is_figure: 
            print("PE_id=",rank,"for nc,na,nu,U_bot[72,:]=", nc, na, atm_nu_scda, U_bot[72,:])
            sv.plot_vector_field(U_bot.copy(), V_bot.copy(), T_bot.copy(), "atm", rank, nc, na)

    cm2_cda_plugs.ocean_step_pre()

    cm2_cda_plugs.ocean_step(SSU, SSV, SST, SSZ, nc, na, nu_scda)

    nu_scda = 0

    if( nc%da_intervals == da_aim_tag%da_intervals ):
        nu_scda = 1
        if rank==0 and is_scda:
            tempp_bot = P_bot
            tempsst = SST
            save_data(time_str=f"{nc}", p_bot=P_bot, sst=SST, data_type="original")

            atm_after, ocn_after = cda_model.do_da(nc, [T_bot, U_bot, V_bot, P_bot], [SST, SSU, SSV, SSZ])
            T_bot, U_bot, V_bot, P_bot = atm_after
            SST, SSU, SSV, SSZ = ocn_after
            save_data(time_str=f"{nc}", p_bot=P_bot, sst=SST, data_type="da")
            save_data(time_str=f"{nc}", p_bot=P_bot-tempp_bot, sst=SST-tempsst, data_type="increment")

        comm.Bcast(T_bot),comm.Bcast(U_bot),comm.Bcast(V_bot),comm.Bcast(P_bot)
        comm.Bcast(SST),comm.Bcast(SSU),comm.Bcast(SSV),comm.Bcast(SSZ)
    
    cm2_cda_plugs.restart_step(nc)
#end of the process

cm2_cda_plugs.main_end()
