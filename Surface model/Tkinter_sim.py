import numpy as np
from Two_D_constants import Two_D_Constants, Two_D_paths
from Two_D_functions import Langrange_multi,dpsidt_func,drdt_func,dzdt_func, gamma
from Two_D_functions import Epsilon_values, c_diff_f,c_diff_g
from Two_D_functions import Epsilon_v2, c_diff, check_constraints_truth
from Two_D_functions import Make_variable_corrections
from two_d_data_processing import check_area, tot_area
import matplotlib.pyplot as plt
import os
import pandas as pd
import progressbar
import time


def Two_d_simulation_stationary_states(
    N:int,k:float,c0:float, dt:float, ds:float
    ,sigma:float,kG:float,tau:float
    ,sim_steps:int, L:float, r0:float
    ,radi:list,z_list:list
    ,Area:list,psi:list
    ,r_unperturb:list ,z_unperturb:list
    ,df_name:str,num_frames:str
    ,data_path:str
    ,save_data:bool = True 
    ,Tolerence = 1e-15
    #,do_correction = True
    ,area_testing = False
    ):
    np.set_printoptions(legacy='1.25')
    start_time = time.time()
    end_sim = False
    print_scale = (sim_steps-2)/1000
    Area_initial = np.sum(Area)
    
    lambs_save, nus_save = [], []
    correct_count_list = np.zeros(sim_steps-1)
    print("Simulation progressbar \n ")
    for t in range(sim_steps-1):
        if end_sim == True:
            break
        if int(t%print_scale) == 0 :
            print(f"completion : {round(t/(print_scale*10),1)}%      time since start = {round((time.time()-start_time)/60,3):.4f} min", end="\r")
        t1,t2 = t%2, (t+1)%2
        
        lambs,nus = Langrange_multi(
                N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                ,Area=Area
                ,psi=psi[t]
                ,radi=radi[t]
                ,z_list=z_list[t]
            )
        lambs_save.append(lambs)
        nus_save.append(nus)
        for i in range(N+1):
            if i == N:
                z_list[t+1][i] = z_list[t][i]
                radi[t+1][i] = radi[t][i]
                #psi[t+1][i] = psi[t][i]
            if i < N:
                z_list[t+1][i] = z_list[t][i] + dt*dzdt_func(i=i,Area=Area,radi=radi[t],nu=nus)

                drdt = drdt_func(
                            i=i
                            ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                            ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                            ,lamb=lambs,nu=nus
                            )
                radi[t+1][i] = radi[t][i] + dt*drdt

                dpsidt = dpsidt_func(
                                i=i
                                ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                ,lamb=lambs,nu=nus
                                )
                psi[t+1][i] = psi[t][i] + dt*dpsidt

            
        
        correction_count = Make_variable_corrections(
            N=N
            ,r=radi[t+1] ,z=z_list[t+1], psi=psi[t+1] 
            ,Area=Area ,Area_init=Area_initial
            ,Tolerence=Tolerence
            ,corr_max=20
        )
        correct_count_list[t] = correction_count
    
    print("\n")
    print(f"\n the simulation time={round((time.time()-start_time)/60,3)} min \n")
    plt.figure()
    plt.plot([i*dt for i in range(sim_steps-1)],correct_count_list,".-")
    plt.title("Number of corrections for each t")
    plt.xlabel("t [s]")

    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi],
            "r": [radi],
            "z": [z_list],
            "r unperturbed": [r_unperturb],
            "z unperturbed": [z_unperturb],
            "area list": [Area],
            'lambs': [lambs_save],
            'nus': [nus_save],
            "L" : L,
            "r0": r0,
            "N": N,
            "c0": c0,
            "k": k,
            "kG": kG,
            "sigma": sigma,
            "tau": tau,
            "sim_steps": sim_steps,
            "dt": dt,
            "ds": ds,
            "gam(i=0)": gamma(0),
            "gam(i>0)": gamma(5),
            "correction count": [correction_count],
            "tolerence":Tolerence,
            "sim completion":True
                        })


        if not os.path.exists(data_path):
            os.makedirs(data_path)
        df.to_pickle(data_path + df_name)

