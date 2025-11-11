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

def Two_D_simulation(
    N:int,k:float,c0:float, dt:float, ds:float
    ,sigma:float,kG:float,tau:float
    ,sim_steps:int, L:float, r0:float
    ,radi:list,z_list:list
    ,Area:list,psi:list
    ,df_name:str,num_frames:str
    ,data_path:str
    ,save_data:bool = True 
    ,condition = 1
    ):

    
    lambs_save = []
    nus_save = []

    print("Simulation progressbar \n ")
    b = progressbar.ProgressBar(maxval=sim_steps-1)
    b.start()
    for t in range(sim_steps-1):
        b.update(t)
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
                psi[t+1][i] = psi[t][i]
            if i < N:
                z_list[t+1][i] = z_list[t][i] + dt*dzdt_func(i=i,Area=Area,radi=radi[t],nu=nus)

                radi[t+1][i] = radi[t][i] + dt*drdt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )

                psi[t+1][i] = psi[t][i] + dt*dpsidt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )
                
            if t%condition==0 and i < N:
                ebf, ebg = Epsilon_values(
                    N=N, r=radi[t+1], z=z_list[t+1] ,psi=psi[t+1] ,Area=Area
                            )
                
                K_r,K_z,K_psi = 0,0,0
                for beta in range(N):
                    ebf_val, ebg_val = ebf[beta] ,ebg[beta]
                    
                    K_r += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="r") 
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="r")
                        )
                    
                    K_z += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="z")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="z")
                        )
                    
                    K_psi += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="psi")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="psi")
                        )
                    
                radi[t+1][i] += K_r
                z_list[t+1][i] += K_z
                psi[t+1][i] += K_psi

    print("\n")
    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi],
            "r":[radi],
            "z":[z_list],
            "area list":[Area],
            'lambs': [lambs_save],
            'nus': [nus_save],
            "L" : L,
            "r0": r0,
            "N": N,
            "c0": c0,
            "k":k,
            "kG":kG,
            "sigma":sigma,
            "tau":tau,
            "sim_steps": sim_steps,
            "dt":dt,
            "ds":ds,
            "gam(i=0)":gamma(0),
            "gam(i>0)":gamma(5)
                        })

        #print(df.info())W


        if not os.path.exists(data_path):
            os.makedirs(data_path)
        df.to_pickle(data_path + df_name)




def Two_D_simulation_V2(
    N:int,k:float,c0:float, dt:float, ds:float
    ,sigma:float,kG:float,tau:float
    ,sim_steps:int, L:float, r0:float
    ,radi:list,z_list:list
    ,Area:list,psi:list
    ,df_name:str,num_frames:str
    ,data_path:str
    ,save_data:bool = True 
    ,condition = 1
    ,Tolerence = 1e-15
    ):

    Area_old = np.sum(Area)
    Area_new = 0
    correction_count = 0
    pos_count = 0
    lambs_save, nus_save = [], []

    print("Simulation progressbar \n ")
    b = progressbar.ProgressBar(maxval=sim_steps-1)
    b.start()
    for t in range(sim_steps-1):
        b.update(t)
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
                psi[t+1][i] = psi[t][i]
            if i < N:
                z_list[t+1][i] = z_list[t][i] + dt*dzdt_func(i=i,Area=Area,radi=radi[t],nu=nus)

                radi[t+1][i] = radi[t][i] + dt*drdt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )

                psi[t+1][i] = psi[t][i] + dt*dpsidt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )
            
        Area_new = tot_area(N=N,r=radi[t+1],z=z_list[t+1])
        dA = np.abs(Area_new - Area_old)
        pos_count += 1
        if Tolerence < dA:
            correction_count += 1
            for i in range(N):
                ebf, ebg = Epsilon_values(
                    N=N, r=radi[t+1], z=z_list[t+1] ,psi=psi[t+1] ,Area=Area
                            )
                K_r,K_z,K_psi = 0,0,0
                for beta in range(N):
                    ebf_val, ebg_val = ebf[beta] ,ebg[beta]
                    
                    K_r += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="r") 
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="r")
                        )
                    
                    K_z += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="z")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="z")
                        )
                    
                    K_psi += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="psi")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="psi")
                        )
                    
                radi[t+1][i] += K_r
                z_list[t+1][i] += K_z
                psi[t+1][i] += K_psi

    print("\n")
    #print(f"correction count={correction_count} of {sim_steps*(N-1)} possible")
    print(f"correction count={correction_count} of {pos_count} possible")

    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi],
            "r": [radi],
            "z": [z_list],
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
            "correction count": correction_count
                        })


        if not os.path.exists(data_path):
            os.makedirs(data_path)
        df.to_pickle(data_path + df_name)




def Two_D_simulation_V3(
    N:int,k:float,c0:float, dt:float, ds:float
    ,sigma:float,kG:float,tau:float
    ,sim_steps:int, L:float, r0:float
    ,radi:list,z_list:list
    ,Area:list,psi:list
    ,df_name:str,num_frames:str
    ,data_path:str
    ,save_data:bool = True 
    ,condition = 1
    ,Tolerence = 1e-15
    ):

    Area_old = np.sum(Area)
    Area_new = 0
    correction_count = 0
    pos_count = 0
    lambs_save, nus_save = [], []

    print("Simulation progressbar \n ")
    b = progressbar.ProgressBar(maxval=sim_steps-1)
    b.start()
    for t in range(sim_steps-1):
        b.update(t)
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
                psi[t+1][i] = psi[t][i]
            if i < N:
                z_list[t+1][i] = z_list[t][i] + dt*dzdt_func(i=i,Area=Area,radi=radi[t],nu=nus)

                radi[t+1][i] = radi[t][i] + dt*drdt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )

                psi[t+1][i] = psi[t][i] + dt*dpsidt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )
            
        Area_new = tot_area(N=N,r=radi[t+1],z=z_list[t+1])
        dA = np.abs(Area_new - Area_old)
        do_correction = False
        if Tolerence < dA:
            do_correction = True
        pos_count += 1
        while do_correction == True:
            correction_count += 1
            for i in range(N):
                ebf, ebg = Epsilon_values(
                    N=N, r=radi[t+1], z=z_list[t+1] ,psi=psi[t+1] ,Area=Area
                            )
                K_r,K_z,K_psi = 0,0,0
                for beta in range(N):
                    ebf_val, ebg_val = ebf[beta] ,ebg[beta]
                    
                    K_r += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="r") 
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="r")
                        )
                    
                    K_z += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="z")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="z")
                        )
                    
                    K_psi += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="psi")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="psi")
                        )
                    
                radi[t+1][i] += K_r
                z_list[t+1][i] += K_z
                psi[t+1][i] += K_psi

            Area_new = tot_area(N=N,r=radi[t+1],z=z_list[t+1])
            dA = np.abs(Area_new - Area_old)
            if Tolerence > dA :
                do_correction = False
                break
    
    print("\n")
    #print(f"correction count={correction_count} of {sim_steps*(N-1)} possible")
    print(f"correction count={correction_count} of {pos_count} possible")

    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi],
            "r": [radi],
            "z": [z_list],
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
            "correction count": correction_count
                        })


        if not os.path.exists(data_path):
            os.makedirs(data_path)
        df.to_pickle(data_path + df_name)






def Two_d_simulation_stationary_states(
    N:int,k:float,c0:float, dt:float, ds:float
    ,sigma:float,kG:float,tau:float, eta:float
    ,sim_steps:int, L:float, r0:float
    ,radi:list,z_list:list
    ,Area:list,psi:list
    ,r_unperturb:list ,z_unperturb:list
    ,df_name:str,num_frames:str
    ,data_path:str
    ,save_data:bool = True 
    ,Tolerence = 1e-10
    ,do_correction = True
    ,area_testing = False
    ):
    np.set_printoptions(legacy='1.25')
    start_time = time.time()
    print_scale = (sim_steps-2)/1000
    Area_initial = np.sum(Area)
    
    #lambs_save, nus_save = [], []
    correct_count_list = np.zeros(sim_steps-1)
    print("Simulation progressbar \n ")
    for t in range(sim_steps-1):
        if int(t%print_scale) == 0 :
            time_since_start = round((time.time()-start_time)/60,3)
            print(f"completion : {round(t/(print_scale*10),1)}%"
                  +f"   time since start = {time_since_start:.4f} min"
                  +f"   estimated time left = {(time_since_start/(t+1))*(sim_steps-t):.2f} min"
                  +f"       {((time_since_start/(t+1))*(sim_steps-t))/60:.2f} hours"
                  , end="\r"
                  )
        #t1,t2 = t%2, (t+1)%2
        
        lambs,nus = Langrange_multi(
                N=N,k=k,c0=c0,sigma=sigma
                ,kG=kG,tau=tau,ds=ds,eta=eta
                ,Area=Area
                ,psi=psi[t]
                ,radi=radi[t]
                ,z_list=z_list[t]
            )
        #lambs_save.append(lambs)
        #nus_save.append(nus)

        for i in range(N+1):
            if i == N:
                z_list[t+1][i] = z_list[t][i]
                radi[t+1][i] = radi[t][i]
                #psi[t+1][i] = psi[t][i]
            if i < N:
                z_list[t+1][i] = z_list[t][i] + dt*dzdt_func(i=i,ds=ds,eta=eta,Area=Area,radi=radi[t],nu=nus)

                drdt = drdt_func(
                            i=i
                            ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,ds=ds,eta=eta
                            ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                            ,lamb=lambs,nu=nus
                            )
                radi[t+1][i] = radi[t][i] + dt*drdt

                dpsidt = dpsidt_func(
                                i=i
                                ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,ds=ds,eta=eta
                                ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                ,lamb=lambs,nu=nus
                                )
                psi[t+1][i] = psi[t][i] + dt*dpsidt
        
        if do_correction == True:
            correction_count = Make_variable_corrections(
                N=N
                ,r=radi[t+1] ,z=z_list[t+1], psi=psi[t+1] 
                ,Area=Area ,Area_init=Area_initial
                ,Tolerence=Tolerence
                ,corr_max=100
                ,t=t
            )
            correct_count_list[t] = correction_count
    
    #b.finish()
    print("\n")
    print(f"\n the simulation time={round((time.time()-start_time)/60,3)} min \n"
          +f"round((time.time()-start_time)/60**2,3) hours"
          )
    

    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi],
            "r": [radi],
            "z": [z_list],
            "r unperturbed": [r_unperturb],
            "z unperturbed": [z_unperturb],
            "area list": [Area],
            #'lambs': [lambs_save],
            #'nus': [nus_save],
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
            "gam(i=0)": gamma(0,ds=ds,eta=eta),
            "gam(i>0)": gamma(2,ds=ds,eta=eta),
            "correction count": [correct_count_list],
            "tolerence":Tolerence,
            "sim completion":True
                        })


        if not os.path.exists(data_path):
            os.makedirs(data_path)
        df.to_pickle(data_path + df_name)

if __name__ == "__main__":
    pass
    exit()