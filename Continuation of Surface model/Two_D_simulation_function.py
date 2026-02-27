import numpy as np
from Two_D_constants import Two_D_Constants, Two_D_paths
from Two_D_functions import Langrange_multi,dpsidt_func,drdt_func,dzdt_func, gamma
from Two_D_functions import Epsilon_values, c_diff_f,c_diff_g
from Two_D_functions import Epsilon_v2, c_diff, check_constraints_truth
from Two_D_functions import Make_variable_corrections
from two_d_data_processing import check_area, tot_area, E_kin,E_pot, Xsqaured_test
from Runge_Kutta import RungeKutta45

import matplotlib.pyplot as plt
import os
import pandas as pd
import progressbar
import time


def Two_d_simulation_stationary_states(
    N:int,k:float,c0:float, dt:float, ds:float
    ,sigma:float,kG:float,tau:float, eta:float
    ,sim_steps:int, L:float, r0:float, dpsi_perturb:float
    ,radi:list,z_list:list
    ,Area:list,psi:list
    ,r_unperturb:list ,z_unperturb:list, psi_unperturb:list
    ,df_name:str,num_frames:str
    ,data_path:str
    ,integration_method:str #= "RK4"
    ,save_data:bool = True 
    ,Tolerence:float = 1e-10
    ,do_correction:bool = True
    ,area_testing:bool = False
    ,print_progress:bool = True
    ):
    np.set_printoptions(legacy='1.25')
    start_time = time.time()
    print_scale = (sim_steps-2)/1000
    Area_initial = np.sum(Area)
    
    Potential_E = np.zeros(sim_steps-1,dtype=float)
    Potential_E_before_correction = np.zeros(sim_steps-1,dtype=float)
    Kinetic_E = np.zeros(sim_steps-1,dtype=float)
    Xsqre = np.zeros(sim_steps-1,dtype=float)
    correct_count_list = np.zeros(sim_steps-1)

    integration_options = ["Euler","RK4"]
    if integration_method not in integration_options:
        print("No integration method choosen correctly")
        exit()
    print("integration_method=",integration_method)
    print("Simulation progressbar \n ")
    for t in range(sim_steps-1):
        if int(t%print_scale) == 0 and print_progress == True:
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
        if integration_method == "Euler":
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

        if integration_method == "RK4":
            kr,kz,kpsi = RungeKutta45(
                N=N,dt=dt,k=k,c0=c0, sigma=sigma
                ,kG=kG ,tau=tau, ds=ds,eta=eta
                ,Area=Area,lamb=lambs,nu=nus
                ,psi_init=psi[t],r_init=radi[t], z_init=z_list[t]
                )
            for i in range(N+1):
                if i == N:
                    z_list[t+1][i] = z_list[t][i]
                    radi[t+1][i] = radi[t][i]
                    #psi[t+1][i] = psi[t][i]
                if i < N:
                    radi[t+1][i] = radi[t][i] + (dt/6)*(kr[1][i] + kr[2][i] + kr[3][i] + kr[4][i])
                    z_list[t+1][i] = z_list[t][i] + (dt/6)*(kz[1][i] + kz[2][i] + kz[3][i] + kz[4][i])
                    psi[t+1][i] = psi[t][i] + (dt/6)*(kpsi[1][i] + kpsi[2][i] + kpsi[3][i] + kpsi[4][i])


        Potential_E_before_correction[t] = E_pot(N=N,k=k,kG=kG,tau=tau,c0=c0,r=radi[t],psi=psi[t],Area=Area)

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
    

        Potential_E[t] = E_pot(N=N,k=k,kG=kG,tau=tau,c0=c0,r=radi[t],psi=psi[t],Area=Area)
        Kinetic_E[t] = E_kin(N=N,t=t,dt=dt,r=radi,z=z_list,Area=Area)
        Xsqre[t] = Xsqaured_test(N=N
                                 ,r_init=r_unperturb,z_init=z_unperturb,psi_init=psi_unperturb
                                 ,r=radi[t],z=z_list[t],psi=psi[t])

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
            "psi unperturbed": [psi_unperturb],
            "area list": [Area],
            "Epot":[Potential_E],
            "Ekin":[Kinetic_E],
            "Ekin before correction": [Potential_E_before_correction],
            "Chi squared test": [Xsqre],
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
            "dpsi perturb": dpsi_perturb,
            "gam(i=0)": gamma(0,ds=ds,eta=eta),
            "gam(i>0)": gamma(2,ds=ds,eta=eta),
            "correction count": [correct_count_list],
            "tolerence":Tolerence,
            "sim completion":True,
            "integration method":integration_method,
            "simulation time [s]":int(time.time()-start_time)
                        })

        print(data_path + df_name)
        if not os.path.exists(data_path):
            os.makedirs(data_path)
        df.to_pickle(data_path + df_name )#+".pkl")
        

if __name__ == "__main__":
    pass
    exit()