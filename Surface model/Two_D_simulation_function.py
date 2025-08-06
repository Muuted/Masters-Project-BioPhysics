import numpy as np
from Two_D_constants import Two_D_Constants, gamma, Two_D_paths
from Two_D_functions import Langrange_multi,dpsidt_func,drdt_func,dzdt_func
from Two_D_functions import Epsilon_values, c_diff_f,c_diff_g
from two_d_data_processing import check_area, tot_area
import matplotlib.pyplot as plt
import os
import pandas as pd
import progressbar


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
                    N=N, r=radi[t], z=z_list[t] ,psi=psi[t] ,Area=Area
                            )
                K_r,K_z,K_psi = 0,0,0
                for beta in range(N):
                    ebf_val, ebg_val = ebf[beta] ,ebg[beta]
                    
                    K_r += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t],psi=psi[t],Area=Area,diff_var="r") 
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t],psi=psi[t],z=z_list[t],Area=Area,diff_var="r")
                        )
                    
                    K_z += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t],psi=psi[t],Area=Area,diff_var="z")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t],psi=psi[t],z=z_list[t],Area=Area,diff_var="z")
                        )
                    
                    K_psi += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t],psi=psi[t],Area=Area,diff_var="psi")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t],psi=psi[t],z=z_list[t],Area=Area,diff_var="psi")
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
            
            Area_new = tot_area(N=N,r=radi[t],z=z_list[t])
            dA = np.abs(Area_new - Area_old)
            pos_count += 1
            if i < N and Tolerence < dA :
                correction_count += 1 
                ebf, ebg = Epsilon_values(
                    N=N, r=radi[t], z=z_list[t] ,psi=psi[t] ,Area=Area
                            )
                K_r,K_z,K_psi = 0,0,0
                for beta in range(N):
                    ebf_val, ebg_val = ebf[beta] ,ebg[beta]
                    
                    K_r += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t],psi=psi[t],Area=Area,diff_var="r") 
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t],psi=psi[t],z=z_list[t],Area=Area,diff_var="r")
                        )
                    
                    K_z += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t],psi=psi[t],Area=Area,diff_var="z")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t],psi=psi[t],z=z_list[t],Area=Area,diff_var="z")
                        )
                    
                    K_psi += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t],psi=psi[t],Area=Area,diff_var="psi")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t],psi=psi[t],z=z_list[t],Area=Area,diff_var="psi")
                        )
                    
                radi[t+1][i] += K_r
                z_list[t+1][i] += K_z
                psi[t+1][i] += K_psi

            Area_old = Area_new
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


if __name__ == "__main__":
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    #sim_steps = 1000

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    Two_D_simulation(
        N=N,k=k,c0=c0,sigma=sigma, dt=dt
        ,kG=kG,tau=tau, sim_steps=sim_steps
        ,L=L, r0=r0
        ,Area=Area_list
        ,psi=psi_list
        ,radi=radi_list
        ,z_list=z_list
        ,df_name = df_name
        ,num_frames = num_frames
        ,data_path = data_path
    )

    df_sim = pd.read_pickle(data_path + df_name)
    r = df_sim['r'][0]        
    Area = df_sim['area list'][0]

    error = False
    for t in range(sim_steps):
        error = check_area(
            t=t,N=N,r=r[t],Area=Area
        )
        if error == True:
            print(
                f"\n ----------------------------- \n"
                +f"error at dt={dt}"
                +f"\n ----------------------------- \n"
                )
            break
    i = 0
    Area_change = []
    for t in range(sim_steps):
        dA = np.pi*( r[t][i+1]**2 - r[t][i]**2 ) - Area[i] 
        Area_change.append(dA)

    #print(Area_change)
    print(len(Area_change))
    plt.figure()
    plt.plot(Area_change[0:len(Area_change)-1],'-o')
    plt.show()
    