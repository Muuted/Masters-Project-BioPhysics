import numpy as np
from Two_D_constants import Two_D_Constants, gamma, Two_D_paths
from Two_D_functions import Langrange_multi,dpsidt_func,drdt_func,dzdt_func
import matplotlib.pyplot as plt
import os
import pandas as pd
import progressbar


def Two_D_simulation(
    N:int,k:float,c0:float
    ,sigma:float,kG:float,tau:float
    ,sim_steps:int
    ,radi:list,z_list:list
    ,Area:list,psi:list
    ,df_name:str,num_frames:str
    ,data_path:str
    ,save_data:bool = True 
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
        ,Area=Area_list
        ,psi=psi_list[t]
        ,radi=radi_list[t]
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
                

    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi_list],
            "r":[radi],
            "z":[z_list],
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
            "Total time [sec]" : sim_steps*dt,
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

    plt.subplots()        
    plt.plot(radi[sim_steps-1],z_list[sim_steps-1])

    plt.show()

if __name__ == "__main__":
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    #sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    Two_D_simulation(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau, sim_steps=sim_steps
        ,Area=Area_list
        ,psi=psi_list
        ,radi=radi_list
        ,z_list=z_list
        ,df_name = df_name
        ,num_frames = num_frames
        ,data_path = data_path
    )