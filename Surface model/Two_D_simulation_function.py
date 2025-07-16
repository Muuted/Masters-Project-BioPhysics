import numpy as np
from Two_D_constants import Two_D_Constants, gamma
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
    ,save_data = True
    ):
    
    for t in range(sim_steps):
        lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area_list
        ,psi=psi_list[t]
        ,radi=radi_list[t]
        ,z_list=z_list[t]
            )
        for i in range(N):
            
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

    sim_steps = 10

    
    Two_D_simulation(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau, sim_steps=sim_steps
        ,Area=Area_list
        ,psi=psi_list
        ,radi=radi_list
        ,z_list=z_list
    )