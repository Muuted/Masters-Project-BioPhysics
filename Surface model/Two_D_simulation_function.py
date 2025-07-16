import numpy as np
from Two_D_constants import Two_D_Constants, gamma
#import os
#import pandas as pd
#import progressbar



def Two_D_simulation(
    N,k,c0,sigma,kG,tau
    ,radi,z_list
    ,Area,psi
    ):
    pass


if __name__ == "__main__":
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]



    Two_D_simulation(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area_list
        ,psi=psi_list[0]
        ,radi=radi_list
        ,z_list=z_list
    )