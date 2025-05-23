import numpy as np
import matplotlib.pyplot as plt
from One_D_Constants import One_D_Constants
from One_D_Functions import Lagran_multi, dPsidt
from plotting_functions import plot_from_psi
import os
import pandas as pd
import progressbar
from Make_movie import Make_movie

def sim_1D_surface(
        save_data=True
                  ):

    args = One_D_Constants(print_val=True)

    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]

    sim_T_tot = int(T/dt)
    multipliers= []

    b = progressbar.ProgressBar(maxval=sim_steps)
    for time in range(sim_steps-1):
        b.update(time)
        multipliers.append(Lagran_multi(
            psi_list=psi_list
            ,t=time,k=k,c0=c0,ds=ds
            ,linalg_lstsq=False
        )
        )

        for link in range(N-2,-1,-1):
            dPsidt(
                i=link, t=time, dt=dt, deltaS=ds
                ,multi=multipliers[time]
                ,psi=psi_list
            )

    x_list = []
    z_list = []
    b2 = progressbar.ProgressBar(maxval=sim_steps)
    for t in range(sim_steps-1):
        b2.update(t + 1)
        x0,z0 = plot_from_psi(psi=psi_list[t],ds=ds,r0=r0)
        x_list.append(x0)
        z_list.append(z0)

    
    if save_data == True:
        df_name= "\\1D surface membrane dynamics"
        df = pd.DataFrame({
            'Psi': [psi_list],
            'x pos': [x_list],
            'z pos': [z_list],
            'multipliers': [multipliers],
            "L" : L,
            "r0": r0,
            "num of chain links": N,
            "Total time [sec]" : T,
            "dt":dt
                        })

        print(df.info())


        if not os.path.exists(data_path):
            os.makedirs(data_path)

        
        df.to_pickle(data_path + df_name)



if __name__ == "__main__":
    sim_1D_surface()
    Make_movie()