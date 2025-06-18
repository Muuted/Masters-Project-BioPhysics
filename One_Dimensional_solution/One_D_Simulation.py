import numpy as np
import matplotlib.pyplot as plt
from One_D_Constants import One_D_Constants
from plotting_functions import plot_from_psi_V2 
from data_processing import cirle_fit
from One_D_sim_function import sim_1D_surface
import os
import pandas as pd
import progressbar
from Make_movie import Make_frames,Make_video



if __name__ == "__main__":
    args = One_D_Constants(print_val=True)

    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    df_name= args[15]

    
    sim_1D_surface(
        L=L,r0=r0,N=N,ds=ds,T=T,dt=dt
        ,psi_list=psi_list,k=k,c0=c0,sim_steps=sim_steps
        ,data_path=data_path ,df_name=df_name
    )
    

    plot_from_psi_V2(
        data_path=data_path
        ,df_name=df_name
    )
    
    Make_frames(
        data_path=data_path
        ,figs_save_path=video_fig_path
        ,df_name="1D surface membrane dynamics"
    )
    

    Make_video(
        output_path = video_save_path
        ,input_path = video_fig_path
        ,video_name = f"dynamics movie links={N}.avi"
        ,fps=8
    )

    x_cen,z_cen,Radius,x_circle ,z_circle = cirle_fit(
        data_path=data_path
        ,df_name=df_name
    )

    print(f"radius={Radius}")