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


def One_iteration_sim():
    args = One_D_Constants(print_val=True)

    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    df_name = args[15]

    
    sim_1D_surface(
        L=L,r0=r0,N=N,ds=ds,T=T,dt=dt
        ,psi_list=psi_list,k=k,c0=c0,sim_steps=sim_steps
        ,data_path=data_path 
        ,df_name=df_name
    )
    
    plot_from_psi_V2(
        data_path=data_path
        ,df_name=df_name
    )
    
    Make_frames(
        data_path=data_path
        ,figs_save_path=video_fig_path + df_name +"\\"
        ,df_name=df_name
    )
    
    Make_video(
        output_path = video_save_path
        ,input_path = video_fig_path + df_name +"\\"
        ,video_name = f"dynamics movie links={N}.avi"
        ,fps=8
    )

    x_cen,z_cen,Radius = cirle_fit(
        data_path=data_path
        ,df_name=df_name
    )

    print(f" \n radius={Radius}")


def Multiple_iteration_sim():
    args = One_D_Constants(print_val=True)

    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    df_name = args[15]

    c0_list = [ c0/2 , c0 , c0*2]

    for i in range(len(c0_list)):
        c0 = c0_list[i]
        df_name_1 = df_name + f" c0={c0}  sim time={T}s"
        
        sim_1D_surface(
            L=L,r0=r0,N=N,ds=ds,T=T,dt=dt
            ,psi_list=psi_list,k=k,c0=c0,sim_steps=sim_steps
            ,data_path=data_path 
            ,df_name=df_name_1
        )
        
        plot_from_psi_V2(
            data_path=data_path
            ,df_name=df_name_1
        )
        
        Make_frames(
            data_path=data_path
            ,figs_save_path=video_fig_path + df_name_1 +"\\"
            ,df_name= df_name_1
        )
        
        Make_video(
            output_path = video_save_path
            ,input_path = video_fig_path + df_name_1 +"\\"
            ,video_name = f"dynamics movie links={N}, c0={c0}.avi"
            ,fps=8
        )

        x_cen,z_cen,Radius = cirle_fit(
            data_path=data_path
            ,df_name=df_name_1
        )

        print(f" \n radius={Radius}")


def Rolling_test():
    args = One_D_Constants(print_val=True)

    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    df_name = args[15]

    #sim_steps = int(sim_steps/100)
    dc0 = c0
    c0 = 1
    iter_count = 0
    Rolling = False
    while Rolling == False:
        print(f"\n \n iter counter = {iter_count} and c0={c0} \n \n")

        sim_1D_surface(
            L=L,r0=r0,N=N,ds=ds,T=T,dt=dt
            ,psi_list=psi_list,k=k,c0=c0,sim_steps=sim_steps
            ,data_path=data_path 
            ,df_name=df_name
        )
        
        plot_from_psi_V2(
            data_path=data_path
            ,df_name=df_name
        )

        df_sim = pd.read_pickle(data_path + df_name)
        x = df_sim["x pos"][0]


        for t in range(sim_steps):
            if x[t][0] > x[t][1]:
                Rolling = True

        if Rolling == True:
            df_sim["c0"] = c0
            df_sim.to_pickle(data_path + df_name)

            Make_frames(
                data_path=data_path
                ,figs_save_path=video_fig_path + df_name +"\\"
                ,df_name=df_name
            )
            
            Make_video(
                output_path = video_save_path
                ,input_path = video_fig_path + df_name +"\\"
                ,video_name = f"dynamics movie links={N}.avi"
                ,fps=8
            )

            x_cen,z_cen,Radius = cirle_fit(
                data_path=data_path
                ,df_name=df_name
            )

            print(f"Rolling was found at c0={c0}")
            break

        if Rolling == False:
            c0 += dc0
            iter_count += 1

if __name__ == "__main__":
    One_iteration_sim()
    #Multiple_iteration_sim()
    #Rolling_test()
    print("dont run unless you want new data")
    pass