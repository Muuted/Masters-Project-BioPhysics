import numpy as np
import matplotlib.pyplot as plt
from One_D_Constants import One_D_Constants,gamma
from One_D_Functions import Lagran_multi, dPsidt, dPsidt_RungeKutta_4, Lagran_multi_V2
from plotting_functions import plot_from_psi_V2 
import os
import pandas as pd
import progressbar
from Make_movie import Make_frames,Make_video

def sim_1D_surface(
        L,r0,N,ds,T,dt
        ,psi_list,k,c0,sim_steps
        ,data_path ,df_name
        ,save_data=True
                  ):

    multipliers= []

    print("Simulation progressbar")
    b = progressbar.ProgressBar(maxval=sim_steps-1)
    for time in range(0,sim_steps-1):
        b.update(time)
        #t1, t2 = time%2 , (time + 1)%2
        x = Lagran_multi_V2(
                psi_list=psi_list
                ,t=time,k=k,c0=c0,ds=ds
                ,linalg_lstsq=False
                ,num_chains=N
                        )
        multipliers.append(x)

        for link in range(len(psi_list[0])-2,-1,-1):
            if link == N :
                psi_list[time+1][link] =  0
                
            if link < N:
                RungeKutta4 = dPsidt_RungeKutta_4(
                    link=link
                    ,N=N, ds=ds, dt=dt
                    ,multipliers=x#multipliers[time]
                    ,psi=psi_list[time][link]
                    )
                
                psi_list[time+1][link] = psi_list[time][link] + (dt/6)*RungeKutta4
                
  
    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi_list],
            #'x pos': [x_list],
            #'z pos': [z_list],
            'multipliers': [multipliers],
            "L" : L,
            "r0": r0,
            "N": N,
            "Total time [sec]" : T,
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