import numpy as np
import matplotlib.pyplot as plt
from One_D_Constants import One_D_Constants,gamma
from One_D_Functions import Lagran_multi, dPsidt, dPsidt_RungeKutta_4, Lagran_multi_V2
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
    df_name= args[15]
    fps_movie ,num_frames = args[16:18]

    multipliers= []
    sim_psi = np.zeros(shape=(2,N+1))
    frames_psi = np.zeros(shape=(num_frames,N+1))
    frame = 0
    print("Simulation progressbar")
    b = progressbar.ProgressBar(maxval=sim_steps-1)
    for time in range(0,sim_steps-1):
        b.update(time)

        t1, t2 = time%2 , (time + 1)%2

        x = Lagran_multi_V2(
                psi_list=sim_psi
                ,t=t1
                ,k=k,c0=c0,ds=ds
                ,linalg_lstsq=False
                ,num_chains=N
                        )


        for link in range(N-2,-1,-1):
            if link == N :
                sim_psi[t2][link] =  0
                print(link)
            if link < N:
                RungeKutta4 = dPsidt_RungeKutta_4(
                    link=link
                    ,N=N, ds=ds, dt=dt
                    ,multipliers=x
                    ,psi=sim_psi[t1][link]
                    )
                
                sim_psi[t2][link] = sim_psi[t1][link] + (dt/6)*RungeKutta4

        if time%int(sim_steps/num_frames )== 0:
            multipliers.append(x)
            for i in range(N+1):
                frames_psi[frame][i] = sim_psi[t1][i]
                
            frame += 1


    
    x_list,z_list = plot_from_psi(
        psi=frames_psi
        ,sim_steps=sim_steps
        ,ds=ds ,r0=r0, L=L
        )
  
    if save_data == True:
        df = pd.DataFrame({
            'Psi': [frames_psi],
            'x pos': [x_list],
            'z pos': [z_list],
            'multipliers': [multipliers],
            "L" : L,
            "r0": r0,
            "num of chain links": N,
            "Total time [sec]" : T,
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
    sim_1D_surface()
    Make_movie()