from Two_D_constants import Two_D_Constants, Two_D_paths
from Two_D_simulation_function import Two_D_simulation
from Make_movie import Make_frames, Make_video
from two_d_data_processing import check_area
from Two_D_functions import Langrange_multi
import pandas as pd

if __name__ == "__main__":
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    #sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    make_movies = False
    run = True
    while run==True:
        Two_D_simulation(
            N=N,k=k,c0=c0,sigma=sigma, dt=dt
            ,kG=kG,tau=tau,L=L,r0=r0
            , sim_steps=sim_steps 
            ,Area=Area
            ,psi=psi_list
            ,radi=radi_list
            ,z_list=z_list
            ,df_name = df_name
            ,num_frames = num_frames
            ,data_path = data_path
        )
        if make_movies == True:
            Make_frames(
            data_path=data_path
            ,figs_save_path=video_fig_path
            ,df_name=df_name
            )

            Make_video(
                output_path = video_save_path
                ,input_path = video_fig_path
                ,video_name = df_name
                ,fps=12
            )
        
        df_sim = pd.read_pickle(data_path + df_name)
        r = df_sim['r'][0]        
        Area = df_sim['area list'][0]

        error = False
        for t in range(sim_steps):
            error = check_area(
                t=t,N=N,r=r,Area=Area
            )
            if error == True:
                print(
                    f"\n ----------------------------- \n"
                    +f"error at dt={dt}"
                    +f"\n ----------------------------- \n"
                    )
                break
        if error == True:
            break
        if error == False:
            dt *= 10
        
        