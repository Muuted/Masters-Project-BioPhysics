import pandas as pd
import time
import matplotlib.pyplot as plt
from Two_D_constants import Two_D_Constants, Two_D_paths, Two_D_Constants_stationary_state
from Two_D_simulation_function import Two_D_simulation, Two_D_simulation_V2, Two_d_simulation_stationary_states
from Make_movie import Make_frames, Make_video
from two_d_data_processing import check_area
from Two_D_functions import Langrange_multi
from two_d_plot_data import plot_Epot_Ekin, plot_tot_area


def Surface_sim_stationary_state_initial_configuration_iterative(
        do_simulation:bool = True
        ,start_from_flat:bool = False
        ,do_perturbation:bool = False
        ,make_movie: bool = True
        ,make_plots: bool = True
    ):
    print("\n Now Running the surface simulation from stationary configurations \n")


    folder_name = "Science stability test\\"
    save_name_list = [" unperturbed", " -perturbed"," +perturbed"]
    dpsi_list = [0, -0.01, 0.01]
    do_perturbation_list = [False,True,True]

    ref_list_len = len(save_name_list)
    for i in range(ref_list_len):
        const_args = Two_D_Constants_stationary_state(
            print_val=True
            ,show_stationary_state=True
            ,start_flat=start_from_flat
            ,perturb=do_perturbation_list[i]
            ,dpsi_perturb= dpsi_list[i]
        )

        L,r0,N,ds,T,dt = const_args[0:6]
        k,c0,sim_steps = const_args[6:9]
        sigma, tau, kG = const_args[9:12]
        Area_list, psi_list = const_args[12:14]
        radi_list,z_list = const_args[14:16]
        r_unperturbed, z_unperturbed = const_args[16:18]
        eta,dpsi_perturb_val,psi_unperturbed,psi2 = const_args[18:22]

        
        path_args = Two_D_paths(folder_names=folder_name + f"sigma,tau,psi2={sigma,tau,psi2}\\")
        data_path, fig_save_path = path_args[0:2]
        video_save_path,figs_for_video_path = path_args[2:4]
        df_name_ref, fps_movie ,num_frames = path_args[4:7]

        df_name = df_name_ref + save_name_list[i] #+ f" N,ds,dt,T,tau,c0,k,kG={N,ds,dt,T,tau,c0,k,kG}"
        #start_time = time.time()
        if do_simulation == True:
            Two_d_simulation_stationary_states(
                N=N ,k=k ,c0=c0 ,sigma=sigma ,dt=dt ,ds=ds,eta=eta
                ,kG=kG ,tau=tau ,sim_steps=sim_steps, dpsi_perturb=dpsi_perturb_val
                ,L=L, r0=r0
                ,Area=Area_list
                ,psi=psi_list
                ,radi=radi_list
                ,z_list=z_list
                ,r_unperturb=r_unperturbed
                ,z_unperturb=z_unperturbed
                ,psi_unperturb=psi_unperturbed
                ,df_name = df_name
                ,num_frames = num_frames
                ,data_path = data_path
                ,Tolerence=1e-5#-5#-10
                ,save_data=True
                #,area_testing=True
            )

        #plt.show()
        #print(f"\n the simulation time={round((time.time()-start_time)/60,3)} min \n")
        if make_movie == True:
            Make_frames(
                data_path=data_path
                ,figs_save_path=figs_for_video_path
                ,df_name=df_name
                ,tot_frames= 500
            )
            Make_video(
                output_path=video_save_path
                ,input_path=figs_for_video_path
                ,video_name= df_name
                ,fps=fps_movie
            )
        
        if make_plots == True:
            plot_Epot_Ekin(
                data_path=data_path
                ,df_name=df_name
                ,output_path=video_save_path
            )
            plot_tot_area(
                data_path=data_path
                ,df_name=df_name
                ,output_path=video_save_path
            )

            #plt.show()
            plt.draw()
            plt.pause(2)
            plt.close("all")

if __name__ == "__main__":
    Surface_sim_stationary_state_initial_configuration_iterative(
        do_simulation = True
        #,start_from_flat = False
        #,do_perturbation = False #True
        ,make_movie = True
        ,make_plots= True
    )
