import pandas as pd
import time
import matplotlib.pyplot as plt
from Two_D_constants import Two_D_Constants, Two_D_paths, Two_D_Constants_stationary_state
from Two_D_simulation_function import Two_D_simulation, Two_D_simulation_V2, Two_d_simulation_stationary_states
from Make_movie import Make_frames, Make_video
from two_d_data_processing import check_area
from Two_D_functions import Langrange_multi
from two_d_plot_data import plot_Epot_Ekin, plot_tot_area
from iterative_surface_sim import Surface_sim_stationary_state_initial_configuration_iterative
import multiprocessing
import os
from Testing_ideas import plot_everytyhing
def Surface_sim_stationary_state_initial_configuration_multiprocessing(
    perturb_index:int , const_index: int 
    ):
    do_simulation:bool = True
    start_from_flat:bool = False
    make_movie: bool = False
    make_plots: bool = False
    
    print("\n Now Running the surface simulation from stationary configurations \n")

    tilde_sigma_list = [
        0.329113924050633
        ,0.29873417721519
        ,0.116455696202532 
        ]
    tilde_tau_list = [
        1.0
        ,4.47368421052632
        ,1.0
        ]
    psi2_list = [
        -1.68533976179446e-8
        ,-2.23344534748962e-6
        ,-2.26474921864332e-8
        ]

    save_name_list = ["-perturbed","+perturbed","unperturbed"]
    dpsi_list = [ -0.01, 0.01 ,0]
    do_perturbation_list = [True ,True ,False]
    if const_index == 0 or const_index == 2:
        show_print_val = [False ,False ,False]
    else:
        show_print_val = [False ,False ,True]

    const_args = Two_D_Constants_stationary_state(
        print_val = show_print_val[perturb_index]
        ,show_stationary_state = show_print_val[perturb_index]
        ,start_flat = start_from_flat
        ,perturb = do_perturbation_list[perturb_index]
        ,dpsi_perturb = dpsi_list[perturb_index]
        ,tilde_sigma = tilde_sigma_list[const_index]
        ,tilde_tau = tilde_tau_list[const_index]
        ,psi_L = psi2_list[const_index]
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]
    r_unperturbed, z_unperturbed = const_args[16:18]
    eta,dpsi_perturb_val,psi_unperturbed = const_args[18:21]
    psi2_init ,alpha = const_args[21:23]

    
    path_args = Two_D_paths(folder_names = f"T,dt={T,dt} sigma,tau,psi2=({sigma:0.1e},{tau:0.1e},{psi2_init:0.1e})\\"+save_name_list[perturb_index] +"\\" )
    data_path, fig_save_path = path_args[0:2]
    video_save_path,figs_for_video_path = path_args[2:4]
    df_name_ref, fps_movie ,num_frames = path_args[4:7]

    df_name = df_name_ref + f" N,ds,dt,T,tau={N,ds,tau}"# + ".pkl"
    
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
            ,Tolerence = 1e-5
            ,save_data = True
            ,print_progress = show_print_val[perturb_index]
        )


def main_multiProcessing(num_cpu,perturb_index_list,var_index_list):
    process = []
    for i in range(num_cpu):
        process.append(multiprocessing.Process(
            target=Surface_sim_stationary_state_initial_configuration_multiprocessing
            ,args=(perturb_index_list[i],var_index_list[i])
            ))

    for p in process:
        p.start()

    for p in process:
        p.join


    plot_everytyhing()

if __name__ == "__main__":
    cpu_6 = True
    const_index_1 = 0
    const_index_2 = 1

    perturb_list = [0,1,2,0,1,2]
    const_vars = [0,0,0,1,1,1]

    main_multiProcessing(num_cpu=6,perturb_index_list=perturb_list,var_index_list=const_vars)

    

    
