import pandas as pd
import time
import matplotlib.pyplot as plt
from Two_D_constants import Two_D_paths, Two_D_Constants_stationary_state
from Two_D_simulation_function import  Two_d_simulation_stationary_states

import multiprocessing
import os

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

    folder_pos_list = [
        "triangle sim\\"
        ,"plus sim\\"
        ,"cross sim\\"
    ]
    save_name_list = ["-perturbed","+perturbed","unperturbed"]
    dpsi_list = [ -0.01, 0.01 ,0]
    do_perturbation_list = [True ,True ,False]
    
    
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
        ,pause_timer=10
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]
    r_unperturbed, z_unperturbed = const_args[16:18]
    eta,dpsi_perturb_val,psi_unperturbed = const_args[18:21]
    psi2_init ,alpha = const_args[21:23]

    
    folder = folder_pos_list[const_index] +f"T,dt,sigma,tau=({T:0.1e},{dt:0.1e},{sigma:0.1e},{tau:0.1e})\\"+save_name_list[perturb_index] +"\\" 
    path_args = Two_D_paths(folder_names = folder)
    data_path, fig_save_path = path_args[0:2]
    video_save_path,figs_for_video_path = path_args[2:4]
    df_name_ref, fps_movie ,num_frames = path_args[4:7]

    df_name = df_name_ref +".pkl"#+ f"N,ds,dt,T,tau=({N},{ds},{tau:0.1e})" + ".pkl"
    
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



def Surface_sim_stationary_state_initial_configuration_multiprocessing_perturb_tau_and_psi(
    perturb_index:int , const_index: int
    ,perturb_psi:bool = True ,perturb_tau:bool = False
    ):

    do_simulation:bool = True
    start_from_flat:bool = False
    make_movie: bool = False
    make_plots: bool = False
    
    print("\n Now Running the surface simulation from stationary configurations \n")

    tilde_sigma_list = [
        0.0253164556962025 #0.329113924050633
        ,0.29873417721519
        ,0.116455696202532 
        ]
    tilde_tau_list = [
        1.31578947368421#1.0
        ,4.47368421052632
        ,1.0
        ]
    psi2_list = [
        -2.83260429562395e-7#-1.68533976179446e-8
        ,-2.23344534748962e-6
        ,-2.26474921864332e-8
        ]

    folder_pos_list = [
        "triangle sims\\"
        ,"plus sims\\"
        ,"cross sims\\"
    ]
    save_name_list = ["smaller","larger","no change"]

    dpsi_list = [ -0.01, 0.01 ,0]

    dtau_list = [1-0.05, 1+0.05 , 1]
    
    if perturb_psi == True:
        do_perturbation_list = [True ,True ,False]
    elif perturb_psi == False:
        do_perturbation_list = [False ,False ,False]

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
        ,pause_timer=10
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]
    r_unperturbed, z_unperturbed = const_args[16:18]
    eta,dpsi_perturb_val,psi_unperturbed = const_args[18:21]
    psi2_init ,alpha = const_args[21:23]

    if perturb_tau == True:
        tau *= dtau_list[perturb_index]
    
    folder = folder_pos_list[const_index] +f"T,dt,sigma,tau=({T:0.1e},{dt:0.1e},{sigma:0.1e},{tau:0.1e})\\"+save_name_list[perturb_index] +"\\" 
    path_args = Two_D_paths(folder_names = folder)
    data_path, fig_save_path = path_args[0:2]
    video_save_path,figs_for_video_path = path_args[2:4]
    df_name_ref, fps_movie ,num_frames = path_args[4:7]

    df_name = df_name_ref +".pkl"#+ f"N,ds,dt,T,tau=({N},{ds},{tau:0.1e})" + ".pkl"
    
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

def main_multiProcessing():
    num_cpu = 6
    perturb_list = [0,1,2,0,1,2]

    const_vars = [0,0,0,1,1,1]#[2,2,2]#,1,1,1] # 0 = triangle, 1 = plus and 2 = cross (in Latex)
    process = []
    for i in range(num_cpu):
        process.append(multiprocessing.Process(
            target=Surface_sim_stationary_state_initial_configuration_multiprocessing
            ,args=(perturb_list[i],const_vars[i])
            ))

    for p in process:
        p.start()

    for p in process:
        p.join


def main_multiProcessing_tau_and_psi():
    num_cpu = 3
    perturb_list = [0,1,2,0,1,2]
    perturb_tau = False
    perturb_psi = True
    const_vars = [0,0,0]#,1,1,1]#[2,2,2]#,1,1,1] # 0 = triangle, 1 = plus and 2 = cross (in Latex)
    process = []
    for i in range(num_cpu):
        process.append(multiprocessing.Process(
            target=Surface_sim_stationary_state_initial_configuration_multiprocessing_perturb_tau_and_psi
            ,args=(perturb_list[i],const_vars[i],perturb_psi,perturb_tau)
            ))

    for p in process:
        p.start()

    for p in process:
        p.join
if __name__ == "__main__":

    #main_multiProcessing()
    main_multiProcessing_tau_and_psi()

    

    
