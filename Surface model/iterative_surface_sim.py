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
        ,make_movie: bool = False#True
        ,make_plots: bool = True
    ):
    print("\n Now Running the surface simulation from stationary configurations \n")

    tau_list = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    psi2_L_list = [
        -0.3264e-4,-0.1634e-4,-0.0522e-4,-0.0322e-4,-0.0207e-4,-0.0139e-4,-0.0096e-4,-0.0050e-4
        ,-0.0037e-4,-0.0028e-4,-0.0021e-4,-0.0017e-4,-0.0013e-4,-0.0010e-4,-0.0008e-4,-0.0007e-4
        ,-0.0005e-4,-0.0004e-4,-0.0003e-4,-0.0003e-4,-0.0002e-4,-0.0002e-4,-0.0002e-4,-0.0001e-4
    ]

    sigma_list = [
        -0.4000,-0.3696,-0.3089,-0.2785,-0.2481,-0.2177,-0.1873,-0.1266,-0.0962
        ,-0.0658,-0.0354,-0.0051,0.0253,0.0557,0.0861,0.1165
        ,0.1468,0.1772,0.2076,0.2380,0.2684,0.2987,0.3291,0.3595
    ]
    
    
    save_name_list = ["-perturbed","+perturbed","unperturbed"]
    dpsi_list = [ -0.01, 0.01 ,0]
    do_perturbation_list = [True ,True ,False]

    for i in range(len(save_name_list)):
        const_args = Two_D_Constants_stationary_state(
            print_val = True
            ,show_stationary_state = True
            ,start_flat = start_from_flat
            ,perturb = do_perturbation_list[i]
            ,dpsi_perturb = dpsi_list[i]
        )

        L,r0,N,ds,T,dt = const_args[0:6]
        k,c0,sim_steps = const_args[6:9]
        sigma, tau, kG = const_args[9:12]
        Area_list, psi_list = const_args[12:14]
        radi_list,z_list = const_args[14:16]
        r_unperturbed, z_unperturbed = const_args[16:18]
        eta,dpsi_perturb_val,psi_unperturbed = const_args[18:21]
        psi2_init ,alpha = const_args[21:23]
        
        path_args = Two_D_paths(folder_names = save_name_list[i] + f" sigma,tau,psi2=({sigma:0.1e},{tau:0.1e},{psi2_init:0.1e})\\")
        data_path, fig_save_path = path_args[0:2]
        video_save_path,figs_for_video_path = path_args[2:4]
        df_name_ref, fps_movie ,num_frames = path_args[4:7]

        df_name = df_name_ref + f" N,ds,dt,T,tau,c0={N,ds,dt,T,tau,c0}"
        

        print(f"\n i={i} of {len(save_name_list)-1} \n ")
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
                ,tot_frames= 250
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
        do_simulation = False#True
        #,start_from_flat = False
        #,do_perturbation = False #True
        ,make_movie = True
        ,make_plots= True
    )
