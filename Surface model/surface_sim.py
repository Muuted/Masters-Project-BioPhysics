import pandas as pd
import time
import matplotlib.pyplot as plt
from Two_D_constants import Two_D_Constants, Two_D_paths, Two_D_Constants_stationary_state
from Two_D_simulation_function import Two_D_simulation, Two_D_simulation_V2, Two_d_simulation_stationary_states
from Make_movie import Make_frames, Make_video
from two_d_data_processing import check_area
from Two_D_functions import Langrange_multi
from two_d_plot_data import plot_Epot_Ekin, plot_tot_area

def surface_sim_find_c0():
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
        


def Surface_sim():
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]


    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,figs_for_video_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    start_time = time.time()
    Two_D_simulation(
        N=N ,k=k ,c0=c0 ,sigma=sigma ,dt=dt ,ds=ds
        ,kG=kG ,tau=tau ,sim_steps=sim_steps
        ,L=L, r0=r0
        ,Area=Area_list
        ,psi=psi_list
        ,radi=radi_list
        ,z_list=z_list
        ,df_name = df_name
        ,num_frames = num_frames
        ,data_path = data_path
    )

    print(f"\n the simulation time={(time.time()-start_time)/60} min \n")
    plot_tot_area()
    exit()
    plot_Epot_Ekin()
    
    
    Make_frames(
        data_path=data_path
        ,figs_save_path=figs_for_video_path
        ,df_name=df_name
    )
    Make_video(
        output_path=video_save_path
        ,input_path=figs_for_video_path
        ,video_name= df_name
        ,fps=fps_movie
    )
    

def Surface_sim_Area_condition():
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]


    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,figs_for_video_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    start_time = time.time()
    Two_D_simulation_V2(
        N=N ,k=k ,c0=c0 ,sigma=sigma ,dt=dt ,ds=ds
        ,kG=kG ,tau=tau ,sim_steps=sim_steps
        ,L=L, r0=r0
        ,Area=Area_list
        ,psi=psi_list
        ,radi=radi_list
        ,z_list=z_list
        ,df_name = df_name
        ,num_frames = num_frames
        ,data_path = data_path
    )

    print(f"\n the simulation time={(time.time()-start_time)/60} min \n")
    plot_tot_area()
    exit()
    plot_Epot_Ekin()
    
    
    Make_frames(
        data_path=data_path
        ,figs_save_path=figs_for_video_path
        ,df_name=df_name
    )
    Make_video(
        output_path=video_save_path
        ,input_path=figs_for_video_path
        ,video_name= df_name
        ,fps=fps_movie
    )
    

def Surface_sim_stationary_state_initial_configuration(
        do_simulation:bool = True
        ,start_from_flat:bool = False
        ,do_perturbation:bool = False
        ,make_movie: bool = True
        ,make_plots: bool = True
    ):
    print("\n Now Running the surface simulation from stationary configurations \n")

    const_args = Two_D_Constants_stationary_state(
        print_val=True
        ,show_stationary_state=True
        ,start_flat=start_from_flat
        ,perturb=do_perturbation
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]
    r_unperturbed, z_unperturbed = const_args[16:18]
    eta = const_args[18]

    
    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,figs_for_video_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    df_name += f" N,ds,dt,T={N,ds,dt,T}" # c0={c0} tau={tau}"#f" ds={dt}and N={N} and ds={ds} c0={c0} kG={kG}"
    #start_time = time.time()
    if do_simulation == True:
        Two_d_simulation_stationary_states(
            N=N ,k=k ,c0=c0 ,sigma=sigma ,dt=dt ,ds=ds,eta=eta
            ,kG=kG ,tau=tau ,sim_steps=sim_steps
            ,L=L, r0=r0
            ,Area=Area_list
            ,psi=psi_list
            ,radi=radi_list
            ,z_list=z_list
            ,r_unperturb=r_unperturbed
            ,z_unperturb=z_unperturbed
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

        plt.show()
    

def Speed_diagnosing():
    import cProfile, pstats, io
    from pstats import SortKey
    pr = cProfile.Profile()
    pr.enable()

    #surface_sim_find_c0()
    #Surface_sim()
    #Surface_sim_Area_condition()
    Surface_sim_stationary_state_initial_configuration(
        #do_simulation=False
        #start_from_flat=False
    )
    pr.disable()
    s = io.StringIO()
    sortby = SortKey.CUMULATIVE#SortKey.NAME#SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats(20)
    print(s.getvalue())
    
if __name__ == "__main__":
    #surface_sim_find_c0()
    #Surface_sim()
    #Surface_sim_Area_condition()
    Surface_sim_stationary_state_initial_configuration(
        do_simulation = True
        ,start_from_flat = False
        ,do_perturbation = False
        ,make_movie = True
        ,make_plots= True
    )

    #Speed_diagnosing()