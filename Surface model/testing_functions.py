from Two_D_constants import Two_D_Constants, Two_D_paths
from Two_D_simulation_function import Two_D_simulation
from Make_movie import Make_frames, Make_video
from two_d_data_processing import tot_area, E_pot, E_kin
from Two_D_functions import Langrange_multi, Epsilon_values
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def test_Lagrange_multi():
    const_args = Two_D_Constants(
        print_val=True
        ,init_rand_psi=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    make_movies = False
    """lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area
        ,psi=psi_list[0]
        ,radi=radi_list[0]
        ,z_list=z_list[0]
        ,print_matrix=True
    )"""
    t = 0
    lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area
        ,psi=psi_list[t]
        ,radi=radi_list[t]
        ,z_list=z_list[t]
            )

def test_make_frames():
    const_args = Two_D_Constants(
        print_val=True
    )
    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]


    Make_frames(
        data_path=data_path
        ,figs_save_path=video_fig_path
        ,df_name=df_name
    )

def test_make_video():
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    Make_video(
        output_path = video_save_path
        ,input_path = video_fig_path
        ,video_name = df_name
        ,fps=12
    )

def test_epsilon_value():
    const_args = Two_D_Constants(
        print_val=False#True
        #,init_rand_psi=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    print(f"r={radi_list[0]}")
    print(f"psi={psi_list[0]}")
    print(f"Area={Area}")
    ef,eg = Epsilon_values(
        N=N,r=radi_list[0],z=z_list[0],psi=psi_list[0],Area=Area
        ,print_matrix=True
    )

    print(f"ef={ef}")
    print(f"eg={eg}")
    
def test_tot_area():

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    dt = df_sim["dt"][0]
    sim_steps = df_sim["sim_steps"][0]

    Area_change= np.zeros(sim_steps)
    time = np.zeros(sim_steps)
    for t in range(sim_steps):
        Area_change[t] += tot_area(N=N,r=r[t],z=z[t])
        time[t] = t*dt

    Amin, Amax = min(Area_change) ,max(Area_change)
    Aratio = Amax/Amin 
    fig, ax=plt.subplots()
    plt.plot(time,Area_change,'.')
    plt.xlabel("time [s]")
    plt.ylabel("total area")
    plt.title(
        f"Ratio of Amax=Amin={Aratio} \n "
        +f"Amax - AMin={Amax-Amin}"
        )
    ax.ticklabel_format(useOffset=False)
    
    dA = np.zeros(sim_steps-1)
    for t in range(sim_steps-1):
        dA[t] = Area_change[t+1] - Area_change[t]

    fig, ax = plt.subplots()
    plt.plot(dA[0:sim_steps-1],'.')
    ax.ticklabel_format(useOffset=False)
    plt.title("Change in Area")
    plt.show()
    
def testing_Epot_Ekin():
    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    const_args = Two_D_Constants(
        print_val=True
    )
    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    #exit()
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    psi = df_sim['psi'][0]
    Area = df_sim['area list'][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    dt = df_sim["dt"][0]
    sim_steps = df_sim["sim_steps"][0]

    T = []
    S = []

    for t in range(sim_steps-1):
        T.append(
            E_kin(N=N,t=t,dt=dt,r=r,z=z,Area=Area)
            )
        S.append(
            E_pot(N=N,k=k,kG=kG,sigma=sigma,tau=tau,c0=c0
                  ,r=r[t],z=z[t],psi=psi[t],Area=Area)
        )


    plt.figure()
    plt.plot(T)
    plt.title("E kin")

    plt.figure()
    plt.plot(S)
    plt.title("E pot")

    plt.show()

if __name__ == "__main__":
    #test_Lagrange_multi()
    #test_make_frames()
    #test_make_video()
    #test_check_area()
    #test_epsilon_value()
    #test_tot_area()
    testing_Epot_Ekin()
