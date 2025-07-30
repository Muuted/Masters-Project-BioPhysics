import numpy as np
from Two_D_constants import Two_D_Constants, gamma, Two_D_paths
from Two_D_functions import Delta_s
import matplotlib.pyplot as plt
import os
import pandas as pd
import progressbar


def check_area(
        t:int,N:int,r:list,Area:list
        ,tolerence:float=1e-10
        ):
    error = False
    #print(np.shape(Area))
    for i in range(N):
        area_change = np.pi*( r[i+1]**2 - r[i]**2 )

        if Area[i] != area_change :
            print(
                f"we had a change in area \n"
                +f"Area[{i}]={Area[i]}  and area_change[{i}]={area_change} \n"
                +f"at time_step={t}  and position i={i} \n"
                f" Delta A = {Area[i] - area_change}"
                )
            error = True
            break

    return error

def check_area_from_data(
        df_name:str
        ,data_path:str
        ):
    
    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    Area = df_sim['area list'][0]
    dt = df_sim['dt'][0]
    L = df_sim["L"][0]
    tau = df_sim["tau"][0]
    sigma = df_sim["sigma"][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    gam2 = df_sim["gam(i>0)"][0]
    sim_steps = df_sim["sim_steps"][0]
    r0 = df_sim["r0"][0]
    #T_tot = df_sim["Total time [sec]"][0]

    area_change = np.zeros(shape=(sim_steps,N),dtype=float)
    tolerence = 1e-10

    error = False
    for t in range(1,sim_steps):
        for i in range(N):
            area_change[t][i] =np.pi*( r[t][i+1]**2 - r[t][i]**2 )

            if Area[i] != area_change[t][i] :
                print(
                    f"we had a change in area \n"
                    +f"Area[{i}]={Area[i]}  and area_change[{i}]={area_change[t][i]} \n"
                    +f"at time_step={t}  and position i={i} \n"
                    f" Delta A = {Area[i] - area_change[t][i]}"
                    )
                error = True
                break#exit()
        if error == True:
            break

    """print(
        f"\n check all the areas and no change in area was found. with dt={dt}"
    )"""
    return error
    
def rotate_coords(
         df_name:str
        ,data_path:str
    ):
    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    N = df_sim["N"][0]
    sim_steps = df_sim["sim_steps"][0]

    r = df_sim['r'][0][sim_steps-1] # use this as reference maybe?

    x = df_sim['r'][0][sim_steps-1] # calling it x, initially
    z = df_sim['z'][0][sim_steps-1]
    y = np.zeros(N+1)

    pos_vec = [np.zeros(3) for i in range(N+1)]# 3 pos for x,y,z and then the number of time steps in the other direction

    for i in range(N+1):
        pos_vec[i][0] = x[i]
        pos_vec[i][1] = y[i]
        pos_vec[i][2] = z[i]

    num_angles = 20
    phi = np.zeros(num_angles)
    dphi = 2*np.pi/num_angles
    for i in range(num_angles):
        phi[i] = i*dphi
    



if __name__ == "__main__":
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    #sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    rotate_coords(
        df_name=df_name,data_path=data_path
    )

    