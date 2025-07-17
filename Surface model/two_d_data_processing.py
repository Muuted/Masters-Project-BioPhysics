import numpy as np
from Two_D_constants import Two_D_Constants, gamma, Two_D_paths
from Two_D_functions import Delta_s
import matplotlib.pyplot as plt
import os
import pandas as pd
import progressbar



def check_area(
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
    T_tot = df_sim["Total time [sec]"][0]

    area_change = []
    for t in range(1,sim_steps):
        pass

    





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