import numpy as np
import matplotlib.pyplot as plt
import progressbar
import pandas as pd

from Two_D_constants import Two_D_Constants
from Two_D_functions import Delta_s


def plot_from_psi_V2(
        #psi:list,sim_steps: int,ds: float,r0:float, L:float
        data_path,df_name
        ):

    df_sim = pd.read_pickle(data_path + df_name)
    
    dt = df_sim['dt'][0]
    N = dt = df_sim['N'][0]
    L = df_sim["L"][0]
    ds = df_sim["ds"][0]
    r0 = df_sim["r0"][0]
    sim_steps = df_sim["sim_steps"][0]
    psi = df_sim["psi"][0]

    
    x = np.full(shape=(sim_steps,N),dtype=float,fill_value=0)
    z = np.zeros(shape=(sim_steps,N),dtype=float) #Because we have N+1 points, but run dynamics on N
    print(np.shape(x))
    x[:,N-1] = L + r0
    print_list = []
    print_list2 = []
    print("\n Getting x,z progressbar")
    b2 = progressbar.ProgressBar(maxval=sim_steps)
    b2.start()
    for t in range(sim_steps):
        b2.update(t)
        tolerance= 1e-10 #1e-4
        for i in range(N-2,-1,-1):
            x[t][i] = x[t][i+1] - ds*np.cos(psi[t][i])
            z[t][i] = z[t][i+1] - ds*np.sin(psi[t][i])
            
            a = np.sqrt((x[t][i+1]-x[t][i])**2 + (z[t][i+1]-z[t][i])**2)
            if ds*(1+tolerance) <= a  <= ds*(1-tolerance) :
                print(f"error on constant length , x[t][i]={x[t][i]} and z[t][i]={x[t][i]}")
                exit()

    
    df_sim["x pos"] = [x]
    df_sim["z pos"]= [z]

    df_sim.to_pickle(data_path + df_name)
    #print("hello")
    #print(df_sim.info())
 


if __name__ == "__main__":
    args = Two_D_Constants
    
    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    df_name = args[15]
    df_sim = pd.read_pickle(data_path + df_name)
    psi_list = df_sim["Psi"][0]
    z = df_sim["z pos"][0]
    x = df_sim["x pos"][0]

    