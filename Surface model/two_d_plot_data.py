from Two_D_constants import Two_D_Constants, Two_D_paths
from Two_D_simulation_function import Two_D_simulation
from Make_movie import Make_frames, Make_video
from two_d_data_processing import tot_area, E_pot, E_kin
from Two_D_functions import Langrange_multi, Epsilon_values
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


    
def plot_tot_area(data_path="",df_name=""):

    if data_path == "" and df_name== "":
        print(f" No paths were given in the plot_tot_area function")
        exit()

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
    plt.plot(time,Area_change,'.-')
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
    plt.plot(dA[0:sim_steps-1],'.-')
    ax.ticklabel_format(useOffset=False)
    plt.title("Change in Area")
    plt.draw()
    #plt.show()
    
def plot_Epot_Ekin(
        data_path="",df_name=""
        ):

    if data_path == "" and df_name== "":
        print(f" No paths were given, in the plot_Epot_Ekin function")
        exit()

    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    #exit()
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    psi = df_sim['psi'][0]
    Area = df_sim['area list'][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    k = df_sim['k'][0]
    kG = df_sim['kG'][0]
    sigma = df_sim['sigma'][0]
    tau = df_sim['tau'][0]
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

    t_vec = [dt*i for i in range(sim_steps-1)]
    fontsize = 15
    
    fig, ax = plt.subplots(2,1)
    ax[0].plot(t_vec,T,".-",label="Kinetic energy")
    ax[0].set_xlabel("time [s]",fontsize=fontsize)
    ax[0].set_ylabel(r"$E_{kin}$",fontsize=fontsize)
    ax[0].set_title("Kinetic energy, not scale properly \n" +f"min(E_kin)={min(T)}",fontsize=fontsize)
    ax[0].legend()
    #plt.figure()
    ax[1].plot(t_vec,S,".-",label="Potential energy")
    ax[1].set_xlabel("time [s]",fontsize=fontsize)
    ax[1].set_ylabel(r"$E_{pot}$",fontsize=fontsize)
    ax[1].set_title("Potential energy \n" +f"min(E_pot)={min(S)} \n" +r"$\Delta E_{pot}$="+f"{max(S)-min(S)}",fontsize=fontsize)
    ax[1].legend()

    plt.draw()
    #plt.show()





if __name__ == "__main__":
    #plot_tot_area()
    plot_Epot_Ekin()