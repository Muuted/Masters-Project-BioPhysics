import numpy as np
#from Two_D_constants import Two_D_Constants, gamma, Two_D_paths, mass
#from Two_D_functions import Delta_s
import matplotlib.pyplot as plt
import os
import pandas as pd
import progressbar

def mass(i:int ,Area:list):
    rho = 1
    m = -1
    if i == 0:
        m = rho*Area[i]/2
    if i > 0:
        m = rho*( Area[i]+ Area[i-1])/2
    if m == -1:
        print("We have got the wrong correct mass, in the mass function")
        exit()
    return m

def check_area(
        t:int,N:int,r:list,z:list,Area:list
        ,tolerence:float=1e-10
        ):
    error = False
    #print(np.shape(Area))
    for i in range(N):
        area_change = np.pi*( r[i+1]+ r[i] )*np.sqrt( (r[i+1]- r[i])**2 + (z[i+1]- z[i])**2 )

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
    

def tot_area(
        N:int,r:list,z:list
        ):
    Area = 0
    for i in range(N):
        Area += np.pi*( r[i+1]+ r[i] )*np.sqrt( (r[i+1]- r[i])**2 + (z[i+1]- z[i])**2 )
    return Area


def E_pot_original(
        N:int, k:float, kG:float ,sigma:float 
        ,tau:float ,c0:float
        ,r:list,z:list,psi:list
        ,Area:list
        ):
    Atot = np.sum(Area)
    Epot = tau*r[0] - sigma*Atot/(2*np.pi)
    for i in range(N):
        if i < N -1:
            Epot += (
                (k*Area[i]/(2*np.pi*(r[i+1]+r[i])))*(
                    np.pi*(psi[i+1]-psi[i])*(r[i+1]+r[i])/Area[i] + np.sin(psi[i])/r[i] - c0
                    )**2*r[i]
                + sigma*Area[i]*r[i]/(np.pi*(r[i+1]+r[i]))
                + kG*(psi[i+1]-psi[i])*np.sin(psi[i])
            )
        if i == N -1 :
            Epot += (
                (k*Area[i]/(2*np.pi*(r[i+1]+r[i])))*(
                    np.pi*(-psi[i])*(r[i+1]+r[i])/Area[i] + np.sin(psi[i])/r[i] - c0
                    )**2*r[i]
                + sigma*Area[i]*r[i]/(np.pi*(r[i+1]+r[i]))
                + kG*(-psi[i])*np.sin(psi[i])
            )

    return Epot

def E_pot(
        N:int, k:float, kG:float
        ,tau:float ,c0:float
        ,r:list,psi:list
        ,Area:list
        )->float:
    Epot = tau*r[0]
    for i in range(N):
        if i < N - 1:
            Epot += (
                (k*Area[i]/(4*np.pi))*(
                    np.pi*(psi[i+1]-psi[i])*(r[i+1]+r[i])/Area[i] + np.sin(psi[i])/r[i] - c0
                    )**2
                + kG*(psi[i+1]-psi[i])*np.sin(psi[i])
            )
        if i == N - 1 :
            Epot += (
                (k*Area[i]/(4*np.pi))*(
                    -np.pi*psi[i]*(r[i+1]+r[i])/Area[i] + np.sin(psi[i])/r[i] - c0
                    )**2
                - kG*psi[i]*np.sin(psi[i])
            )

    return Epot

def E_kin(
        N:int, t:int, dt:float
        ,r:list ,z:list ,Area:list
        )->float:
    Ekin = 0
    for i in range(N):
        m = mass(i=i,Area=Area)

        dot_r = (r[t+1][i] - r[t][i])/dt
        dot_z = (z[t+1][i] - z[t][i])/dt

        Ekin += m*( dot_r**2 + dot_z**2 )

    return Ekin



def Xsqaured_test(N:int
                  ,r_init:list,z_init:list,psi_init:list
                  ,r:list,z:list,psi:list)->float:
    X = 0
    for i in range(N):
        X += (r_init[i] - r[i])**2 + (z_init[i] - z[i])**2 #+ (psi_init[i] - psi[i])**2
        #The angles are removed as they are between -pi and pi, which could completely dominate
        #the values for X given that r and z could potentially be much smaller.
    return X
    
def Excess_Area(rmax:float,Area_tot:float):
    return Area_tot - np.pi*rmax**2


if __name__ == "__main__":
    exit()
