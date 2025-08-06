import numpy as np
import random

def gamma(i):
    gam = 1e0# standard 1e0
    if i==0:
        return gam/2
    if i > 0:
        return gam

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

def Two_D_paths():
    """------ paths ---------"""
    
    save_path = "2D sim results\\" + "Testing\\"
    data_path = save_path
    fig_save_path = save_path + "figures and video\\"
    video_save_path = save_path +"figures and video\\"
    figs_for_video_path = save_path +"figures for video\\"


    """------ Saved files names  ---------"""
    df_name= "2D surface"
    fps_movie = 24
    num_frames = 100

    path_args=[
        data_path, fig_save_path
        ,video_save_path,figs_for_video_path
        ,df_name, fps_movie ,num_frames
    ]

    return path_args

def Two_D_Constants(
        print_val=False
        ,init_rand_psi = False
        ):
    """------ constants ---------"""
    L = 100 #1e-6 # micrometers  :  Total length of line
    ds =  1e-1 # 0.1  e-9 #L/(N-1) # micrometers  :  Length of each chain
    r0 = 5 #50 #0.5e-6 # micrometer  :   radius of hole
    c0 = 0.25e0# 0.25e8 # 1/m   : 
    k = 1 #1e-12#  8e-20 # J    :  Mean curvature modulus
    sigma = 0.1#1 # k*c0**2 #
    tau = 0.1#1 # 0.1 #
    kG = 1 #   :  Guassian curvature modulus
     

    N = 10#25 #int(L/ds) # 99 + 1 # Number of chain links
    #m = 1e-6 # grams  :   Mass of each chain link
    T = 1 #5.45#s  : total time simulated
    dt = 2e-5 # s time step.
    sim_steps =int(1e4)# int(T/dt) # : number of simulation steps
    

    """------ variables list ---------"""
    # list of variables
    # we are gonna assume for now that the membrane is just initially flat.
    psi_list = np.zeros(shape=(sim_steps,N+1),dtype=float) # all the angles are just flat
    r_list =  np.zeros(shape=(sim_steps,N+1),dtype=float)
    z_list =  np.zeros(shape=(sim_steps,N+1),dtype=float)
    Area_list = np.zeros(N,dtype=float)
    
    for i in range(N+1):
        r_list[0][i] = r0 + i*ds
    
    for i in range(N):
        Area_list[i] =  np.pi*( r_list[0][i+1]**2 - r_list[0][i]**2 )
        if Area_list[i] == 0 :
            print(f"Area[{i}]=0")
            exit()
            
    if init_rand_psi == True:
        for i in range(N+1):
            if i%2 == 0:
                a = -1
            else:
                a = 1
            psi_list[0][i] = (3.14/10)*a
    
    
    if print_val == True:
        print(
            f" \n \n"
            + "------------- Constant used in Simulation -------------- \n "
            + f"    Total length of surface L={L} sim units \n "
            + f"    number of chain links : {N} \n " 
            + f"    r0 = {r0} sim units \n "
            + f"    k = {k}  sim units \n "
            + f"    c0 = {c0}  sim units \n "
            + f"    kG = {kG}  sim units \n "
            + f"    sigma = {sigma}  sim units \n "
            + f"    tau = {tau}  sim units \n "
            + f"    ds = {ds:0.1e} sim units \n "
            + f"    dt = {dt:0.1e} s \n "
            + f"    gamma(i!=0) = {gamma(2)} unit?  \n "
            + f"    Total sim time = {T} s \n "
            + f"    Sim steps = {sim_steps:0.1e} \n "
            + f" ------------------------------------------------------ \n \n "
        )
        
    args = [
        L,r0,N,ds,T,dt
        ,k,c0,sim_steps
        ,sigma,tau,kG
        ,Area_list
        ,psi_list, r_list,z_list
        ]

    return args

if __name__ == "__main__":   
    Two_D_Constants(print_val=True)
