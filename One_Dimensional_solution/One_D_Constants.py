import numpy as np
import random

def gamma(i,N=31):
    if N == 0:
        gam = 100
        if i==0:
            return gam/2
        else:
            return gam
    if N >0:
        b = 50
        a= 50/N
        return a*i + b

def One_D_Constants(
        print_val=False
        ):
    """------ paths ---------"""
    save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\skole\\1 Tidligere semestre\\Kandidat speciale\\Sim data\\1D surface sim data\\"
    data_path =  save_path 
    fig_save_path = save_path + "figures and video\\"
    video_save_path = save_path +"figures and video\\"
    video_fig_path = save_path +"figures for video\\"

    
    """------ constants ---------"""
    L = 10e-6 # micrometers  :  Total length of line
    r0 = 0.5e-6 # micrometer  :   radius of hole
    N = 30 + 1 # Number of chain links
    m = 1e-6 # grams  :   Mass of each chain link
    ds =  L/(N-1) # micrometers  :  Length of each chain
    T = 1e0 # s  : total time simulated
    dt = 1e-4 # s time step.
    sim_steps = int(T/dt) # : number of simulation steps
    k = 1e-5 #    :  Mean curvature modulus
    kG = 1 #   :  Guassian curvature modulus
    c0 = 1e-6 #   :  
    
    
    """------ variables list ---------"""
    # list of variables
    # we are gonna assume for now that the membrane is just initially flat.
    psi_list = np.zeros(shape=(sim_steps,N)) # all the angles are just flat
    psi_list[0][0] = 3.14/10
    for i in range(len(psi_list[0])):
        if i%2 == 0:
            a = -1
        else:
            a = 1
        psi_list[0][i] = (3.14/10)*a

    if print_val == True:
        print(
            f" \n \n"
            + "------------- Constant used in Simulation -------------- \n "
            + f"    Total length of surface L={L:e} meter \n "
            + f"    number of chain links : {N} \n " 
            + f"    r0 = {r0} \n "
            + f"    ds = {ds:e} m \n "
            + f"    dt = {dt:e} s \n "
            + f"    Total sim time = {T} s \n "
            + f"    Sim steps = {sim_steps} \n "
            + f" ------------------------------------------------------ \n \n "
        )
        
    args = [
        L,r0,N,ds,T,dt
        ,psi_list 
        ,k,c0,sim_steps
        ,save_path, data_path, fig_save_path
        ,video_save_path,video_fig_path
        ]


    return args

if __name__ == "__main__":   
    One_D_Constants()
