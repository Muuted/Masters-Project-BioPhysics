import numpy as np
import random

def gamma(i):
    gam = 1
    if i==0:
        return gam/2
    else:
        return gam

def One_D_Constants():
    """------ paths ---------"""
    save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\skole\\1 Tidligere semestre\\Kandidat speciale\\Sim data\\1D surface sim data"
    data_path =  save_path 
    fig_save_path = save_path + "\\figures and video"
    video_save_path = save_path +"\\figures and video"
    video_fig_path = save_path +"\\figures for video"

    
    """------ constants ---------"""
    L = 10e-6 # micrometers  :  Total length of line
    r0 = 0.5e-6 # micrometer  :   radius of hole
    N = 3 + 1 # Number of chain links
    m = 1e-6 # grams  :   Mass of each chain link
    gamma = 1#e-6 # units=   : The drag coefficient
    ds =  L/(N-1) # micrometers  :  Length of each chain
    T = 10 # s  : total time simulated
    dt = 1e-7 # s time step.
    k = 1 #    :  Mean curvature modulus
    kG = 1 #   :  Guassian curvature modulus
    c0 = 1 #   :  
    
    """------ variables list ---------"""
    # list of variables
    # we are gonna assume for now that the membrane is just initially flat.
    psi_list = np.zeros(shape=(T,N)) # all the angles are just flat
    psi_list[0][0] = 3.14/10
    for i in range(len(psi_list[0])):
        if i%2 == 0:
            a = -1
        else:
            a = 1
        psi_list[0][i] = (3.14/10)*a


    args = [
        L,r0,N,ds,T,dt
        ,psi_list 
        ,k,c0
        ,save_path, data_path, fig_save_path
        ,video_save_path,video_fig_path
        ]


    return args

if __name__ == "__main__":   
    One_D_Constants()
