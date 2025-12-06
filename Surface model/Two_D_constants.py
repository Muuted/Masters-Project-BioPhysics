from two_d_continues_integration import find_init_stationary_state
from Two_D_functions import Perturbation_of_inital_state,gamma
import numpy as np
import random
import matplotlib.pyplot as plt
np.set_printoptions(legacy='1.25')
"""def gamma(i):
    gam = 1e0# standard 1e0
    if i==0:
        return gam/2
    if i > 0:
        return gam"""



def Two_D_paths(folder_names=""):
    """------ paths ---------"""
    save_path =  "2D sim results\\" + "Data for thesis\\really long\\"#Data simulation\\"#"Verification\\"
    if folder_names == "":
        save_path = save_path +  "c0=0 tau=0\\"
    else:
        save_path = save_path + folder_names
    
    data_path = save_path
    fig_save_path = save_path + "figures and video\\"
    video_save_path = save_path +"figures and video\\"
    figs_for_video_path = save_path +"figures for video\\"

    
    """------ Saved files names  ---------"""
    df_name= "2D Suface"
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
    kG = 1
    sigma = k*c0**2 #
    tau = 1 #

    N = 10#25 #int(L/ds) # 99 + 1 # Number of chain links
    #m = 1e-6 # grams  :   Mass of each chain link
    T = 1 #5.45#s  : total time simulated
    dt = 1e-5 # s time step.
    sim_steps = int(1e3)# int(T/dt) # : number of simulation steps
    

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



def Two_D_Constants_stationary_state(
        print_val:bool=False
        ,show_stationary_state:bool = True
        ,pause_timer:float = 4
        ,start_flat:bool = False
        ,perturb:bool = False
        ,dpsi_perturb="",tilde_sigma=""
        ,tilde_tau="",psi_L="",N="",ds = ""
        ,perturb_tau =   ""
        )->list:
    np.set_printoptions(legacy='1.25')

    """------ constants ---------"""
    if N == "":
        N = 40#20#20 #60#20#80 #int(L/ds) # 99 + 1 # Number of chain links
    #m = 1e-6 # grams  :   Mass of each chain link
    T = 1e-7 #1e-7#1e-7#1e-7#1e-8 #5e-7 #1e-6 #3e-7# 0.3e-6# 20e-7 #10 #5.45#s  : total time simulated seconds
    dt = 1.0e-13#0.125e-12 #1e-11 #5e-11 #s time step. 
    sim_steps = int(T/dt) # : number of simulation steps
    L = 100.0 #1e-6 # micrometers  :  Total length of line
    if ds == "":
        ds = 1.5e-2/2#*(2/3)#8#0.3#1.5#e-8#e-2 #e-8 # 1.5/2#/3 #0.3 #1e-1 # 0.1  e-9 #L/(N-1) # micrometers  :  Length of each chain
    r0 = 5.0 #50 #0.5e-6 # micrometer  :   radius of hole
    if dpsi_perturb == "":
        dpsi_perturb = -0.02
    num_pertub = 10#int((N-1)/2)
    
    #Base variables
    eta = 1.0#e-3 #e-3 # SI: kg /(ms)
    c0 = 25 #0.25e8 #25#0.25#e8#0.25e0# 0.25e8 # 1/m : 
    k = 80 #8.0e-20 #1 # 8e-20 # Joule    :  Mean curvature modulus
    kG = None #-0.75*k# Joule

    # scaling parameters
    lc = 1/c0
    sigma_c = k*c0**2
    tau_c = k*c0

    #Dimless variables
    if tilde_sigma == "":
        tilde_sigma = 0.0253164556962025 #0.329113924050633  #0.29873417721519  #
    if tilde_tau == "":
        tilde_tau = 1.31578947368421 #4.47368421052632       #1 

    #Converted variables
    sigma = tilde_sigma*sigma_c
    tau = tilde_tau*tau_c 
    
    rs2 = 20*lc 
    zs2 = 0
    s0, sN = 0, 50*lc
    if psi_L == "":
        psi_L = -2.83260429562395e-7 #6.531116e-8   #-1.68533976179446e-8  

    #print(f"n={N}, ds={ds:e} , sigma={sigma} , psi2={psi_L} ,tau={tau}")
    #Creating lists for the variables.
    psi_list = np.zeros(shape=(sim_steps,N),dtype=float) # all the angles are just flat
    r_list =  np.zeros(shape=(sim_steps,N+1),dtype=float)
    z_list =  np.zeros(shape=(sim_steps,N+1),dtype=float)
    Area_list = np.zeros(N,dtype=float)


    #Initiating the inital state of the membrane
    psi,r,z, r_contin, z_contin, alpha = find_init_stationary_state(
            sigma=sigma ,k=k ,c0=c0 ,tau=tau ,ds=ds
            ,psi_L=psi_L ,r_L=rs2 ,z_L=zs2 ,s0=s0 ,sN=sN
            ,total_points = N
        )
    alpha = (c0 - (psi[1] - psi[0])/ds )*r[0]/np.sin(psi[0]) - 1
    kG = k*alpha
    

    if start_flat == True:
        for i in range(N+1):
           r_list[0][i] = r0 + i*ds
    else:
        """------ variables list ---------"""
        for i in range(N+1):
            if i < N :
                psi_list[0][i] = psi[i]
            r_list[0][i] = r[i]
            z_list[0][i] = z[i]
    
    for i in range(N+1):
        if i < N :
            Area_list[i] =  np.pi*(r_list[0][i+1] + r_list[0][i])*np.sqrt( (r_list[0][i+1] - r_list[0][i])**2 + (z_list[0][i+1] - z_list[0][i])**2 )
            if Area_list[i] == 0 :
                print(f"Area[{i}]=0")
                exit()
    
    r_unperturb = [i for i in r_list[0]]
    z_unperturb = [i for i in z_list[0]]
    psi_unperturb = [i for i in psi_list[0]]

    if perturb == True:
        Perturbation_of_inital_state(
            points_perturbed=num_pertub 
            ,ds=ds, N=N
            ,r=r_list[0]
            ,z=z_list[0]
            ,psi=psi_list[0]
            ,Area=Area_list
            ,delta_psi= dpsi_perturb
            ,show_initial_condi=True
        )
    if perturb == False:
        dpsi_perturb = 0

    if show_stationary_state==True:
        plt.figure()
        font_size = 10
        #plt.plot(r_contin,z_contin,marker="",label="integration")
        plt.plot(r_list[0],z_list[0],"o-",label="Discreet")
        plt.plot(r_list[0][0],z_list[0][0],"o",color="k",label="s1")
        plt.plot(r_list[0][len(r_list[0])-1],z_list[0][len(r_list[0])-1],"o",color="y",label="s2")
        if start_flat == False:
            plt.plot(r_contin,z_contin,linestyle="--",marker="",color="k",label="integration")
        if perturb == True:
            plt.plot(r_unperturb,z_unperturb,"-o",label="unperturbed")
        plt.xlim(min(r_list[0])-ds, max(r_list[0])+ds)
        ceil = max(r_list[0]) - min(r_list[0]) + 2*ds
        plt.ylim(-ceil/10, 9*ceil/10)
        #plt.xlim(3,83)
        #plt.ylim(-10,70)
        plt.xlabel("r",fontsize=font_size)
        plt.ylabel("z",fontsize=font_size)
        plt.title(
            f"Quick peak at the neck configuration before dynanic simulation "
            ,fontsize=font_size
            )
        plt.legend()
        #plt.xlim(min(r)*0.95, max(r)*1.05)
        #plt.ylim(-5,max(r)-min(r)-5)
        #plt.show()
        #exit()
        plt.draw()
        plt.pause(pause_timer)
        plt.close("all")


    
    #c0 = 0
    #tau = 0

    if perturb_tau != "":
        tau *= perturb_tau
    
    #print(f"end pos =(r,z)={round(r_list[0][len(r_list[0])-1],5),round(z_list[0][len(r_list[0])-1],5)}")
    if print_val == True:
        print(
            f" \n \n"
            + "------------- Constant used in Simulation -------------- \n "
            #+ f"    Total length of surface L={L} sim units \n "
            + f"    number of chain links N: {N} \n " 
            #+ f"    r0 = {r0} sim units \n "
            + f"    eta = {eta:.1e} sim units \n "
            + f"    k = {k:0.1e}  sim units \n "
            + f"    kG = {kG:.1e}  sim units \n "
            + f"    c0 = {c0:.1e}  sim units \n "
            + f"    tau = {tau:0.1e}  sim units \n "
            + f"    sigma = {sigma:0.1e}  sim units \n "
            + f"    ds = {ds:0.1e} sim units \n "
            + f"    dt = {dt:0.1e} s \n "
            + f"    eta = {eta:0.1e} s \n "
            + f"    gamma(i!=0) = {gamma(i=2,ds=ds,eta=eta)} unit?  \n "
            + f"    Total sim time = {T} s \n "
            + f"    Sim steps = {sim_steps:0.1e} \n "
            + f"    dpsi = {dpsi_perturb:0.1e} \n "
            + f"    alpha = {alpha:0.1e} \n "
            + f" ------------------------------------------------------ \n \n "
        )
    
   
    args = [
        L,r0,N,ds,T,dt
        ,k,c0,sim_steps
        ,sigma,tau,kG
        ,Area_list
        ,psi_list, r_list ,z_list
        ,r_unperturb,z_unperturb
        ,eta, dpsi_perturb , psi_unperturb
        ,psi_L, alpha
        ]

    return args

if __name__ == "__main__":   
    #Two_D_Constants(print_val=True)
    Two_D_Constants_stationary_state(
        show_stationary_state=True
        ,pause_timer=300
        ,start_flat=False
        ,perturb=False
        ,print_val=True
    )
