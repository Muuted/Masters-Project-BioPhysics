import numpy as np
from Two_D_constants import Two_D_Constants, gamma, Two_D_paths
from Two_D_functions import Langrange_multi,dpsidt_func,drdt_func,dzdt_func
from Two_D_functions import Epsilon_values, c_diff_f,c_diff_g
from Two_D_functions import Epsilon_v2, c_diff
from two_d_data_processing import check_area, tot_area
import matplotlib.pyplot as plt
import os
import pandas as pd
import progressbar
import time

def Two_D_simulation(
    N:int,k:float,c0:float, dt:float, ds:float
    ,sigma:float,kG:float,tau:float
    ,sim_steps:int, L:float, r0:float
    ,radi:list,z_list:list
    ,Area:list,psi:list
    ,df_name:str,num_frames:str
    ,data_path:str
    ,save_data:bool = True 
    ,condition = 1
    ):

    
    lambs_save = []
    nus_save = []

    print("Simulation progressbar \n ")
    b = progressbar.ProgressBar(maxval=sim_steps-1)
    b.start()
    for t in range(sim_steps-1):
        b.update(t)
        lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area
        ,psi=psi[t]
        ,radi=radi[t]
        ,z_list=z_list[t]
            )
        lambs_save.append(lambs)
        nus_save.append(nus)
        for i in range(N+1):
            if i == N:
                z_list[t+1][i] = z_list[t][i]
                radi[t+1][i] = radi[t][i]
                psi[t+1][i] = psi[t][i]
            if i < N:
                z_list[t+1][i] = z_list[t][i] + dt*dzdt_func(i=i,Area=Area,radi=radi[t],nu=nus)

                radi[t+1][i] = radi[t][i] + dt*drdt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )

                psi[t+1][i] = psi[t][i] + dt*dpsidt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )
                
            if t%condition==0 and i < N:
                ebf, ebg = Epsilon_values(
                    N=N, r=radi[t+1], z=z_list[t+1] ,psi=psi[t+1] ,Area=Area
                            )
                
                K_r,K_z,K_psi = 0,0,0
                for beta in range(N):
                    ebf_val, ebg_val = ebf[beta] ,ebg[beta]
                    
                    K_r += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="r") 
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="r")
                        )
                    
                    K_z += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="z")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="z")
                        )
                    
                    K_psi += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="psi")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="psi")
                        )
                    
                radi[t+1][i] += K_r
                z_list[t+1][i] += K_z
                psi[t+1][i] += K_psi

    print("\n")
    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi],
            "r":[radi],
            "z":[z_list],
            "area list":[Area],
            'lambs': [lambs_save],
            'nus': [nus_save],
            "L" : L,
            "r0": r0,
            "N": N,
            "c0": c0,
            "k":k,
            "kG":kG,
            "sigma":sigma,
            "tau":tau,
            "sim_steps": sim_steps,
            "dt":dt,
            "ds":ds,
            "gam(i=0)":gamma(0),
            "gam(i>0)":gamma(5)
                        })

        #print(df.info())W


        if not os.path.exists(data_path):
            os.makedirs(data_path)
        df.to_pickle(data_path + df_name)




def Two_D_simulation_V2(
    N:int,k:float,c0:float, dt:float, ds:float
    ,sigma:float,kG:float,tau:float
    ,sim_steps:int, L:float, r0:float
    ,radi:list,z_list:list
    ,Area:list,psi:list
    ,df_name:str,num_frames:str
    ,data_path:str
    ,save_data:bool = True 
    ,condition = 1
    ,Tolerence = 1e-15
    ):

    Area_old = np.sum(Area)
    Area_new = 0
    correction_count = 0
    pos_count = 0
    lambs_save, nus_save = [], []

    print("Simulation progressbar \n ")
    b = progressbar.ProgressBar(maxval=sim_steps-1)
    b.start()
    for t in range(sim_steps-1):
        b.update(t)
        lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area
        ,psi=psi[t]
        ,radi=radi[t]
        ,z_list=z_list[t]
            )
        lambs_save.append(lambs)
        nus_save.append(nus)
        for i in range(N+1):
            if i == N:
                z_list[t+1][i] = z_list[t][i]
                radi[t+1][i] = radi[t][i]
                psi[t+1][i] = psi[t][i]
            if i < N:
                z_list[t+1][i] = z_list[t][i] + dt*dzdt_func(i=i,Area=Area,radi=radi[t],nu=nus)

                radi[t+1][i] = radi[t][i] + dt*drdt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )

                psi[t+1][i] = psi[t][i] + dt*dpsidt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )
            
        Area_new = tot_area(N=N,r=radi[t+1],z=z_list[t+1])
        dA = np.abs(Area_new - Area_old)
        pos_count += 1
        if Tolerence < dA:
            correction_count += 1
            for i in range(N):
                ebf, ebg = Epsilon_values(
                    N=N, r=radi[t+1], z=z_list[t+1] ,psi=psi[t+1] ,Area=Area
                            )
                K_r,K_z,K_psi = 0,0,0
                for beta in range(N):
                    ebf_val, ebg_val = ebf[beta] ,ebg[beta]
                    
                    K_r += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="r") 
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="r")
                        )
                    
                    K_z += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="z")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="z")
                        )
                    
                    K_psi += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="psi")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="psi")
                        )
                    
                radi[t+1][i] += K_r
                z_list[t+1][i] += K_z
                psi[t+1][i] += K_psi

    print("\n")
    #print(f"correction count={correction_count} of {sim_steps*(N-1)} possible")
    print(f"correction count={correction_count} of {pos_count} possible")

    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi],
            "r": [radi],
            "z": [z_list],
            "area list": [Area],
            'lambs': [lambs_save],
            'nus': [nus_save],
            "L" : L,
            "r0": r0,
            "N": N,
            "c0": c0,
            "k": k,
            "kG": kG,
            "sigma": sigma,
            "tau": tau,
            "sim_steps": sim_steps,
            "dt": dt,
            "ds": ds,
            "gam(i=0)": gamma(0),
            "gam(i>0)": gamma(5),
            "correction count": correction_count
                        })


        if not os.path.exists(data_path):
            os.makedirs(data_path)
        df.to_pickle(data_path + df_name)




def Two_D_simulation_V3(
    N:int,k:float,c0:float, dt:float, ds:float
    ,sigma:float,kG:float,tau:float
    ,sim_steps:int, L:float, r0:float
    ,radi:list,z_list:list
    ,Area:list,psi:list
    ,df_name:str,num_frames:str
    ,data_path:str
    ,save_data:bool = True 
    ,condition = 1
    ,Tolerence = 1e-15
    ):

    Area_old = np.sum(Area)
    Area_new = 0
    correction_count = 0
    pos_count = 0
    lambs_save, nus_save = [], []

    print("Simulation progressbar \n ")
    b = progressbar.ProgressBar(maxval=sim_steps-1)
    b.start()
    for t in range(sim_steps-1):
        b.update(t)
        lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area
        ,psi=psi[t]
        ,radi=radi[t]
        ,z_list=z_list[t]
            )
        lambs_save.append(lambs)
        nus_save.append(nus)
        for i in range(N+1):
            if i == N:
                z_list[t+1][i] = z_list[t][i]
                radi[t+1][i] = radi[t][i]
                psi[t+1][i] = psi[t][i]
            if i < N:
                z_list[t+1][i] = z_list[t][i] + dt*dzdt_func(i=i,Area=Area,radi=radi[t],nu=nus)

                radi[t+1][i] = radi[t][i] + dt*drdt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )

                psi[t+1][i] = psi[t][i] + dt*dpsidt_func(
                                                    i=i
                                                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                                    ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                                    ,lamb=lambs,nu=nus
                                                    )
            
        Area_new = tot_area(N=N,r=radi[t+1],z=z_list[t+1])
        dA = np.abs(Area_new - Area_old)
        do_correction = False
        if Tolerence < dA:
            do_correction = True
        pos_count += 1
        while do_correction == True:
            correction_count += 1
            for i in range(N):
                ebf, ebg = Epsilon_values(
                    N=N, r=radi[t+1], z=z_list[t+1] ,psi=psi[t+1] ,Area=Area
                            )
                K_r,K_z,K_psi = 0,0,0
                for beta in range(N):
                    ebf_val, ebg_val = ebf[beta] ,ebg[beta]
                    
                    K_r += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="r") 
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="r")
                        )
                    
                    K_z += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="z")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="z")
                        )
                    
                    K_psi += (
                        ebf_val*c_diff_f(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],Area=Area,diff_var="psi")
                        + ebg_val*c_diff_g(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="psi")
                        )
                    
                radi[t+1][i] += K_r
                z_list[t+1][i] += K_z
                psi[t+1][i] += K_psi

            Area_new = tot_area(N=N,r=radi[t+1],z=z_list[t+1])
            dA = np.abs(Area_new - Area_old)
            if Tolerence > dA :
                do_correction = False
                break
    
    print("\n")
    #print(f"correction count={correction_count} of {sim_steps*(N-1)} possible")
    print(f"correction count={correction_count} of {pos_count} possible")

    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi],
            "r": [radi],
            "z": [z_list],
            "area list": [Area],
            'lambs': [lambs_save],
            'nus': [nus_save],
            "L" : L,
            "r0": r0,
            "N": N,
            "c0": c0,
            "k": k,
            "kG": kG,
            "sigma": sigma,
            "tau": tau,
            "sim_steps": sim_steps,
            "dt": dt,
            "ds": ds,
            "gam(i=0)": gamma(0),
            "gam(i>0)": gamma(5),
            "correction count": correction_count
                        })


        if not os.path.exists(data_path):
            os.makedirs(data_path)
        df.to_pickle(data_path + df_name)






def Two_d_simulation_stationary_states(
    N:int,k:float,c0:float, dt:float, ds:float
    ,sigma:float,kG:float,tau:float
    ,sim_steps:int, L:float, r0:float
    ,radi:list,z_list:list
    ,Area:list,psi:list
    ,df_name:str,num_frames:str
    ,data_path:str
    ,save_data:bool = True 
    ,Tolerence = 1e-15
    ,do_correction = True
    ,area_testing = False
    ):
    np.set_printoptions(legacy='1.25')
    start_time = time.time()

    Area_initial = np.sum(Area)
    Area_new = 0
    
    lambs_save, nus_save = [], []
    correct_count_list = np.zeros(sim_steps-1)
    print("Simulation progressbar \n ")
    b = progressbar.ProgressBar(maxval=sim_steps-1)
    b.start()
    for t in range(sim_steps-1):
        b.update(t)
        t1,t2 = t%2, (t+1)%2
        
        lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area
        ,psi=psi[t]
        ,radi=radi[t]
        ,z_list=z_list[t]
            )
        lambs_save.append(lambs)
        nus_save.append(nus)
        for i in range(N+1):
            if i == N:
                z_list[t+1][i] = z_list[t][i]
                radi[t+1][i] = radi[t][i]
                #psi[t+1][i] = psi[t][i]
            if i < N:
                z_list[t+1][i] = z_list[t][i] + dt*dzdt_func(i=i,Area=Area,radi=radi[t],nu=nus)

                drdt =drdt_func(
                            i=i
                            ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                            ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                            ,lamb=lambs,nu=nus
                            )
                radi[t+1][i] = radi[t][i] + dt*drdt

                dpsidt = dpsidt_func(
                                i=i
                                ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                                ,Area=Area,psi=psi[t],radi=radi[t],z_list=z_list[t]
                                ,lamb=lambs,nu=nus
                                )
                psi[t+1][i] = psi[t][i] + dt*dpsidt
            
        Area_new = tot_area(N=N,r=radi[t+1],z=z_list[t+1])
        dA = Area_new - Area_initial #+1
        
        if area_testing == True:
            """ start: Lists and variables for problem finding"""
            Area_compare = [dA]
            #Area_compare.append(dA)
            total_Area_change = [Area_initial, Area_new]
        
            r_before = [i for i in radi[t+1]]
            z_before = [i for i in z_list[t+1]]

            r_change_values_list = [[] for i in range(N)]
            z_change_values_list = [[] for i in range(N)]
            psi_change_values_list = [[] for i in range(N)]

            """ end: Lists and variables for problem finding"""
        
        correction_count = 0
        while Tolerence < abs(dA) and do_correction == True:
            correction_count += 1
            #print(f"correction count={correction_count}",end="\r")
            epsilon = Epsilon_v2(
                    N=N, r=radi[t+1], z=z_list[t+1] ,psi=psi[t+1] ,Area=Area
                            )
            scaleing = 1
            for i in range(N):      
                K_r,K_z,K_psi = 0,0,0
                for beta in range(2*N):
                    
                    K_r += (
                        epsilon[beta]*c_diff(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="r")
                        )*scaleing
                    
                    K_z += (
                        epsilon[beta]*c_diff(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="z")
                        )*scaleing
                    
                    K_psi += (
                        epsilon[beta]*c_diff(i=beta,j=i,N=N,r=radi[t+1],psi=psi[t+1],z=z_list[t+1],Area=Area,diff_var="psi")
                        )*scaleing
                    
                radi[t+1][i] += K_r
                z_list[t+1][i] += K_z
                psi[t+1][i] += K_psi

                if area_testing == True:
                    r_change_values_list[i].append(K_r)
                    z_change_values_list[i].append(K_z)
                    psi_change_values_list[i].append(K_psi)                

                    if correction_count == 1:
                        r_after = [ i for i in radi[t+1]]
                        z_after = [ i for i in z_list[t+1]]

            Area_new = tot_area(N=N,r=radi[t+1],z=z_list[t+1])
            dA = Area_new - Area_initial

            if area_testing == True:
                total_Area_change.append(Area_new)
                Area_compare.append(dA)

            corr_max = 100
            if correction_count >= corr_max:
                print(f"{corr_max} corrections, is too many corrections, we close the program. ")
                print(
                    f"\n At interuption: the simulation time={round((time.time()-start_time)/60,3)} min \n"
                    +f"And t={t} of {sim_steps} \n"
                    +f"and time={t*dt} [s] at stopping")
                plt.figure()
                plt.plot(radi[t],z_list[t],"o-",label="at t")
                plt.plot(radi[t+1],z_list[t+1],"o-",label="at t+1")
                plt.legend()
                plt.show()
                #break
                exit()
            
        
        correct_count_list[t] = correction_count

        if area_testing == True:
            print(f"Shape of r change values list = {np.shape(r_change_values_list)}")
            print(f"dA[end] ={Area_compare[len(Area_compare)-1]} and num correct count ={correction_count}")
            print(f"Area_compare[0]={Area_compare[0]}")
            #print(f"list Area compare = {Area_compare}")
            print(f" Checking the dA \n"
                +f"dA before = {tot_area(N=N,r=r_before,z=z_before) -Area_initial} \n"
                +f"dA After = {tot_area(N=N,r=r_after,z=z_after)-Area_initial} \n"
                +f" Checking the total area \n"
                +f"Initial Area = {Area_initial} \n"
                +f"Area no correction = {tot_area(N=N,r=r_before,z=z_before) } \n"
                +f"Area one correction = {tot_area(N=N,r=r_after,z=z_after)}"
                )

            plt.figure()
            font_size = 15
            #plt.plot(0,Area_old,"o",label=r"$Area_{init}$")
            plt.plot(Area_compare,".-")
            
            plt.title("testing of correction dA = Area_new - Area_initial \n" +f"tolerence={Tolerence} and dA_final={dA}",fontsize=font_size)
            plt.xlabel(f"number of corrections for specific t={t}",fontsize=font_size)
            plt.ylabel("dA",fontsize=font_size)


            plt.figure()
            plt.plot(r_before,z_before,"o-",label="No correction")
            plt.plot(r_after,z_after,".-",label="1 correction")
        
            plt.title("check correction")
            plt.legend()

            plt.figure()
            plt.plot(total_Area_change,".-",label="total area")
            plt.title(
                f"Total area of membrane during corrections \n"
                +f"max={max(total_Area_change)} ,min={min(total_Area_change)} \n " 
                +f"and max/min={max(total_Area_change)/min(total_Area_change)}"
                )
            
            fig, ax = plt.subplots(1,3)
            for i in range(N):
                ax[0].plot(r_change_values_list[i],".-",label=f"linknum={i}")
                ax[0].set_title(f"K_r")
                ax[0].legend()

                ax[1].plot(z_change_values_list[i],".-",label=f"linknum={i}")
                ax[1].set_title(f"K_z")
                ax[1].legend()

                ax[2].plot(psi_change_values_list[i],".-",label=f"linknum={i}")
                ax[2].set_title(f"K_psi")
                ax[2].legend()
            plt.show()
            #exit()
    
    b.finish()
    print("\n")
    print(f"\n the simulation time={round((time.time()-start_time)/60,3)} min \n")
    plt.figure()
    plt.plot(correct_count_list,".-")
    plt.show()
    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi],
            "r": [radi],
            "z": [z_list],
            "area list": [Area],
            'lambs': [lambs_save],
            'nus': [nus_save],
            "L" : L,
            "r0": r0,
            "N": N,
            "c0": c0,
            "k": k,
            "kG": kG,
            "sigma": sigma,
            "tau": tau,
            "sim_steps": sim_steps,
            "dt": dt,
            "ds": ds,
            "gam(i=0)": gamma(0),
            "gam(i>0)": gamma(5),
            "correction count": [correction_count]
                        })


        if not os.path.exists(data_path):
            os.makedirs(data_path)
        df.to_pickle(data_path + df_name)

if __name__ == "__main__":
    pass
    exit()
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    #sim_steps = 1000

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    Two_D_simulation(
        N=N,k=k,c0=c0,sigma=sigma, dt=dt
        ,kG=kG,tau=tau, sim_steps=sim_steps
        ,L=L, r0=r0
        ,Area=Area_list
        ,psi=psi_list
        ,radi=radi_list
        ,z_list=z_list
        ,df_name = df_name
        ,num_frames = num_frames
        ,data_path = data_path
    )

    df_sim = pd.read_pickle(data_path + df_name)
    r = df_sim['r'][0]        
    Area = df_sim['area list'][0]

    error = False
    for t in range(sim_steps):
        error = check_area(
            t=t,N=N,r=r[t],Area=Area
        )
        if error == True:
            print(
                f"\n ----------------------------- \n"
                +f"error at dt={dt}"
                +f"\n ----------------------------- \n"
                )
            break
    i = 0
    Area_change = []
    for t in range(sim_steps):
        dA = np.pi*( r[t][i+1]**2 - r[t][i]**2 ) - Area[i] 
        Area_change.append(dA)

    #print(Area_change)
    print(len(Area_change))
    plt.figure()
    plt.plot(Area_change[0:len(Area_change)-1],'-o')
    plt.show()
    