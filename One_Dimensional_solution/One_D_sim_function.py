import numpy as np
import matplotlib.pyplot as plt
from One_D_Constants import gamma
from One_D_Functions import Lagran_multi, dPsidt_RungeKutta_4, dPsidt
import os
import pandas as pd
import progressbar


def sim_1D_surface(
        L,r0,N,ds,T,dt
        ,psi_list,k,c0,sim_steps
        ,data_path ,df_name
        ,save_data=True
                  ):

    multipliers= []

    print("Simulation progressbar \n ")
    b = progressbar.ProgressBar(maxval=sim_steps-1)
    b.start()
    for time in range(0,sim_steps-1):
        b.update(time)
        #t1, t2 = time%2 , (time + 1)%2
        x = Lagran_multi(
                psi_list=psi_list
                ,t=time,k=k,c0=c0,ds=ds
                ,linalg_lstsq=False
                ,num_chains=N
                        )
        multipliers.append(x)

        for link in range(len(psi_list[0])-2,-1,-1):
            if link == N :
                psi_list[time+1][link] =  0
                
            if link < N:
                RungeKutta4 = dPsidt_RungeKutta_4(
                    link=link
                    ,N=N, ds=ds, dt=dt
                    ,multipliers=x#multipliers[time]
                    ,psi=psi_list[time][link]
                    )
                
                psi_list[time+1][link] = psi_list[time][link] + (dt/6)*RungeKutta4

                """dpsidt = dPsidt(
                    i=link, N=N, deltaS=ds
                    ,multi=x
                    ,psi=psi_list[time][link]
                        )
                psi_list[time+1][link] = psi_list[time][link] + dt*dpsidt"""
  
    if save_data == True:
        df = pd.DataFrame({
            'psi': [psi_list],
            'multipliers': [multipliers],
            "L" : L,
            "r0": r0,
            "N": N,
            "c0": c0,
            "k":k,
            "Total time [sec]" : T,
            "sim_steps": sim_steps,
            "dt":dt,
            "ds":ds,
            "gam(i=0)":gamma(0),
            "gam(i>0)":gamma(5)
                        })

        #print(df.info())W

        print("\n")
        if not os.path.exists(data_path):
            os.makedirs(data_path)
        df.to_pickle(data_path + df_name)


if __name__ == "__main__":
    pass