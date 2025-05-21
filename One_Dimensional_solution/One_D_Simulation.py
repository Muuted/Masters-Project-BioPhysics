import numpy as np
import matplotlib.pyplot as plt
from One_D_Constants import One_D_Constants,gamma
from One_D_Functions import *
from plotting_functions import plot_from_psi

def sim_1D_surface():

    args = One_D_Constants()

    L,r0,N,ds,T,dt = args[0:6]
    x_list ,z_list ,psi_list =args[6:9]
    lambda_list ,nu_list = args[9:11]
    k,c0 = args[11:13]
    sim_T_tot = int(T/dt)

    for time in range(T-1):
        print("time:",time)
        multipliers = Lagran_multi(
            psi_list=psi_list
            ,t=time,k=k,c0=c0,ds=ds
        )[0]
        
        for link in range(N-2,-1,-1):
            dPsidt(
                i=link, t=time, dt=dt, deltaS=ds
                ,multi=multipliers
                ,psi=psi_list
            )


    x0,y0 = plot_from_psi(psi=psi_list[0],ds=ds,r0=r0)
    x1,y1 = plot_from_psi(psi=psi_list[1],ds=ds,r0=r0)
    xT,yT = plot_from_psi(psi=psi_list[T-1],ds=ds,r0=r0)
    plt.figure()
    plt.plot(x0,y0,'-*',label="t=0")
    plt.plot(x1,y1,'-.',label="t=1")
    plt.plot(xT,yT,'--',label="t=T")

    plt.legend()
    plt.show()




if __name__ == "__main__":
    sim_1D_surface()