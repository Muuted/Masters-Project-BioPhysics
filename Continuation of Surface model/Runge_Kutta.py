from Two_D_functions import dpsidt_func,drdt_func,dzdt_func
import numpy as np


def RungeKutta45(
    N:int,dt:float,k:float,c0:float, sigma:float
    ,kG:float, tau:float, ds:float,eta:float
    ,Area:list,psi_init:list,r_init:list, z_init:list
    ,lamb:list , nu:list
        ) -> list:
    
    k_psi = np.zeros(shape=(5,len(psi_init)),dtype=float)
    k_r = np.zeros(shape=(5,len(r_init)),dtype=float)
    k_z = np.zeros(shape=(5,len(z_init)),dtype=float)

    r = np.zeros(shape=(len(r_init)),dtype=float)
    z = np.zeros(shape=(len(z_init)),dtype=float)
    psi = np.zeros(shape=(len(psi_init)),dtype=float)

    for j in range(1,5):
        for i in range(N):
            if  j < 4:
                r[i] =  r_init[i] + (dt/2)*k_r[j-1][i] 
                z[i] = z_init[i] + (dt/2)*k_z[j-1][i] 
                psi[i] = psi_init[i] + (dt/2)*k_psi[j-1][i] 
            if j == 4:
                r[i] = r_init[i] + dt*k_r[j-1][i] 
                z[i] = z_init[i] + dt*k_z[j-1][i]
                psi[i] =  psi_init[i] + dt*k_psi[j-1][i] 

        for i in range(N):
            k_r[j][i] = (
                drdt_func(
                    i=i,N=N,k=k,c0=c0, sigma=sigma, kG=kG, tau=tau, ds=ds, eta=eta
                    ,Area=Area,lamb=lamb , nu=nu
                    ,z_list=z ,psi=psi ,radi=r
                    )
            )
            k_z[j][i] = (
                dzdt_func(
                    i=i,Area=Area,nu=nu ,eta=eta,ds=ds
                    ,radi=r
                    )
            )
            k_psi[j][i] = (
                dpsidt_func(
                    i=i
                    ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,eta=eta,ds=ds
                    ,Area=Area,lamb=lamb,nu=nu
                    ,radi=r,z_list=z,psi=psi 
                )
            )

    return k_r,k_z,k_psi




if __name__ == "__main__":
    pass