from Two_D_functions import dpsidt_func,drdt_func,dzdt_func
import numpy as np



def RungeKutta45(
    t:int,N:int,dt:float,k:float,c0:float, sigma:float
    ,kG:float, tau:float, ds:float,eta:float
    ,Area:list,psi:list,r:list, z:list
    ,lamb:list , nu:list
        ):
    
    k_psi =[[],[],[],[]]
    k_r =  [[],[],[],[]]
    k_z = [[],[],[],[]]
    for i in range(N):
        pass






















def drdt_RK4(
        i:int,N:int,dt:float,k:float,c0:float, sigma:float, kG:float, tau:float, ds:float,eta:float
        ,Area:list,psi:list,r:list, z:list
        ,lamb:list , nu:list
        ):
    k1 = drdt_func(
        i=i,N=N,k=k,c0=c0, sigma=sigma, kG=kG, tau=tau, ds=ds, eta=eta
        ,Area=Area,psi=psi,lamb=lamb , nu=nu, z_list=z
        ,raid=r
        )
    k2 = drdt_func(
        i=i,N=N,k=k,c0=c0, sigma=sigma, kG=kG, tau=tau, ds=ds, eta=eta
        ,Area=Area,psi=psi,lamb=lamb , nu=nu, z_list=z
        ,raid=[j if j!=r[i] else j + dt*k1/2 for j in r ]
        )    
    k3 = drdt_func(
        i=i,N=N,k=k,c0=c0, sigma=sigma, kG=kG, tau=tau, ds=ds, eta=eta
        ,Area=Area,psi=psi,lamb=lamb , nu=nu, z_list=z
        ,raid=[j if j!=r[i] else j + dt*k2/2 for j in r ]
        )
    k4 = drdt_func(
        i=i,N=N,k=k,c0=c0, sigma=sigma, kG=kG, tau=tau, ds=ds, eta=eta
        ,Area=Area,psi=psi,lamb=lamb , nu=nu, z_list=z
        ,raid=[j if j!=r[i] else j + dt*k3 for j in r ]
        )    
    return (dt/6)*(k1 + k2 + k3 + k4)


def  dzdt_RK4(i:int,dt:float,Area:list,radi:list, nu:list):
    k1 = dzdt_func(i=i,Area=Area,radi=radi,nu=nu)
    return dt*k1


def dPsidt_RK4(
        i:int,N:int,k:float,c0:float,sigma:float,kG:float,tau:float,dt:float
        ,Area:list,radi:list,z_list:list
        ,lambs:list,nus:list
        ,psi:list
        ):    
    k1 = dpsidt_func(
        i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi=psi 
                )
    k2 = dpsidt_func(i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi= [j if j!=psi[i] else j + dt*k1/2 for j in psi ]
    )
    k3 = dpsidt_func(i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi=[j if j!=psi[i] else j + dt*k2/2 for j in psi ]
    )
    k4 = dpsidt_func(i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi= [j if j!=psi[i] else j + dt*k3 for j in psi ]
    )
    return (dt/6)*(k1 + k2 + k3 + k4)
