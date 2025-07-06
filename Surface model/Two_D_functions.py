import numpy as np
from Two_D_constants import Two_D_Constants, gamma
#import os
#import pandas as pd
#import progressbar


def Delta_s(A:list,r:list,i:int):
    return A[i]/(np.pi(r[i+1] +r[i]))


def Q_function(
        i,N,k,c0, sigma, kG, tau
        ,Area:list,psi:list,radi:list
        ):
    a,a1,a2,a3,a4,a5,a6,a7 = 0,0,0,0,0,0,0,0
    if i == 0:
        
        a11 = k*Area[i]/(2*np.pi*(radi[i+1]+radi[i])**2)
        a12 = (
            np.pi*(psi[i+1]-psi[i])*(radi[i+1]+radi[i])/Area[i]
            +np.sin(psi[i])/radi[i]
            - c0
        )
        a1 = a11*a12**2

        a2 = - k*(psi[i+1]-psi[i])/(radi[i+1] +radi[i])*a12

        a31 = k*Area[i]*np.sin(psi[i])/(np.pi*(radi[i+1]+radi[i])*radi[i]**2)
        a3 = a31*a12 - tau

        a4 = - sigma*Area[i]*radi[i+1]/(np.pi*(radi[i+1]+radi[i])**2)

        a5 = kG*(psi[i+1]-psi[i])*np.sin(psi[i])/radi[i]**2


    if 0 < i < N-1:
        
        a11 = (k*Area[i-1])/(2*np.pi*(radi[i]+radi[i-1]**2))
        a12 =(
            np.pi*(psi[i]-psi[i-1])*(radi[i]+radi[i-1])/Area[i-1]
            +np.sin(psi[i-1])/radi[i-1] 
            -c0
        )
        a1 = a11*a12**2

        a2 = -k*(psi[i] - psi[i-1])/(radi[i] + radi[i-1])*a12

        a31 = k*Area[i]/(2*np.pi*(radi[i+1]+radi[i])**2)
        a32 =(
            np.pi*(psi[i+1]-psi[i])*(radi[i+1]+radi[i])/Area[i]
            + np.sin(psi[i])/radi[i]
            - c0
        )
        a3 = a31*a32**2

        a4 = -k*(psi[i+1]-psi[i])/(radi[i+1]+radi[i])*a32

        a5 = k*Area[i]*np.sin(psi[i])/(np.pi*(radi[i+1]+radi[i])*radi[i]**2)*a32

        a6 = -sigma*(
            Area[i]*radi[i+1] / (radi[i+1] + radi[i])**2
            -Area[i-1]*radi[i-1] / (radi[i] + radi[i-1])**2
        )
        
        a7 = kG*(psi[i+1]-psi[i])*np.sin(psi[i])/radi[i]**2

    if i == N-1:
        a11 = (k*Area[i-1])/(2*np.pi*(radi[i]+radi[i-1]**2))
        a12 =(
            np.pi*(psi[i]-psi[i-1])*(radi[i]+radi[i-1])/Area[i-1]
            +np.sin(psi[i-1])/radi[i-1] 
            -c0
        )**2
        a1 = a11*a12

        a21 = - k*(psi[i]-psi[i-1])/(radi[i] + radi[i-1])
        a22 = (
            np.pi*(psi[i]-psi[i-1])*(radi[i] + radi[i-1]) 
            +np.sin(psi[i-1])/radi[i-1] 
            -c0
        )
        a2 = a21*a22 
        

        a31 = (k*Area[i])/(2*np.pi*(radi[i+1]+radi[i]**2))
        a32 = (
            -np.pi*psi[i]*(radi[i+1]+radi[i])/Area[i]
            +np.sin(psi[i])/radi[i] 
            -c0
        )
        a3 = a31*a32**2

        a41 = k*psi[i]/(radi[i+1]+radi[i])
        a42 =a32
        a4 = a41*a42

        a51 = k*Area[i]*np.sin(psi[i])/(np.pi*radi[i]**2*(radi[i+1]+radi[i]))
        a52 = a32
        a5 = a51*a52

        a6 = -sigma*(
            Area[i]*radi[i+1]/(radi[i+1]+radi[i])**2
            -Area[i-1]*radi[i-1]/(radi[i]+radi[i-1])**2
        )

        a7 = - kG*psi[i]*np.sin(psi[i])/radi[i]**2

    a = a1 + a2 + a3 + a4 + a5 + a6 + a7
    
    return a/gamma(i)



def drdt_func(
        i,N,k,c0, sigma, kG, tau
        ,Area:list,psi:list,radi:list
        ,lamb:list , nu:list, z_list:list
        ):
    drdt,a1,a2,a3 = 0,0,0,0

    if i == 0:
        a11 = -2*np.pi*radi[i]*lamb[i]/Area[i]
        a12 = np.pi*nu[i]*(z_list[i+1] - z_list[i])/Area[i]

        a1 = a11+ a12

    if 0 < i < N-1:
        a21 = 2*np.pi*radi[i]*(
            lamb[i-1]/Area[i-1] - lamb[i]/Area[i] 
        )

        a22 = np.pi*(
            nu[i-1]*(z_list[i] - z_list[i-1])/Area[i-1]
            +nu[i]*(z_list[i+1] - z_list[i])/Area[i]
        )
        a2 = a21 + a22

    if i == N - 1:
        a31 = 2*np.pi*radi[i]*(
            lamb[i-1]/Area[i-1] - lamb[i]/Area[i] 
        )

        a32 = np.pi*(
            nu[i-1]*(z_list[i] - z_list[i-1])/Area[i-1] - nu[i]*z_list[i]/Area[i]
        )
        a3 = a31 + a32
    
    drdt = (a1 + a2 + a3) + Q_function(
                                        i=i,N=N,k=k,c0=c0
                                        ,sigma=sigma,kG=kG,tau=tau
                                        ,Area=Area,psi=psi,radi=radi
                                    )
    
    return drdt





if __name__ == "__main__":
    
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    psi_list, Area_list = const_args[12:14]
    
    radi_list = [L + r0 - i*ds for i in range(N,-1,-1)]

    Q = Q_function(
        i=N-1,N=N,k=k,c0=c0
        ,sigma=sigma,tau=tau,kG=kG
        ,Area=Area_list,psi=psi_list[0],radi=radi_list
    )

    drdt = drdt_func(
        i=N-1,N=N,k=k,c0=c0
        ,sigma=sigma,tau=tau,kG=kG
        ,Area=Area_list,psi=psi_list[0],radi=radi_list
        ,lamb=np.zeros(N), nu=np.zeros(N)
        ,z_list=np.zeros(N)
    )
    print(Q)
    print(drdt)
