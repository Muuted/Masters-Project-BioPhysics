import numpy as np
from Two_D_constants import Two_D_Constants, gamma
#import os
#import pandas as pd
#import progressbar


def Delta_s(A:list,r:list,i:int):
    return A[i]/(np.pi(r[i+1] +r[i]))

def Kronecker(i,j):
    if i == j:
        return 1
    if i != j:
        return 0

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
        a1 = a11*(a12**2)

        a2 = - k*(psi[i+1]-psi[i])/(radi[i+1] +radi[i])*a12

        a31 = k*Area[i]*np.sin(psi[i])/(np.pi*(radi[i+1]+radi[i])*radi[i]**2)
        a3 = a31*a12 - tau

        a4 = - sigma*Area[i]*radi[i+1]/(np.pi*(radi[i+1]+radi[i])**2)

        a5 = 0#kG*(psi[i+1]-psi[i])*np.sin(psi[i])/radi[i]**2


    if 0 < i < N-1:
        
        a11 = (k*Area[i-1])/(2*np.pi*(radi[i]+radi[i-1]**2))
        a12 =(
            np.pi*(psi[i]-psi[i-1])*(radi[i]+radi[i-1])/Area[i-1]
            +np.sin(psi[i-1])/radi[i-1] 
            -c0
        )
        a1 = a11*(a12**2)

        a2 = -k*(psi[i] - psi[i-1])/(radi[i] + radi[i-1])*a12

        a31 = k*Area[i]/(2*np.pi*(radi[i+1]+radi[i])**2)
        a32 =(
            np.pi*(psi[i+1]-psi[i])*(radi[i+1]+radi[i])/Area[i]
            + np.sin(psi[i])/radi[i]
            - c0
        )
        a3 = a31*(a32**2)

        a4 = -k*(psi[i+1]-psi[i])/(radi[i+1]+radi[i])*a32

        a5 = k*Area[i]*np.sin(psi[i])/(np.pi*(radi[i+1]+radi[i])*radi[i]**2)*a32

        a6 = -sigma*(
            Area[i]*radi[i+1] / (radi[i+1] + radi[i])**2
            -Area[i-1]*radi[i-1] / (radi[i] + radi[i-1])**2
        )
        
        a7 = 0#kG*(psi[i+1]-psi[i])*np.sin(psi[i])/radi[i]**2

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

        a7 = 0#- kG*psi[i]*np.sin(psi[i])/radi[i]**2

    a = a1 + a2 + a3 + a4 + a5 + a6 + a7
    
    return a/gamma(i)



def drdt_func(
        i,N,k,c0, sigma, kG, tau
        ,Area:list,psi:list,radi:list
        ,lamb:list , nu:list, z_list:list
        ):
    a1,a2,a3 = 0,0,0

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


def drdt_RungeKutta_4(
        i,dt,N,k,c0, sigma, kG, tau
        ,Area:list,psi:list,radi:list
        ,lamb:list , nu:list, z_list:list
        ):
    
    pass

def dzdt_func(
        i,Area:list,radi:list, nu:list
        ):
    
    if i == 0 :
        dzdt = - np.pi*nu[i]*(radi[i+1] + radi[i])/(gamma(i)*Area[i])
    if i > 0 :
        dzdt = np.pi*(
            nu[i-1]*(radi[i] + radi[i-1])/Area[i-1]
            -nu[i]*(radi[i+1] + radi[i])/Area[i]
        )/gamma(i)
    
    return dzdt

def  dzdt_RungeKutta_4(i,dt,Area:list,radi:list, nu:list):
    
    dzdt_1 = dzdt_func(i=i,Area=Area,radi=radi,nu=nus)

    dzdt_2 = dzdt_func(i=i,Area=Area,radi=radi,nu=nus)

def dpsidt_func(  i,N,k,c0, sigma, kG, tau
        ,Area:list,psi:list,radi:list
        ,lamb:list , nu:list, z_list:list
        ):
    
    if i < N-1 :
        dzdt_i_next = dzdt_func(i=i+1,Area=Area,radi=radi,nu=nu)
        dzdt_i = dzdt_func(i=i,Area=Area,radi=radi,nu=nu)

        drdt_i_next = drdt_func(
            i=i+1
            ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
            ,Area=Area,psi=psi,radi=radi,z_list=z_list
            ,lamb=lamb,nu=nu
            )
        
        drdt_i = drdt_func(
            i=i
            ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
            ,Area=Area,psi=psi,radi=radi,z_list=z_list
            ,lamb=lamb,nu=nu
            )
        
        a1 = (radi[i+1] + radi[i])*np.cos(psi[i])
        a2 = (z_list[i+1]-z_list[i])*np.cos(psi[i]) - 2*radi[i+1]*np.sin(psi[i])
        a3 =(z_list[i+1]-z_list[i])*np.cos(psi[i]) - 2*radi[i]*np.sin(psi[i])

        dpsidt = np.pi*(   a1*(dzdt_i_next - dzdt_i) + a2*drdt_i_next + a3*drdt_i  )/Area[i]

    if i == N-1:
        dzdt_i = dzdt_func(i=i,Area=Area,radi=radi,nu=nu)
        drdt_i = drdt_func(
            i=i
            ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
            ,Area=Area,psi=psi,radi=radi,z_list=z_list
            ,lamb=lamb,nu=nu
            )
        
        a1 = (radi[i+1] + radi[i])*np.cos(psi[i])
        a3 =(z_list[i+1]-z_list[i])*np.cos(psi[i]) - 2*radi[i]*np.sin(psi[i])

        dpsidt = np.pi*( -a1*dzdt_i  + a3*drdt_i )/Area[i]

    return dpsidt

def dPsidt_RungeKutta_4(
        i,N,k,c0,sigma,kG,tau
        ,Area,radi,z_list
        ,lambs,nus
        ,psi
        ):
    
    dpdt_1 = dpsidt_func(
        i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi=psi 
                )

    dpdt_2 = dpsidt_func(i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi= psi + (dt/2)*dpdt_1
    )

    dpdt_3 = dpsidt_func(i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi= psi + (dt/2)*dpdt_2
    )

    dpdt_4 = dpsidt_func(i=i
        ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area,radi=radi,z_list=z_list
        ,lamb=lambs,nu=nus
        ,psi= psi + dt*dpdt_3
    )

    Runge_kutta = dpdt_1 + 2*dpdt_2 + 2*dpdt_3 + dpdt_4
    return Runge_kutta


def Langrange_multi(
         N ,k ,c0 ,sigma ,kG ,tau #, ds
        ,Area:list,psi:list
        ,radi:list
        ,z_list:list
        ,print_matrix = False
        ):

    A = np.zeros(shape=(2*N,2*N))
    b = np.zeros(2*N)

    for i in range(2*N):
        for j in range(2*N):
            l = i%N 
            n = i%N + N
            nu_i_next ,nu_i ,nu_i_before = 0,0,0
            lamb_i_next ,lamb_i ,lamb_i_before = 0,0,0
            lambs, nus = 0,0
            # The first N equations
            if i == 0:
                #---------------- calc lambda vals -------------------------
                if i+1==j:
                    lamb_i_next = (
                        -2*np.pi**2*(
                            radi[l+1]*(z_list[l+1] - z_list[l])*np.sin(psi[l])
                            +2*radi[l+1]**2*np.cos(psi[l])
                        )/(gamma(l+1)*Area[l+1]*Area[l])
                    )#*Kronecker(i+1,j)
                
                if i == j:
                    lamb_i = (
                        -2*np.pi**2*(
                            (z_list[l+1] - z_list[l] )*radi[l]*np.sin(psi[l]) 
                            + 2*radi[l]**2*np.cos(psi[l])
                        )/(gamma(l)*Area[l]**2)
                    )#*Kronecker(i,j)

                #lambs = lamb_i_next + lamb_i

                #---------------- calc nu vals -------------------------
                if n+1== j:
                    nu_i_next_1 = (
                        (   (z_list[l+1] - z_list[l] )*np.sin(psi[l])  
                            +2*radi[l+1]*np.cos(psi[l]) 
                        )*(z_list[l+2]-z_list[l+1])
                    )#/(gamma(l+1)*Area[l+1]*Area[l])

                    nu_i_next_2 = -(
                        (radi[l+1] + radi[l])*(radi[l+2] + radi[l+1])*np.sin(psi[l])
                    )#/(gamma(l+1)*Area[l+1]*Area[l])

                    nu_i_next = np.pi**2*(nu_i_next_1 + nu_i_next_2)/(gamma(l+1)*Area[l+1]*Area[l])#*Kronecker(n+1,j)

                if n == j :
                    nu_i_1 = (radi[l+1] + radi[l])**2*np.sin(psi[l])

                    nu_i_2 =(
                        (z_list[l+1]-z_list[l])*np.sin(psi[l]) + 2*radi[l]*np.cos(psi[l])
                    )*(z_list[l+1]-z_list[l])

                    nu_i = np.pi**2*(nu_i_1 + nu_i_2)/(gamma(l)*Area[l]**2) # *Kronecker(n,j)

                #nus = nu_i_next + nu_i
            
            if 0 < i < N-2:
                if i+1==j :
                    #print(f"l={l} and j={j} and i ={i}")
                    lamb_i_next = -(
                        (2*np.pi**2*radi[l]/(gamma(l+1)*Area[l+1]*Area[l]))*(
                            (z_list[l+1]-z_list[l])*np.sin(psi[l]) +2*radi[l+1]*np.cos(psi[l])
                        )
                    )#*Kronecker(i+1,j)

                if i == j:
                    lamb_i = (2*np.pi**2/Area[l]**2)*(
                        (z_list[l+1]-z_list[l])*(radi[l+1]-radi[l])*np.sin(psi[l])
                        +2*(radi[l+1]**2 - radi[l]**2) * np.cos(psi[l])
                    )#*Kronecker(i,j)
                    lamb_i = 2*np.pi**2*(
                        (z_list[l+1] - z_list[l])*np.sin(psi[l])*( radi[l+1]/gamma(l+1) - radi[l]/gamma(l))
                        + 2*np.cos(psi[l])*( radi[l+1]**2/gamma(l+1) - radi[l]**2/gamma(l))
                    )/Area[l]**2

                if i-1 ==j:
                    lamb_i_before = 2*np.pi**2*radi[l]*(
                        (z_list[l+1]-z_list[l])*np.sin(psi[l]) + 2*radi[l]*np.cos(psi[l])
                    )/(gamma(l)*Area[l]*Area[l-1]) #*Kronecker(i-1,j)

                #lambs = lamb_i_next + lamb_i + lamb_i_before
                
                #---------- nu ----------------------------------
                if n+1 == j:
                    nu_i_next = np.pi**2*(
                        2*radi[l+1]*(z_list[l+2] - z_list[l+1])*np.cos(psi[l])
                        +(z_list[l+1] - z_list[l] - (radi[l+1] + radi[l])*(radi[l+2] + radi[l+1]))*np.sin(psi[l])
                    )/(gamma(l+1)*Area[l+1]*Area[l]) # *Kronecker(n+1,j)
                
                if n == j:
                    nu_i = (np.pi**2/Area[l]**2)*(
                        (1/gamma(l+1) + 1/gamma(l))*( (radi[l+1]+ radi[l])**2 + (z_list[l+1] - z_list[l])**2 )*np.sin(psi[l])
                        + 2*(z_list[l+1] - z_list[l])*( radi[l+1]/gamma(l+1) + radi[l]/gamma(l)  )*np.cos(psi[l])
                    )#*Kronecker(n,j)

                if n-1 == j:
                    nu_i_before = np.pi**2/(gamma(l)*Area[l]*Area[l-1])*(
                        (
                            (z_list[l+1] - z_list[l])*np.sin(psi[l]) + 2*radi[l]*np.cos(psi[l])
                        )*(z_list[l] - z_list[l-1])
                        -(radi[l+1] + radi[l])**2*np.sin(psi[l])
                    )#*Kronecker(n-1,j)

                #nus = nu_i_next + nu_i + nu_i_before
            
            if i == N-2:#### do this ########## do this ############ do this ######################################################          
                if i+1 == j:
                    lamb_i_next = -2*np.pi**2*(
                                (z_list[l+1] - z_list[l])*radi[l+1]*np.sin(psi[l])
                                + 2*radi[l+1]**2*np.cos(psi[l])
                    )/(gamma(l+1)*Area[l+1]*Area[l])
                
                if i == j:
                    lamb_i = 2*np.pi**2*(
                        (radi[l+1]/gamma(l+1) -  radi[l]/gamma(l))*(z_list[l+1] - z_list[l])*np.sin(psi[l])
                        + 2*( radi[l+1]**2/gamma(l+1) -  radi[l]**2/gamma(l))*np.cos(psi[l])
                    )/Area[l]**2

                if i-1 == j:
                    lamb_i_before = 2*np.pi**2*(
                        (z_list[l+1] - z_list[l])*radi[l]*np.sin(psi[l]) + 2*radi[l]**2*np.cos(psi[l])
                    )/(gamma(l)*Area[l]*Area[l-1])

                #lambs = lamb_i_next + lamb_i + lamb_i_before

                if n+1 == j:
                    nu_i_next_1 =(radi[l+1] + radi[l])*(radi[l+2] + radi[l+1])*np.sin(psi[l])
                    
                    nu_i_next_2 = ((z_list[l+1]- z_list[l])*np.sin(psi[l]) + 2*radi[l+1]*np.cos(psi[l]))*z_list[l+1]

                    nu_i_next = -np.pi**2*( nu_i_next_1 + nu_i_next_2)/(gamma(l+1)*Area[l+1]*Area[l])

                if n == j:
                    nu_i_1 = (
                        ((radi[l+1] + radi[l])**2 +(z_list[l+1]- z_list[l])**2)*np.sin(psi[l])
                        )*(1/gamma(l+1) +1/gamma(l))
                    
                    nu_i_2 = 2*(z_list[l+1]- z_list[l])*np.cos(psi[l])*(radi[l+1]/gamma(l+1) + radi[l]/gamma(l))

                    nu_i = np.pi**2*(nu_i_1 + nu_i_2 )/Area[l]**2

                if n-1 == j:

                    nu_i_before_1 = (
                        (z_list[l+1]-z_list[l])*np.sin(psi[l]) + 2*radi[l]*np.cos(psi[l]) 
                    )*(z_list[l]-z_list[l-1])

                    nu_i_before_2 = -(radi[l+1] + radi[l])*(radi[l] + radi[l-1])*np.sin(psi[l])
                    
                    nu_i_before = np.pi**2*(nu_i_before_1 + nu_i_before_2)/(gamma(l)*Area[l]*Area[l-1])

                #nus = nu_i_next + nu_i + nu_i_before

            if i == N-1:     
                if i == j:
                    lamb_i = -2*np.pi**2*radi[l]*(
                        2*radi[l]*np.cos(psi[l]) - z_list[l]*np.sin(psi[l])
                    )/( gamma(l)*Area[l]**2  )#*Kronecker(i,j)

                if i-1== j:
                    lamb_i_before = 2*np.pi**2*radi[l]*(
                            2*radi[l]*np.cos(psi[l]) - z_list[l]*np.sin(psi[l])
                        )/(gamma(l)*Area[l]*Area[l-1])  #*Kronecker(i-1,j)
                
                #lambs = lamb_i_next + lamb_i + lamb_i_before
                #---------- nu ----------------------------------
                if n == j:
                    nu_i = (np.pi**2/(gamma(l)*Area[l]**2))*(
                        2*radi[l]*z_list[l]*np.cos(psi[l])
                       - ((radi[l+1] + radi[l])**2 + z_list[l]**2)*np.sin(psi[l])
                    )#*Kronecker(n,j)
                if n-1 == j:
                    nu_i_before = (np.pi**2/gamma(l))*(
                        (2*radi[l]*np.cos(psi[l]) - z_list[l]*np.sin(psi[l]))*(z_list[l] - z_list[l-1])
                        -(radi[l+1]+radi[l])*(radi[l]+radi[l-1])*np.sin(psi[l])
                    )/(Area[l]*Area[l-1]) # *Kronecker(n-1,j)

                #nus = nu_i_next + nu_i + nu_i_before
            # The 2nd N equations. where l is the shifted index

            if i > N - 1:
                if i%N == j:
                    lamb_i = np.sin(psi[l])*Kronecker(i%N,j)
                #lambs = lamb_i_next + lamb_i + lamb_i_before
                if n==j:
                    nu_i = -np.cos(psi[l])*Kronecker(n,j)
                #nus = nu_i_next + nu_i + nu_i_before

            if print_matrix == True:
                A[i][j] =  '{:.0e}'.format(lamb_i_next + lamb_i + lamb_i_before + nu_i_next + nu_i + nu_i_before) #'{:.0e}'.format(lambs+nus)
            else:
                A[i][j] = lamb_i_next + lamb_i + lamb_i_before + nu_i_next + nu_i + nu_i_before #lambs + nus
            

            #---------------- calc b vals -------------------------     
        if i < N :   
            Q_i_next = Q_function(
                i=i+1
                ,N=N,k=k,c0=c0,sigma=sigma,kG=kG
                ,tau=tau,Area=Area,psi=psi,radi=radi
            )
            Q_i = Q_function(
                i=i
                ,N=N,k=k,c0=c0,sigma=sigma,kG=kG
                ,tau=tau,Area=Area,psi=psi,radi=radi
            )

            b1 = (
                -Q_i_next*(
                    np.pi*(
                        ( z_list[i+1]-z_list[i])*np.sin(psi[i] ) 
                        +2*radi[i+1]*np.cos(psi[i])
                        )/(gamma(i+1)*Area[i])
                )
                -Q_i*(
                    np.pi*(
                        ( z_list[i+1]-z_list[i])*np.sin(psi[i] ) 
                        +2*radi[i]*np.cos(psi[i])
                        )/(gamma(i)*Area[i])
                )
            )
        if i > N - 1:
            if i == N :
                b11 = -k*(
                    np.pi*(psi[i%N+1] - psi[i%N])*(radi[i%N+1]+radi[i%N])/Area[i%N]
                    +np.sin(psi[i%N])/radi[i%N] - c0
                )*(
                    1 - Area[i%N]*np.cos(psi[i%N])/(np.pi*radi[i%N]*(radi[i%N+1] + radi[i%N]))
                )
                b12 = -kG*(  np.sin(psi[i%N]) - ( psi[i%N+1] - psi[i%N] )*np.cos(psi[i%N])  )
                
                b1 = b11 + b12
            if N < i < 2*N -1:
                b21 = -k*(
                    np.pi*(psi[i%N+1] - psi[i%N])*(radi[i%N+1]+radi[i%N])/Area[i%N]
                    +np.sin(psi[i%N])/radi[i%N] - c0
                )*(
                    1 - Area[i%N]*np.cos(psi[i%N])/( np.pi*radi[i%N]*( radi[i%N+1] + radi[i%N] ) )
                )
            
                b22 =  k*(
                    np.pi*( psi[i%N] - psi[i%N-1] )*( radi[i%N] + radi[i%N-1] )/Area[i%N-1]
                    +np.sin(psi[i%N-1])/radi[i%N-1] - c0
                )
                
                b23 = -kG*(
                    np.sin(psi[i%N]) - np.sin(psi[i%N-1])
                    - (psi[i%N+1] - psi[i%N] )*np.cos(psi[i%N])
                )

                b1 = b21 + b22 + b23
            
            if i == 2*N -1:
                b21 = -k*(
                    -np.pi*psi[i%N]*(radi[i%N+1]+radi[i%N])/Area[i%N]
                    + np.sin(psi[i%N])/radi[i%N] - c0
                )*(
                    1 - Area[i%N]*np.cos(psi[i%N])/(np.pi*radi[i%N]*(radi[i%N+1] + radi[i%N]))
                )
            
                b22 =  k*(
                    np.pi*(psi[i%N] - psi[i%N-1])*(radi[i%N]+radi[i%N-1])/Area[i%N-1]
                    + np.sin(psi[i%N-1])/radi[i%N-1] - c0
                )
                
                b23 = -kG*(
                    np.sin(psi[i%N]) - np.sin(psi[i%N-1])
                    + psi[i%N]*np.cos(psi[i%N])
                )

                b1 = b21 + b22 + b23
        
        if print_matrix == True:
            b[i] = '{:.0e}'.format(b1)
        else:
            b[i] = b1

    if print_matrix == True:
        print(f"A: {np.shape(A)[0]}x{np.shape(A)[1]}\n ",A)
        print("b:",b)
        #print("x:",x)
        print("psi:",psi)
        x = np.linalg.solve(A,b)
        lamb_return = x[0:N]
        nu_return = x[N:2*N]
        print("len(lamb_return):",np.shape(lamb_return))
        print("len(nus_return):",np.shape(nu_return))
    else:
        x = np.linalg.solve(A,b)

        lamb_return = x[0:N]
        nu_return = x[N:2*N]

    return lamb_return,nu_return


def constraint_f(i:int,N:int,r:list,psi:list,Area:list):
    if i > N:
        print(f"the value of i is to large, in the constraint equation")
        exit()
    else:
        f = np.pi*(r[i+1]**2 - r[i]**2)/Area[i] - np.cos(psi[i])

    return f

def constraint_g(i:int,N:int,r:list,z:list,psi:list,Area:list):
    if i > N:
        print(f"the value of i is to large, in the constraint equation")
        exit()
    else:
        g = np.pi*(z[i+1]-z[i])*(r[i+1] + r[i])/Area[i] - np.sin(psi[i])
    
    return g


def c_diff_f(
        i:int,j:int,N:int
        ,r:list,psi:list,Area:list
        ,diff_var:str =""
        ):
    diff_var_list = ["r","z","psi"]
    df = 0
    if diff_var == "r":
        df = 2*np.pi*(r[i+1]*Kronecker(i+1,j) - r[i]*Kronecker(i,j))/Area[i]
        
    if diff_var == "z":
        df = 0

    if diff_var == "psi":
        df = np.sin(psi[i])*Kronecker(i,j)
        
    if diff_var not in diff_var_list:
            #Error handling
            print(f"wrong diff variable of either (r,z,psi) error in c_diff_f")
            exit()

    return df


def c_diff_g(
        i:int,j:int,N:int
        ,r:list,z:list,psi:list,Area:list
        ,diff_var:str =""
        ):
    diff_var_list = ["r","z","psi"]
    dg = 0
    if diff_var == "r":
        dg = np.pi*(z[i+1]-z[i])*(  Kronecker(i+1,j) +Kronecker(i,j)    )/Area[i]
        
    if diff_var == "z":
        dg = np.pi*(r[i+1]+r[i])*(Kronecker(i+1,j) -Kronecker(i,j))/Area[i]

    if diff_var == "psi":
        dg = -np.cos(psi[i])*Kronecker(i,j)

    if diff_var not in diff_var_list:
            #Error handling
            print(f"wrong diff variable of either (r,z,psi) error in c_diff_g")
            exit()

    return dg


def Epsilon_values(
        N:int,r:list
        ,z:list,psi:list,Area:list
        ,print_matrix = False
        ):
    
    A = np.zeros(shape=(2*N,2*N))
    b = np.zeros(2*N)

    for alpha in range(2*N):
        for beta in range(2*N):
            l = alpha%N 
            n = beta%N # + N
            K = 0
            if alpha < N and beta < N:
                K = 0
                for j in range(N):
                    K += (
                        c_diff_f(i=l,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="r")*c_diff_f(i=n,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="r")
                        +c_diff_f(i=l,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="z")*c_diff_f(i=n,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="z")
                        +c_diff_f(i=l,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="psi")*c_diff_f(i=n,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="psi")
                    )
                   
            if alpha < N and beta >= N:
                K = 0
                for j in range(N):
                    K += (
                        c_diff_f(i=l,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="r")*c_diff_g(i=n,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="r")
                        +c_diff_f(i=l,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="z")*c_diff_g(i=n,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="z")
                        +c_diff_f(i=l,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="psi")*c_diff_g(i=n,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="psi")
                    )

            if alpha >= N  and beta < N:
                K = 0
                for j in range(N):
                    K += (
                        c_diff_g(i=l,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="r")*c_diff_f(i=n,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="r")
                        +c_diff_g(i=l,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="z")*c_diff_f(i=n,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="z")
                        +c_diff_g(i=l,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="psi")*c_diff_f(i=n,j=j,N=N,r=r,psi=psi,Area=Area,diff_var="psi")
                    )

            if alpha >= N and beta >= N:
                K = 0
                for j in range(N):
                    K += (
                        c_diff_g(i=l,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="r")*c_diff_g(i=n,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="r")
                        + c_diff_g(i=l,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="z")*c_diff_g(i=n,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="z")
                        + c_diff_g(i=l,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="psi")*c_diff_g(i=n,j=j,N=N,r=r,psi=psi,z=z,Area=Area,diff_var="psi")
                    )
            
            A[alpha][beta] = K
            
        if alpha < N:
            b[alpha] = -constraint_f(i=alpha%N,N=N,r=r,psi=psi,Area=Area)
        
        if alpha >= N:
            b[alpha] = -constraint_g(i=alpha%N,N=N,r=r,z=z,psi=psi,Area=Area)

    if print_matrix == True:
        print(f"A: {np.shape(A)[0]}x{np.shape(A)[1]}\n ",A)
        print("b:",b)
        x = np.linalg.solve(A,b)
        epsilon_f = x[0:N]
        epsilon_g = x[N:2*N]
    else:
        x = np.linalg.solve(A,b)
        epsilon_f = x[0:N]
        epsilon_g = x[N:2*N]

    return epsilon_f, epsilon_g


if __name__ == "__main__":
    a = np.zeros(5)
    print(a)
    exit()
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    
    radi_list = [L + r0 - i*ds for i in range(N,-1,-1)]
    zz_list = [0 for i in range(N,-1,-1)]

    print(
        f"radi shape = {np.shape(radi_list)} \n"
        +f"z shape = {np.shape(zz_list)} \n"
        +f"psi shape = {np.shape(psi_list)} \n"
        f"Area shape = {np.shape(Area_list)} \n"
    )
    #exit()
    lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area_list
        ,psi=psi_list[0]
        ,radi=radi_list
        ,z_list=zz_list
        ,print_matrix=True
    )

    