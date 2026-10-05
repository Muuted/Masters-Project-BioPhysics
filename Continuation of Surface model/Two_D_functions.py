import numpy as np
import matplotlib.pyplot as plt
from two_d_data_processing import tot_area
#from Two_D_constants import gamma
#import os
#import pandas as pd
#import progressbar

def gamma(i,ds,eta):
    #eta = 1 # (micr gram)/(micro meter seconds)  #1e-3 Pa*s
    if isinstance(ds,float) or isinstance(ds,int) == True:
        gam = eta*ds# standard 1e0
        if i==0:
            return gam/2
        if i > 0:
            return gam
    elif isinstance(ds,list) == True:
        print("ds is a list, but I havent implemented the code yet to use that")
        exit()
    else:
        print("ds was neither integer, float or a list, check agian. Closing program")
        exit()


def Delta_s(A:list,r:list,i:int):
    return A[i]/(np.pi(r[i+1] +r[i]))

def Kronecker(i,j):
    if i == j:
        return 1
    if i != j:
        return 0

def B_function(
        i:int,N:int,c0:float,Area:list,psi:list,radi:list
        ):
    B = ""
    if 0 <= i < N-1:
        B = np.pi*(psi[i+1]-psi[i])*(radi[i+1]+radi[i])/Area[i] + np.sin(psi[i])/radi[i] - c0

    elif i == N-1:
        B = -np.pi*psi[i]*(radi[i+1] + radi[i])/Area[i] + np.sin(psi[i])/radi[i] - c0
    else:
        print("i never took a correct value")
        exit()
    if B == "":
        print(f"B function never took a value")
        exit()
    return B

def Q_function(
        i:int,N:int,k:float,c0:float, sigma:float, kG:float, tau:float
        ,Area:list,psi:list,radi:list
        ):
    
    a,a1,a2,a3,a4,a5,a6,a7 = 0,0,0,0,0,0,0,0
    if i == 0:
        B = B_function(i=i,N=N,c0=c0,Area=Area,psi=psi,radi=radi)
        
        a11 = -k*Area[i]/(2*np.pi)
        a12 = np.pi*(psi[i+1] - psi[i])/Area[i] - np.sin(psi[i])/(radi[i]**2)
        a = a11*B*a12 - tau

        return a 

    elif 0 < i < N-1:
        B_before = B_function(i=i-1,N=N,c0=c0,Area=Area,psi=psi,radi=radi)
        B = B_function(i=i,N=N,c0=c0,Area=Area,psi=psi,radi=radi)
        
        a11 = -k*(psi[i] - psi[i-1])/2
        a1 = a11*B_before 

        a21 = -k*Area[i]/(2*np.pi)
        a22 = np.pi*(psi[i+1] - psi[i])/Area[i] - np.sin(psi[i])/(radi[i]**2)
        a2 = a21*B*a22

        a = a1 + a2
        return a

    elif i == N-1:
        B_before = B_function(i=i-1 ,N=N,c0=c0,Area=Area,psi=psi,radi=radi)
        B = B_function(i=i,N=N,c0=c0,Area=Area,psi=psi,radi=radi)

        a11 = -k*(psi[i] - psi[i-1])/2
        a1 = a11*B_before

        a21 = k*Area[i]/(2*np.pi)
        a22 = np.pi*psi[i]/Area[i] + np.sin(psi[i])/(radi[i]**2)
        a2 = a21*B*a22

        a = a1 + a2
        return a

    else:
        print(f"\n index i={i} and is to larges as max N={N} this takes the max value of N-1={N-1}"
              +f"\n this is freom the Q_function")
        exit()

def drdt_func(
        i:int,N:int,k:float,c0:float, sigma:float, kG:float, tau:float, ds:float, eta:float
        ,Area:list
        ,psi:list,radi:list, z_list:list
        ,lamb:list , nu:list
        ):
    
    Q = Q_function(
            i=i,N=N,k=k,c0=c0
            ,sigma=sigma,kG=kG,tau=tau
            ,Area=Area,psi=psi,radi=radi
    )
    
    a1,a2,a3 = 0,0,0
    a11 ,a12 ,a21 ,a22 ,a31 ,a32 =0,0,0,0,0,0

    if i == 0:
        a11 = -2*np.pi*radi[i]*lamb[i]/Area[i]
        a12 = np.pi*nu[i]*(z_list[i+1] - z_list[i])/Area[i]

        return  (Q + a11 + a12)/gamma(i,ds=ds,eta=eta)

    elif 0 < i < N-1:
        a21 = 2*np.pi*radi[i]*(
            lamb[i-1]/Area[i-1] - lamb[i]/Area[i] 
        )

        a22 = np.pi*(
            nu[i-1]*(z_list[i] - z_list[i-1])/Area[i-1]
          + nu[i]*(z_list[i+1] - z_list[i])/Area[i]
        )
        return (Q + a21 + a22 )/gamma(i,ds=ds,eta=eta)

    elif i == N - 1:
        a31 = 2*np.pi*radi[i]*(
            lamb[i-1]/Area[i-1] - lamb[i]/Area[i] 
        )

        a32 = np.pi*(
            nu[i-1]*(z_list[i] - z_list[i-1])/Area[i-1] - nu[i]*z_list[i]/Area[i]
        )
        return (Q + a31 + a32)/gamma(i,ds=ds,eta=eta)

    else:
        print(f"\n in drdt_function, non of the available i, as was i={i} where max is N-1={N-1}")
        exit()
    #drdt = (a1 + a2 + a3) + Q
                                    
    
    #return drdt/gamma(i,ds=ds,eta=eta)



def dzdt_func(
        i:int,ds:float,eta:float,Area:list,radi:list, nu:list
        ):
    
    if i == 0 :
        dzdt = - np.pi*nu[i]*(radi[i+1] + radi[i])/(Area[i])
        return dzdt/gamma(i,ds=ds,eta=eta)
    
    elif 0 < i  :
        dzdt = np.pi*(
            nu[i-1]*(radi[i] + radi[i-1])/Area[i-1]
          - nu[i]*(radi[i+1] + radi[i])/Area[i]
        )
        return dzdt/gamma(i,ds=ds,eta=eta)
    
    else:
        print(f"\n in dzdt_func i was not i either 0 or greater than 0,the value of i={i} \n error program is terminated")
        exit()
    #return dzdt/gamma(i,ds=ds,eta=eta)



def dpsidt_func(  i:int,N:int,k:float,c0:float, sigma:float, kG:float, tau:float, ds:float ,eta:float
        ,Area:list,psi:list,radi:list
        ,lamb:list , nu:list, z_list:list
        ):
    a1 ,a2 ,a3 = 0,0,0
    if 0 <= i < N-1 :
        dzdt_i_next = dzdt_func(i=i+1,ds=ds,eta=eta,Area=Area,radi=radi,nu=nu)
        dzdt_i = dzdt_func(i=i,ds=ds,eta=eta,Area=Area,radi=radi,nu=nu)

        drdt_i_next = drdt_func(
            i=i+1
            ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau, ds=ds,eta=eta
            ,Area=Area,psi=psi,radi=radi,z_list=z_list
            ,lamb=lamb,nu=nu
            )
        
        drdt_i = drdt_func(
            i=i
            ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau, ds=ds,eta=eta
            ,Area=Area,psi=psi,radi=radi,z_list=z_list
            ,lamb=lamb,nu=nu
            )
        
        a1 = (radi[i+1] + radi[i])*np.cos(psi[i])
        a2 = (z_list[i+1] - z_list[i])*np.cos(psi[i]) - 2*radi[i+1]*np.sin(psi[i])
        a3 = (z_list[i+1] - z_list[i])*np.cos(psi[i]) + 2*radi[i]*np.sin(psi[i])

        dpsidt = np.pi*(  a1*(dzdt_i_next - dzdt_i) + a2*drdt_i_next + a3*drdt_i  )/Area[i]
    
        return dpsidt

    elif i == N-1:
        dzdt_i = dzdt_func(i=i,ds=ds,eta=eta,Area=Area,radi=radi,nu=nu)
        drdt_i = drdt_func(
            i=i
            ,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau, ds=ds,eta=eta
            ,Area=Area,psi=psi,radi=radi,z_list=z_list
            ,lamb=lamb,nu=nu
            )
        
        a1 = -(radi[i+1] + radi[i])*np.cos(psi[i])
        a3 = -z_list[i]*np.cos(psi[i]) + 2*radi[i]*np.sin(psi[i])

        dpsidt = np.pi*( a1*dzdt_i  + a3*drdt_i )/Area[i]

        return dpsidt
    
    else:
        print(f"\n In dpsidt_func. Non of the values for i was chosen, value of i={i}  \n but has to be in range 0 to N-1={N-1}. Program termined")
        exit()



def dSdpsi_func(i:int,N:int,c0:float,k:float,kG:float,r:list,psi:list,Area:list)->float:
    return_val = ""
    if i == 0:
        a11 = k*Area[i]/(2*np.pi)
        a12 = B_function(i=i,N=N,c0=c0,Area=Area,psi=psi,radi=r)
        a13 = (
            -np.pi*(r[i+1]+r[i])/Area[i] + np.cos(psi[i])/r[i]
            )
        a21 = kG*(
            -np.sin(psi[i]) + (psi[i+1]-psi[i])*np.cos(psi[i])
        )
        return a11*a12*a13 + a21
    elif 0 < i < N-1 :
        a11 = k*(r[i]+r[i-1])/2
        a12 = B_function(i=i-1 ,N=N,c0=c0,Area=Area,psi=psi,radi=r)
        
        a21 = k*Area[i]/(2*np.pi)
        a22 = B_function(i=i,N=N,c0=c0,Area=Area,psi=psi,radi=r)
        a23 = (
            -np.pi*(r[i+1]+r[i])/Area[i] + np.cos(psi[i])/r[i]
            )
        a31 = kG*(
            np.sin(psi[i-1]) - np.sin(psi[i]) + (psi[i+1]-psi[i])*np.cos(psi[i])
        )
        return a11*a12 + a21*a22*a23 + a31        
    elif i == N - 1 :
        a11 = k*(r[i]+r[i-1])/2
        a12 = B_function(i=i-1,N=N,c0=c0,Area=Area,psi=psi,radi=r)
        
        a21 = k*Area[i]/(2*np.pi)
        a22 = B_function(i=i,N=N,c0=c0,Area=Area,psi=psi,radi=r)
        a23 = (
            -np.pi*(r[i+1]+r[i])/Area[i] + np.cos(psi[i])/r[i]
            )
        a31 = kG*(
            np.sin(psi[i-1]) - np.sin(psi[i]) - psi[i]*np.cos(psi[i])
        )
        return a11*a12 + a21*a22*a23 + a31        
    else:
        print(f"\n i never took a value from the possible set. 0<= i={i} <= N-1 ")
    if return_val == "":
        print("\n return val never took value \n")
        exit()

    #return return_val


def constraint_f(i:int,N:int,r:list,psi:list,Area:list) -> float:
    f = ""
    if 0 <= i <= N-1:
        f = np.pi*(r[i+1]**2 - r[i]**2)/Area[i] - np.cos(psi[i])
    else:#elif i > N - 1:
        print(f"the value of i is to large, in the constraint equation \n"
              +f"i={i} and N={N}")
        exit()
    if f == "":
        print(f"value of f never took a number")
        exit()
    return f

def constraint_g(i:int,N:int,r:list,z:list,psi:list,Area:list)-> float:
    g = ""
    if 0 <= i < N-1:
        g = np.pi*(z[i+1]-z[i])*(r[i+1] + r[i])/Area[i] - np.sin(psi[i])
    elif i == N-1:
        g = -np.pi*z[i]*(r[i+1] + r[i])/Area[i] - np.sin(psi[i])
    else:#elif i >= N:
        print(f"the value of i is to large, in the constraint equation")
        exit()
    if g == "":
        print(f"value of f new took a number")
        exit()
    return g


def check_constraints_truth(N:int,r:list,z:list,psi:list,Area:list,tol:float)->bool:
    err = False
    for i in range(N):
        f = constraint_f(i=i,N=N,r=r,psi=psi,Area=Area)
        g = constraint_g(i=i,N=N,r=r,z=z,psi=psi,Area=Area)

        if tol < abs(f) or tol < abs(g):
            err = True
            break

    return err




if __name__ == "__main__":
    pass

    