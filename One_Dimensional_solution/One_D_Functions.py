import numpy as np
from One_D_Constants import One_D_Constants,gamma
#import os
#import pandas as pd
#import progressbar

    
def Kroncker(i,j):
    if i == j :
        return 1
    else:
        return 0

def dPsidt(i,N,multi,psi,deltaS):
    if int(len(multi))%2 != 0:
        print("wrong dimensions on lagrange multiplier list")
        exit()
    l = i + int(len(multi)/2)
    
    if i == 0:
        a1 = (multi[i+1]/gamma(i+1) - multi[i]/gamma(i))*np.sin(psi)/deltaS**2
        a2 = (multi[l]/gamma(i) - multi[l+1]/gamma(i+1) )*np.cos(psi)/deltaS**2

    if 0 < i < N - 1:
        a11 = (multi[i-1] - multi[i])/gamma(i)
        a12 = (multi[i] - multi[i+1])/gamma(i+1) 
        a1 = (a11 - a12)*np.sin(psi)/deltaS**2

        a21 = (multi[l] - multi[l+1])/gamma(i+1) 
        a22 = (multi[l-1] - multi[l])/gamma(i)
        a2 = (a21 - a22)*np.cos(psi)/deltaS**2

    if i == N - 1:
        a1 = ( multi[i-1] - multi[i])*np.sin(psi)/(gamma(i)*deltaS**2)
        a2 = -(multi[l-1]- multi[l])*np.cos(psi)/(gamma(i)*deltaS**2)

    # updating the new psi.
    return a1 + a2
    #psi[t+1][i] = psi[t][i] + dt*(a1 +a2)
    
def dPsidt_RungeKutta_4(link,N,ds,dt,multipliers,psi):
    dpdt_1 = dPsidt(
        i=link, N=N, deltaS=ds
        ,multi=multipliers
        ,psi= psi
            )

    dpdt_2 = dPsidt(
        i=link, N=N, deltaS=ds
        ,multi=multipliers
        ,psi= psi + (dt/2)*dpdt_1
    )

    dpdt_3 = dPsidt(
        i=link, N=N, deltaS=ds
        ,multi=multipliers
        ,psi= psi + (dt/2)*dpdt_2
    )

    dpdt_4 = dPsidt(
        i=link, N=N, deltaS=ds
        ,multi=multipliers
        ,psi= psi + dt*dpdt_3
    )

    Runge_kutta = dpdt_1 + 2*dpdt_2 + 2*dpdt_3 + dpdt_4
    return Runge_kutta


def Lagran_multi(
        psi_list,t,k,c0,ds
        ,num_chains
        ,linalg_lstsq =True
        ,print_matrix = False
                 ):
    
    N = num_chains
    b = np.full(shape=(2*N),fill_value=0,dtype=float)
    A = np.full(shape=(2*N,2*N),fill_value=10,dtype=float)

    for i in range(2*N):
        a,b1 = 0,0
        for j in range(2*N):
            l = i%N + N
            if i == 0:
                b1 = 0

                a1 = Kroncker(i,j)/gamma(i%N) - Kroncker(i+1,j)/gamma((i+1)%N)

                a2 = Kroncker(l,j)/gamma(i%N) - Kroncker(l+1,j)/gamma((i+1)%N)

                a = a1*np.cos(psi_list[t][j%N]) + a2*np.sin(psi_list[t][j%N])

            if 0 < i < N - 1 :
                b1 = 0
                
                a11 = (Kroncker(i,j) - Kroncker(i+1,j))/gamma((i+1)%N)
                a12 = (Kroncker(i-1,j) - Kroncker(i,j))/gamma(i%N)
                a1 = (a11 - a12)*np.cos(psi_list[t][j%N])
                
                a21 = (Kroncker(l,j) - Kroncker(l+1,j))/gamma((i+1)%N)
                a22 = (Kroncker(l-1,j) - Kroncker(l,j))/gamma(i%N)
                a2 = (a21 - a22)*np.sin(psi_list[t][j%N])
                
                a = a1 + a2
            if i == N - 1:
                b1 = 0
                a1 = (Kroncker(i-1,j) - Kroncker(i,j))*np.cos(psi_list[t][j%N])/gamma(i)
                a2 = (Kroncker(l-1,j) - Kroncker(l,j))*np.sin(psi_list[t][j%N])/gamma(i)
                a = a1 + a2

            if i == N:
                b1 = -(k/ds)*(psi_list[t][(i+1)%N] - psi_list[t][i%N]) + k*c0 

                a1 = Kroncker(i%N,j)*np.sin(psi_list[t][j%N])
                a2 = Kroncker(l,j)*np.cos(psi_list[t][j%N])
                a = a1 - a2

            if N < i < 2*N - 1 :
                b1 = -(k/ds)*(psi_list[t][(i+1)%N] + psi_list[t][(i-1)%N] - 2*psi_list[t][i%N]) 

                a1 = Kroncker(i%N,j)*np.sin(psi_list[t][j%N])
                a2 = Kroncker(l,j)*np.cos(psi_list[t][j%N])
                a = a1 - a2

            if i == 2*N-1:
                b1 = -(k/ds)*(psi_list[t][(i-1)%N] - 2*psi_list[t][i%N]) 

                a1 = Kroncker(i%N,j)*np.sin(psi_list[t][j%N])
                a2 = Kroncker(l,j)*np.cos(psi_list[t][j%N])
                a = a1 - a2

            A[i][j] = a
            b[i] = b1

    if print_matrix == True:
        print(f"A: {np.shape(A)[0]}x{np.shape(A)[1]}\n",A)
        print("b:",b)
       #print("x:",x)
        print("psi:",psi_list[t])

    x = np.linalg.solve(A,b)
    return x




if __name__ == "__main__":
    args = One_D_Constants(
        init_rand_psi=True
    )
    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0  =args[6:9]

    
    print(
        "\n --------------------  \n"
    )

    x = Lagran_multi(
        psi_list=psi_list
        ,t=0,k=k,c0=c0,ds=ds
        ,print_matrix=True
        ,linalg_lstsq=False
        ,num_chains=N
        )
    
 
    








    
    