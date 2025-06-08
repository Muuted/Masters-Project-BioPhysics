import numpy as np
from One_D_Constants import One_D_Constants,gamma


    
def Kroncker(i,j):
    if i == j :
        return 1
    else:
        return 0
"""
def Lagran_multi(
        psi_list,t,k,c0,ds
        ,num_chains
        ,linalg_lstsq =True
        ,print_matrix = False
                 ):
    
    N = num_chains #+ 1
    NN = N +1
    b = np.full(shape=(2*N),fill_value=10,dtype=float)
    A = np.full(shape=(2*N,2*N),fill_value=10,dtype=float)

    for i in range(2*N):
        l = i%N + N
        b1 = 0
        for j in range(2*N):
            a1,a11,a12 = 0,0,0
            a2,a21,a22 = 0,0,0
            #b1 = 0
            if i < N:
                if i == 0:
                    if j < N:
                        a1 = (  Kroncker(i,j)/gamma(i%N)  -  Kroncker(i+1,j)/gamma((i+1)%N)  )*np.cos(psi_list[t][j%NN])
                    if j >= N:
                        a2 = (  Kroncker(l,j)/gamma(l%N) - Kroncker(l+1,j)/gamma((l+1)%N) )*np.sin(psi_list[t][j%NN])
                
                if 0 < i < N:
                    if j < N:
                        a11 = ( Kroncker(i,j) - Kroncker(i+1,j) )/gamma((i+1)%N)
                        a12 = -( Kroncker(i-1,j) - Kroncker(i,j) )/gamma(i%N)
                        a1 = ( a11 + a12 )*np.cos(psi_list[t][j%NN]) 
                    if j >= N:
                        a21 = ( Kroncker(l,j) - Kroncker(l+1,j) )/gamma((l+1)%N)            
                        a22 = -( Kroncker(l-1,j) - Kroncker(l,j) )/gamma(l%N)
                        a2 = (a21 + a22 )*np.sin(psi_list[t][j%NN])                           
            
            if i>=N:
                if i==N:
                    b1 = 1*(-(k/ds)*(psi_list[t][(i+1)%NN] - psi_list[t][i%NN]) + k*c0 )
                    
                    a1 = np.sin(psi_list[t][j%NN])*Kroncker(i%N,j)
                    a2 = -np.cos(psi_list[t][j%NN])*Kroncker(l,j)

                if i > N:
                    b1 = 1*(-(k/ds)*(psi_list[t][(i+1)%NN] + psi_list[t][(i-1)%NN] - 2*psi_list[t][i%NN]) )

                    a1 = np.sin(psi_list[t][j%NN])*Kroncker(i%N,j)
                    a2 = -np.cos(psi_list[t][j%NN])*Kroncker(l,j)
            
                
            if round(a1+a2,1) == 2.9:
                print("(i,j)=",i,j)
                print(f"a11={a11}")
                print(f"a12={a12}")
                print(f"a1={a1}")
                print(f"a2={a2}")
            A[i][j] = a1+a2 #round(a1+a2,1)
        b[i] = b1 #round(b1,1)

    if linalg_lstsq == True:
        x = np.linalg.lstsq(A,b,rcond=None)[0]
    else:
        x = np.linalg.solve(A,b)
    
    if print_matrix == True:
        print(f"A: {np.shape(A)[0]}x{np.shape(A)[1]}\n",A)
        print("b:",b)
        
    return x

"""
def dPsidt(i,N,multi,psi,deltaS):
    l = i + N #+ 1
    if i == 0:
        a1 = (multi[i+1]/gamma(i+1) - multi[i]/gamma(i))*np.sin(psi)/deltaS**2
        a2 = (multi[l]/gamma(i) - multi[l+1]/gamma(i+1) )*np.cos(psi)/deltaS**2

    if i > 0:
        a11 = (multi[i-1] - multi[i])/gamma(i)
        a12 = (multi[i] - multi[i+1])/gamma(i+1) 
        a1 = (a11 - a12)*np.sin(psi)/deltaS**2

        a21 = (multi[l] - multi[l+1])/gamma(i+1) 
        a22 = (multi[l-1] - multi[l])/gamma(i)
        a2 = (a21 - a22)*np.cos(psi)/deltaS**2

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
    NN = N + 1
    b = np.full(shape=(2*NN),fill_value=0,dtype=float)
    A = np.full(shape=(2*NN,2*NN),fill_value=10,dtype=float)

    for i in range(2*(N+1)):
        for j in range(2*(N+1)):
            l = i%NN + NN
            if i == 0:
                b[i] = 0

                if j < N + 1:
                    a1 = Kroncker(i,j)/gamma(i%NN) - Kroncker(i+1,j)/gamma((i+1)%NN)
                    a2 = 0
                if j >= N +1 :
                    a2 = Kroncker(l,j)/gamma(i%NN) - Kroncker(l+1,j)/gamma((i+1)%NN)
                    a1 = 0
                A[i][j] = a1*np.cos(psi_list[t][j%NN]) + a2*np.sin(psi_list[t][j%NN])

            if 0 < i < N+1 :
                b[i] = 0
                
                if j < N + 1:
                    a11 = (Kroncker(i,j) - Kroncker(i+1,j))/gamma((i+1)%NN)
                    a12 = (Kroncker(i-1,j) - Kroncker(i,j))/gamma(i%NN)
                    a1 = (a11 - a12)*np.cos(psi_list[t][j%NN])
                    a2 = 0
                if j >= N + 1:
                    a21 = (Kroncker(l,j) - Kroncker(l+1,j))/gamma((i+1)%NN)
                    a22 = (Kroncker(l-1,j) - Kroncker(l,j))/gamma(i%NN)
                    a2 = (a21 - a22)*np.sin(psi_list[t][j%NN])
                    a1 = 0

                A[i][j] = a1 + a2
            
            if i == N+1:
                b[i] = 1*( -(k/ds)*(psi_list[t][(i+1)%NN] - psi_list[t][i%NN]) + k*c0 )

                a1 = Kroncker(i%NN,j)*np.sin(psi_list[t][j%NN])
                a2 = Kroncker(l,j)*np.cos(psi_list[t][j%NN])
                A[i][j] = a1 - a2

            if i > N+1 :
                b[i] = - (k/ds)*(psi_list[t][(i+1)%NN] + psi_list[t][(i-1)%NN] - 2*psi_list[t][i%NN]) 

                a1 = Kroncker(i%NN,j)*np.sin(psi_list[t][j%NN])
                a2 = Kroncker(l,j)*np.cos(psi_list[t][j%NN])
                A[i][j] = a1 - a2


    x = np.linalg.solve(A,b)

    if print_matrix == True:
        print(f"A: {np.shape(A)[0]}x{np.shape(A)[1]}\n",A)
        print("b:",b)
        print("x:",x)
        print("psi:",psi_list[t])

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

    x = Lagran_multi_V2(
        psi_list=psi_list
        ,t=0,k=k,c0=c0,ds=ds
        ,print_matrix=True
        ,linalg_lstsq=False
        ,num_chains=2
        )
    exit()
    x1 = Lagran_multi(
        psi_list=psi_list
        ,t=0,k=k,c0=c0,ds=ds
        ,print_matrix=True
        ,linalg_lstsq=False
        ,num_chains=3
        )
 
    








    
    