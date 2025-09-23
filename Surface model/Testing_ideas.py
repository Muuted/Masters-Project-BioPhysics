import numpy as np
from Two_D_functions import Kronecker, c_diff_f,c_diff_g,constraint_f,constraint_g

def c_diff(
        i:int,j:int,N:int
        ,r:list,z:list,psi:list,Area:list
        ,diff_var =""
        ):
    
    diff_var_list = ["r","z","psi",0,1,2]
    if diff_var not in diff_var_list:
        print(f"c diff error in diff_var")
        exit()
    c_diff_val = ""
    if i < N :
        c_diff_val =c_diff_f(
            i=i%N,j=j%N,N=N
            ,r=r,psi=psi,Area=Area
            ,diff_var=diff_var
            )
    if i >= N :
        c_diff_val =c_diff_g(
            i=i%N,j=j%N,N=N
            ,r=r,z=z,psi=psi,Area=Area
            ,diff_var=diff_var
            )

    if c_diff_val == "":
        print(f"Error c_diff_val didnt take a value because i={i}")
        exit()
    
    return c_diff_val




def Epsilon_v2(
        N:int,r:list,z:list,psi:list,Area:list
        ,print_matrix = False
        ):
    A = np.zeros(shape=(2*N,2*N))
    b = np.zeros(2*N)
    for alpha in range(2*N):
        for beta in range(2*N):
            a = 0
            for n in range(N):
                for variables in range(3):
                    A[alpha][beta] += (
                        c_diff(
                            i=alpha,j=n,N=N,r=r,z=z,psi=psi,Area=Area,diff_var=variables
                        )*c_diff(
                            i=beta,j=n,N=N,r=r,z=z,psi=psi,Area=Area,diff_var=variables
                            )
                        )
            

        if alpha < N :
            b[alpha] = -constraint_f(i=alpha%N,N=N,r=r,psi=psi,Area=Area)
        if alpha >= N :
            b[alpha] = - constraint_g(i=alpha%N,N=N,r=r,z=z,psi=psi,Area=Area)

    
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

    return epsilon_f,epsilon_g
