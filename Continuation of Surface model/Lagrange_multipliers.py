import  numpy as np

from Two_D_functions import dSdpsi_func, gamma

def find_lambda(i:int,j:int,N:int,r:list,z:list,psi:list,A:list):
    if i >= N or j >= N:
        print(f"\n The either of i,j is larger the max index value (i,j)={i,j}, program terminates \n ")
        exit()

    if j == i + 1:
        pass
    elif j == i:
        pass
    elif j == i-1:
        pass
    else:
        return 0

    pass


def Lagrange_multipliers(
        N:int ,k:float ,c0:float ,sigma:float 
        ,kG:float ,tau:float , ds:float, eta:float
        ,Area:list,psi:list
        ,r:list
        ,z:list
        ,print_matrix = False
    ):
    A = np.zeros(shape=(2*N,2*N))
    b = np.zeros(2*N)

    for n in range(2*N):
        for m in range(2*N):
            i,j = n%N , m%N
            

            """---------- For the constraint eqations matrix part -------------------------------"""
            if n == 0:
                if 0 <= m < N:
                    """---------- Calculate the lambda values -------------------------------"""
                    if j == i + 1:
                        a1 = (z[j+1] - z[j])*np.sin(psi[j])
                        a2 = 2*r[j+1]*np.cos(psi[j])
                        a3 = gamma(i=j,ds=ds,eta=eta)*Area[j+1]*Area[j]

                        A[n,m] = 2*np.pi**2 * ( a1 + a2) * r[j+1] /a3

                    if j == i:
                        a1 = (z[j+1] - z[j])*np.sin(psi[j])/Area[j]
                        a2 =  r[j+1]/ gamma(i=j+1,ds=ds,eta=eta) - r[j]/ gamma(i=j,ds=ds,eta=eta)
                        a3 = 2*np.cos(psi[j])/Area[j]
                        a4 = r[j+1]**2/gamma(i=j+1,ds=ds,eta=eta) + r[j]**2/ gamma(i=j,ds=ds,eta=eta)

                        A[n,m] = 2*np.pi**2 * ( a1*a2 + a3*a4 )/Area[j]
                
                
                if N <= m < 2*N:
                    """---------- Calculate the nu values -------------------------------"""
                    if j == i + 1:
                        a11 = (z[j+1] - z[j])*np.sin(psi[j])
                        a12 = 2*r[j+1]*np.cos(psi[j])
                        a1 = (a11 + a12 )*(z[j+2] - z[j+1])

                        a21 = r[j+1] + r[j]
                        a22 = r[j+2] + r[j+1]
                        a2 = a21*a22*np.sin(psi[j])

                        a3 = gamma(i=j,ds=ds,eta=eta)*Area[j+1]*Area[j]

                        A[n,m] = np.pi**2* (a1 - a2)/a3
                    if j == i:
                        a11 = ( r[j+1] + r[j] )**2
                        a12 = (z[j+1] - z[j])**2
                        a13 = 1/ gamma(i=j+1,ds=ds,eta=eta) + 1/ gamma(i=j,ds=ds,eta=eta)
                        a1 = ( a11 + a12 )*a13*np.sin(psi[j])

                        a21 = 2* (z[j+1] - z[j]) * np.cos(psi[j])
                        a22 = r[j+1]/ gamma(i=j+1,ds=ds,eta=eta) - r[j]/ gamma(i=j,ds=ds,eta=eta)
                        a2 = a21*a22

                                         
            

            if 0 < n < N-2:
                if 0 <= m < N:
                    """---------- Calculate the lambda values -------------------------------"""
                    if j == i + 1:
                        pass
                    if j == i:
                        pass
                    if j == i-1:
                        pass
                    
                if N <= m < 2*N:
                    """---------- Calculate the nu values -------------------------------"""
                    if j == i + 1:
                        pass
                    if j == i:
                        pass
                    if j == i-1:
                        pass

            
            if  n == N-2:
                if 0 <= m < N:
                    """---------- Calculate the lambda values -------------------------------"""
                    if j == i + 1:
                        pass
                    if j == i:
                        pass
                    if j == i - 1:
                        pass


                if N <= m < 2*N:
                    """---------- Calculate the nu values -------------------------------"""
                    if j == i + 1:
                        pass
                    if j == i:
                        pass
                    if j == i - 1:
                        pass
            

            if  n == N-1:
                if 0 <= m < N:
                    """---------- Calculate the lambda values -------------------------------"""
                    if j == i + 1:
                        pass
                    if j == i:
                        pass
                    if j == i - 1:
                        pass


                if 0 <= m < N:
                    """---------- Calculate the nu values -------------------------------"""
                    if j == i + 1:
                        pass
                    if j == i:
                        pass
                    if j == i - 1:
                        pass



            if 0 <= n < N:
                """---------- For the constraint eq's results part -------------------------------"""
                if i==0:
                    pass
                if 0 < i < N-2:
                    pass
                if i == N-2:
                    pass
                if i == N-1:
                    pass
        
            
            if N <= n < 2*N:
                """---------- For the Lagrangian derivative matrix part -------------------------------"""
                if 0<= m < N:
                    if i == j:
                        A[n,m] = np.sin(psi[j])
                if N <= m < 2*N:
                    if i == j:
                        A[n,m] = -np.cos(psi[j])


        if N <= n < 2*N:
            """---------- For the Lagrangian derivative results part -------------------------------"""
            b[n] = dSdpsi_func(i=i,N=N,c0=c0,k=k,kG=kG,r=r,psi=psi,Area=Area)


            

if __name__ == "__main___":
    pass