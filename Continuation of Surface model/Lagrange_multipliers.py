import  numpy as np

from Two_D_functions import dSdpsi_func, gamma, Q_function

def lambda_ip1(i:int,N:int,ds:float,eta:float,r:list,z:list,psi:list,Area:list):
    if 0 <= i < N - 1:
        a1 = (z[i+1] - z[i])*np.sin(psi[i])
        a2 = 2*r[i+1]*np.cos(psi[i])
        a3 = gamma(i=i+1,ds=ds,eta=eta)*Area[i+1]*Area[i]

        return -2*np.pi**2*(a1 + a2)*r[i+1]/a3
    elif i == N - 1:
        return 0.0
    else:
        print(f"\n i not in correct range 0 <= i={i} < N -1" +f"\n lambda_ip1 func")
        exit()

def lambda_i(i:int,N:int,ds:float,eta:float,r:list,z:list,psi:list,Area:list):
    if 0 <= i < N - 1 :
        a11 = (z[i+1] - z[i])*np.sin(psi[i])
        a12 = 2*r[i+1]*np.cos(psi[i])
        a1  = (a11 + a12 )*r[i+1]/gamma(i=i+1,ds=ds,eta=eta)

        a21 = (z[i+1] - z[i])*np.sin(psi[i])
        a22 = -2*r[i]*np.cos(psi[i])
        a2  = -( a21 + a22 )*r[i]/gamma(i=i,ds=ds,eta=eta)

        return 2*np.pi**2*( a1 + a2 )/Area[i]**2
    
    elif i == N - 1:
        a21 = z[i]*np.sin(psi[i])
        a22 = 2*r[i]*np.cos(psi[i])
        a2  = ( a21 + a22 )*r[i]/gamma(i=i,ds=ds,eta=eta)

        return 2*np.pi**2*( a2 )/Area[i]**2
    else:
        print(f"\n i not in correct range 0 <= i={i} < N"+f"\n lambda_i func")
        exit()

def lambda_im1(i:int,N:int,ds:float,eta:float,r:list,z:list,psi:list,Area:list):
    if i == 0:
        return 0.0
    elif 0 < i < N -1 :
        a1 = (z[i+1] - z[i])*np.sin(psi[i])
                            
        a2 = -2*r[i]*np.cos(psi[i])

        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i-1]*Area[i]

        return 2*np.pi**2* ( a1 + a2 )*r[i]/a3
    
    elif i == N - 1:
        a1 = z[i]*np.sin(psi[i])
                            
        a2 = 2*r[i]*np.cos(psi[i])

        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i-1]*Area[i]

        return -2*np.pi**2* ( a1 + a2 )*r[i]/a3
    else:
        print(f"\n i not in correct range 1 <= i={i} < N"+f"\n lambda_im1 func")
        exit()


def nus_ip1(i:int,N:int,ds:float,eta:float,r:list,z:list,psi:list,Area:list):
    if 0 <= i < N - 2:
        a11 = (r[i+1] +r[i])*np.sin(psi[i])
        a12 = r[i+2] + r[i+1] 
        a1 = -a11*a12

        a21 = (z[i+1]-z[i])*np.sin(psi[i])
        a22 = 2*r[i+1]*np.cos(psi[i])
        a23 = z[i+2] - z[i+1]
        a2 = (a21 + a22)*a23

        a3 = gamma(i=i+1,ds=ds,eta=eta)*Area[i+1]*Area[i]
        return np.pi**2*(a1 + a2 )/a3    
    elif i == N-2:
        a11 = (r[i+1] +r[i])*np.sin(psi[i])
        a12 = r[i+2] + r[i+1] 
        a1 = -a11*a12

        a21 = (z[i+1]-z[i])*np.sin(psi[i])
        a22 = 2*r[i+1]*np.cos(psi[i])
        a23 = - z[i+1]
        a2 = (a21 + a22)*a23

        a3 = gamma(i=i+1,ds=ds,eta=eta)*Area[i+1]*Area[i]

        return np.pi**2*(a1 + a2 )/a3
    elif i == N - 1:
        return 0.0
    else:
        print(f"\n i not in correct range 1 <= i={i} < N-1"+f"\n nus_ip1 func")
        exit()

def nus_i(i:int,N:int,ds:float,eta:float,r:list,z:list,psi:list,Area:list):
    if 0 <= i < N - 1:
        a11 = (r[i+1] +r[i])**2*np.sin(psi[i])
        a12 =  1/gamma(i=i+1,ds=ds,eta=eta) + 1/gamma(i=i,ds=ds,eta=eta) 
        a1 =  a11*a12

        a21 = (z[i+1] - z[i])*np.sin(psi[i])
        a22 = 2*r[i+1]*np.cos(psi[i])
        a23 = (z[i+1] - z[i])
        a2 = (a21 + a22)*a23/gamma(i=i+1,ds=ds,eta=eta)

        a31 = (z[i+1] - z[i])*np.sin(psi[i])
        a32 = -2*r[i]*np.cos(psi[i])
        a33 = (z[i+1] - z[i])
        a3 = (a31 + a32 )*a33/gamma(i=i,ds=ds,eta=eta)

        return np.pi**2*(a1 + a2 + a3)/Area[i]**2
    elif i == N - 1:
        a11 = (r[i+1] +r[i])**2*np.sin(psi[i])
        a12 = 1/gamma(i=i,ds=ds,eta=eta) 
        a1 =  a11*a12

        a31 = z[i]*np.sin(psi[i])
        a32 = 2*r[i]*np.cos(psi[i])
        a33 =  z[i]
        a3 = (a31 + a32 )*a33/gamma(i=i,ds=ds,eta=eta)

        return np.pi**2*(a1 + a3)/Area[i]**2 
    else:
        print(f"\n i not in correct range 1 <= i={i} < N"+f"\n nus_i func")
        exit()

def nus_im1(i:int,N:int,ds:float,eta:float,r:list,z:list,psi:list,Area:list):
    if i == 0:
        return 0.0
    if 0 < i < N - 1:
        a11 = (r[i+1] + r[i] )* np.sin(psi[i])
        a12 = r[i] + r[i-1]
        a1 = - a11*a12

        a21 = (z[i+1]-z[i])*np.sin(psi[i])
        a22 = -2*r[i]*np.cos(psi[i])
        a2 = (a21 + a22)*(z[i] - z[i-1])

        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i]*Area[i-1]

        return np.pi**2*( a1 + a2 )/a3    
    elif i == N - 1:
        a11 = (r[i+1] + r[i] )* np.sin(psi[i])
        a12 = r[i] + r[i-1]
        a1 = - a11*a12

        a21 = z[i]*np.sin(psi[i])
        a22 = 2*r[i]*np.cos(psi[i])
        a2 = -(a21 + a22)*(z[i] - z[i-1])

        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i]*Area[i-1] 

        return np.pi**2*( a1 + a2 )/a3
    
    else:
        print(f"\n i not in correct range 1 <= i={i} < N"+f"\n nus_im1 func")
        exit()



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

            #---------- For the constraint eqations matrix part -------------------------------
            if 0 <= n < N:
                if 0 <= m < N:
                    #---------- Calculate the lambda values -------------------------------
                    if j == i + 1:
                        A[n,m] = lambda_ip1(i=i,N=N,ds=ds,eta=eta,r=r,z=z,psi=psi,Area=Area)

                    if j == i:
                        A[n,m] = lambda_i(i=i,N=N,ds=ds,eta=eta,r=r,z=z,psi=psi,Area=Area)

                    if j == i - 1:
                        A[n,m] = lambda_im1(i=i,N=N,ds=ds,eta=eta,r=r,z=z,psi=psi,Area=Area)
                    
                if N <= m < 2*N:
                    #---------- Calculate the nu values -------------------------------
                    if j == i + 1:
                        A[n,m] = nus_ip1(i=i,N=N,ds=ds,eta=eta,r=r,z=z,psi=psi,Area=Area)
                        
                    if j == i:                        
                        A[n,m] = nus_i(i=i,N=N,ds=ds,eta=eta,r=r,z=z,psi=psi,Area=Area)

                    if j == i-1:
                        A[n,m] = nus_im1(i=i,N=N,ds=ds,eta=eta,r=r,z=z,psi=psi,Area=Area)
            
            #---------- For the Lagrangian derivative matrix part -------------------------------
            if N <= n < 2*N:
                if 0 <= m < N:
                    if i == j:
                        A[n,m] = np.sin(psi[i])
                if N <= m < 2*N:
                    if i == j:
                        A[n,m] = -np.cos(psi[i])

        
        
        if 0 <= n < N:
            """---------- For the constraint eq's results part -------------------------------"""
            if 0 <= i < N-1:
                b11 = (z[i+1] - z[i] )*np.sin(psi[i]) 
                b12 = 2*r[i+1]*np.cos(psi[i]) 
                b1 = (b11 + b12 )/Area[i]

                Qp1 = Q_function(i=i+1,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,Area=Area,psi=psi,radi=r)

                b21 = 2*r[i]*np.cos(psi[i]) 
                b2 = (b11 - b21 )/Area[i]
                Q = Q_function(i=i,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,Area=Area,psi=psi,radi=r)

                b[n] = -np.pi*( b1*Qp1 + b2*Q )

            if i == N-1:
                b11 = z[i]*np.sin(psi[i]) 
                b21 = 2*r[i]*np.cos(psi[i]) 
                b2 = (b11 + b21 )/Area[i]

                Q = Q_function(i=i,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,Area=Area,psi=psi,radi=r)

                b[n] =  np.pi*b2*Q

        if N <= n < 2*N:
            """---------- For the Lagrangian derivative results part -------------------------------"""
            b[n] = dSdpsi_func(i=i,N=N,c0=c0,k=k,kG=kG,r=r,psi=psi,Area=Area)


    
    if print_matrix == True:
        np.set_printoptions(precision=1)
        print(f"A: {np.shape(A)[0]}x{np.shape(A)[1]}\n ",A)
        print("b:",b)
        x = np.linalg.solve(A,b)
        lamb_return = x[0:N]
        nu_return = x[N:2*N]
        print("len(lamb_return):",np.shape(lamb_return))
        print("len(nus_return):",np.shape(nu_return))
        print("lamb_return:",lamb_return)
        print("nus_return:",nu_return)
    else:
        x = np.linalg.solve(A,b)

        lamb_return = x[0:N]
        nu_return = x[N:2*N]

    return lamb_return,nu_return


def Lagrange_multipliers_V1(
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
                        a1 = (z[i+1] - z[i])*np.sin(psi[i])
                        a2 = 2*r[i+1]*np.cos(psi[i])
                        a3 = gamma(i=i+1,ds=ds,eta=eta)*Area[i+1]*Area[i]

                        A[n,m] = -2*np.pi**2 * ( a1 + a2) * r[i+1] /a3

                    if j == i:
                        """a1 = (z[i+1] - z[i])*np.sin(psi[i])
                        a2 =  r[i+1]/gamma(i=i+1,ds=ds,eta=eta) - r[i]/gamma(i=i,ds=ds,eta=eta)
                        a3 = 2*np.cos(psi[i])
                        a4 = r[i+1]**2/gamma(i=i+1,ds=ds,eta=eta) + r[i]**2/gamma(i=i,ds=ds,eta=eta)

                        A[n,m] = 2*np.pi**2 * ( a1*a2 + a3*a4 )/Area[i]**2"""

                        a11 = (z[i+1] -z[i])*np.sin(psi[i])
                        a12 = 2*r[i+1]*np.cos(psi[i])
                        a1 = (a11 + a12 )*r[i]/gamma(i=i+1,ds=ds,eta=eta)

                        a21 = (z[i+1] -z[i])*np.sin(psi[i])
                        a22 = -2*r[i]*np.cos(psi[i])
                        a2 = -(a21 +a22 )*r[i]/gamma(i=i,ds=ds,eta=eta)

                        A[n,m] = 2*np.pi**2*(a1 + a2)/Area[i]**2
                
                
                if N <= m < 2*N:
                    """---------- Calculate the nu values -------------------------------"""
                    if j == i + 1:
                        a11 = (z[i+1] - z[i])*np.sin(psi[i])
                        a12 = 2*r[i+1]*np.cos(psi[i])
                        a1 = (a11 + a12 )*(z[i+2] - z[i+1])

                        a21 = r[i+1] + r[i]
                        a22 = r[i+2] + r[i+1]
                        a2 = a21*a22*np.sin(psi[i])

                        a3 = gamma(i=i+1,ds=ds,eta=eta)*Area[i+1]*Area[i]

                        A[n,m] = np.pi**2* (a1 - a2)/a3
                        
                    if j == i:
                        a11 = ( r[i+1] + r[i] )**2
                        a12 = (z[i+1] - z[i])**2
                        a13 = 1/ gamma(i=i+1,ds=ds,eta=eta) + 1/ gamma(i=i,ds=ds,eta=eta)
                        a1 = ( a11 + a12 )*a13*np.sin(psi[i])

                        a21 = 2* (z[i+1] - z[i]) * np.cos(psi[i])
                        a22 = r[i+1]/ gamma(i=i+1,ds=ds,eta=eta) - r[i]/ gamma(i=i,ds=ds,eta=eta)
                        a2 = a21*a22

                        A[n,m] = np.pi**2*(a1 + a2)/Area[i]**2


            if 0 < n < N-2:
                if 0 <= m < N:
                    """---------- Calculate the lambda values -------------------------------"""
                    if j == i + 1:
                        """a1 = (z[i+1] - z[i])*np.sin(psi[i])
                        a2 = 2*r[i+1]*np.cos(psi[i])
                        a3 = gamma(i=i+1,ds=ds,eta=eta)*Area[i+1]*Area[i]

                        A[n,m] = -2*np.pi**2 * ( a1 + a2) * r[i+1] /a3"""

                        a1 = (z[i+1] - z[i])*np.sin(psi[i])
                        a2 = 2*r[i+1]*np.cos(psi[i])
                        a3 = gamma(i=i+1,ds=ds,eta=eta)*Area[i+1]*Area[i]
                        A[n,m] = -2*np.pi**2 * ( a1 + a2) * r[i+1] /a3


                    if j == i:
                        """a11 = (z[i+1] - z[i])*np.sin(psi[i])
                        a12 = r[i+1]/gamma(i=i+1,ds=ds,eta=eta) - r[i]/gamma(i=i,ds=ds,eta=eta) 
                        a1 =  a11*a12

                        a21 = 2*np.cos(psi[i])
                        a22 = r[i+1]**2/gamma(i=i+1,ds=ds,eta=eta) + r[i]**2/gamma(i=i,ds=ds,eta=eta) 
                        a2 = a21*a22

                        A[n,m] = 2*np.pi**2*( a1 + a2 )/Area[i]**2"""

                        a11 = (z[i+1] - z[i])*np.sin(psi[i])
                        a12 = 2*r[i+1]*np.cos(psi[i])
                        a1  = (a11 + a12 )*r[i+1]/gamma(i=i+1,ds=ds,eta=eta)

                        a21 = (z[i+1] - z[i])*np.sin(psi[i])
                        a22 = -2*r[i]*np.cos(psi[i])
                        a2  = -(a21 + a22 )*r[i]/gamma(i=i,ds=ds,eta=eta)

                        A[n,m] = 2*np.pi**2*( a1 + a2 )/Area[i]**2


                    if j == i - 1:
                        a1 = (z[i+1] - z[i])*np.sin(psi[i])
                        
                        a2 = -2*r[i]*np.cos(psi[i])

                        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i-1]*Area[i]

                        A[n,m] = 2*np.pi**2* ( a1 + a2 )*r[i]/a3

                    
                if N <= m < 2*N:
                    """---------- Calculate the nu values -------------------------------"""
                    if j == i + 1:
                        """a11 = 2*r[i+1]*np.cos(psi[i])*( z[i+1] - z[i])
                        a12 = ( z[i+1] - z[i])*np.sin(psi[i])                         
                        a1 = (a11 + a12 )*(z[i+2] - z[i+1])

                        a22 = -( r[i+1] + r[i] )*( r[i+2] + r[i+1] )
                        a2 = a22 *np.sin(psi[i]) 

                        a3 = gamma(i=i+1,ds=ds,eta=eta)*Area[i+1]*Area[i]

                        A[n,m] = np.pi**2*(a1 + a2 )/a3"""

                        a11 = (r[i+1] +r[i])*np.sin(psi[i])
                        a12 = r[i+2] + r[i+1] 
                        a1 = -a11*a12

                        a21 = (z[i+1]-z[i])*np.sin(psi[i])
                        a22 = 2*r[i+1]*np.cos(psi[i])
                        a23 = z[i+2] - z[i+1]
                        a2 = (a21 + a22)*a23

                        a3 = gamma(i=i+1,ds=ds,eta=eta)*Area[i+1]*Area[i]

                        A[n,m] = np.pi**2*(a1 + a2 )/a3

                        
                    if j == i:                        
                        """a11 = 1/gamma(i=i+1,ds=ds,eta=eta) + 1/gamma(i=i,ds=ds,eta=eta) 
                        a12 =  (r[i+1] + r[i])**2 + (z[i+1] - z[i])**2
                        a1 =  a11*a12*np.sin(psi[i])

                        a21 = 2*(z[i+1] - z[i])*np.cos(psi[i])
                        a22 = r[i+1]/gamma(i=i+1,ds=ds,eta=eta) - r[i]/gamma(i=i,ds=ds,eta=eta) 
                        a2 = a21*a22

                        A[n,m] = np.pi**2*(a1 + a2)/Area[i]**2"""

                        a11 = (r[i+1] +r[i])**2*np.sin(psi[i])
                        a12 =  1/gamma(i=i+1,ds=ds,eta=eta) + 1/gamma(i=i,ds=ds,eta=eta) 
                        a1 =  a11*a12

                        a21 = (z[i+1] - z[i])*np.sin(psi[i])
                        a22 = 2*r[i+1]*np.cos(psi[i])
                        a23 = (z[i+1] - z[i])
                        a2 = (a21 + a22)*a23/gamma(i=i+1,ds=ds,eta=eta)

                        a31 = (z[i+1] - z[i])*np.sin(psi[i])
                        a32 = -2*r[i]*np.cos(psi[i])
                        a33 = (z[i+1] - z[i])
                        a3 = (a31 + a32 )*a33/gamma(i=i,ds=ds,eta=eta)

                        A[n,m] = np.pi**2*(a1 + a2 + a3)/Area[i]**2

                    if j == i-1:
                        """a11 = (z[i+1] - z[i])*np.sin(psi[i])
                        a12 = -2*r[i]*np.cos(psi[i])
                        a13 = z[i] - z[i-1]
                        a1 = (a11 + a22)*a13

                        a2 = -(r[i+1] + r[i])*( r[i] + r[i-1] )*np.sin(psi[i])

                        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i]*Area[i-1]

                        A[n,m] = np.pi**2*( a1 + a2 )/a3"""

                        a11 = (r[i+1] + r[i] )* np.sin(psi[i])
                        a12 = r[i] + r[i-1]
                        a1 = -a11*a12

                        a21 = (z[i+1]-z[i])*np.sin(psi[i])
                        a22 = -2*r[i]*np.cos(psi[i])
                        a2 = (a21 + a22)*(z[i] - z[i-1])

                        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i]*Area[i-1]

                        A[n,m] = np.pi**2*( a1 + a2 )/a3


            if  n == N-2:
                if 0 <= m < N:
                    """---------- Calculate the lambda values -------------------------------"""
                    if j == i + 1:
                        a1 = (z[i+1] - z[i])*np.sin(psi[i])
                        a2 = 2*r[i+1]*np.cos(psi[i])
                        a3 = gamma(i=i+1,ds=ds,eta=eta)*Area[i+1]*Area[i]
                        A[n,m] = -2*np.pi**2 * ( a1 + a2) * r[i+1] /a3

                    if j == i:
                        a11 = (z[i+1] - z[i])*np.sin(psi[i])
                        a12 = 2*r[i+1]*np.cos(psi[i])
                        a1 =   (a11 + a12 )*r[i+1]/gamma(i=i+1,ds=ds,eta=eta)

                        a21 = (z[i+1] - z[i])*np.sin(psi[i])
                        a22 = -2*r[i]*np.cos(psi[i])
                        a2 =   -(a11 + a12 )*r[i]/gamma(i=i,ds=ds,eta=eta)

                        A[n,m] = 2*np.pi**2*( a1 + a2 )/Area[i]**2

                    if j == i - 1:
                        a1 = (z[i+1] - z[i])*np.sin(psi[i])
                        
                        a2 = -2*r[i]*np.cos(psi[i])

                        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i-1]*Area[i]

                        A[n,m] = 2*np.pi**2* ( a1 + a2 )*r[i]/a3


                if N <= m < 2*N:
                    """---------- Calculate the nu values -------------------------------"""
                    if j == i + 1:
                        a11 = (r[i+1] +r[i])*np.sin(psi[i])
                        a12 = r[i+2] + r[i+1] 
                        a1 = -a11*a12

                        a21 = (z[i+1]-z[i])*np.sin(psi[i])
                        a22 = 2*r[i+1]*np.cos(psi[i])
                        a23 = z[i+2] - z[i+1]
                        a2 = (a21 + a22)*a23

                        a3 = gamma(i=i+1,ds=ds,eta=eta)*Area[i+1]*Area[i]

                        A[n,m] = np.pi**2*(a1 + a2 )/a3

                    if j == i:
                        a11 = (r[i+1] +r[i])**2*np.sin(psi[i])
                        a12 =  1/gamma(i=i+1,ds=ds,eta=eta) + 1/gamma(i=i,ds=ds,eta=eta) 
                        a1 =  a11*a12

                        a21 = (z[i+1] - z[i])*np.sin(psi[i])
                        a22 = 2*r[i+1]*np.cos(psi[i])
                        a23 = (z[i+1] - z[i])
                        a2 = (a21 + a22)*a23/gamma(i=i+1,ds=ds,eta=eta)

                        a31 = (z[i+1] - z[i])*np.sin(psi[i])
                        a32 = -2*r[i]*np.cos(psi[i])
                        a33 = (z[i+1] - z[i])
                        a3 = (a31 + a32 )*a33/gamma(i=i,ds=ds,eta=eta)

                        A[n,m] = np.pi**2*(a1 + a2 + a3)/Area[i]**2

                    if j == i - 1:
                        a11 = (r[i+1] + r[i] )* np.sin(psi[i])
                        a12 = r[i] + r[i-1]
                        a1 = -a11*a12

                        a21 = (z[i+1]-z[i])*np.sin(psi[i])
                        a22 = -2*r[i]*np.cos(psi[i])
                        a2 = (a21 + a22)*(z[i] - z[i-1])

                        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i]*Area[i-1]

                        A[n,m] = np.pi**2*( a1 + a2 )/a3
            

            if  n == N-1:
                if 0 <= m < N:
                    """---------- Calculate the lambda values -------------------------------"""
                    if j == i:
                        a1 = z[i]*np.sin(psi[i])
                        
                        a2 = 2*r[i]*np.cos(psi[i])

                        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i]**2
                        
                        A[n,m] = 2*np.pi**2*( a1 + a2 )*r[i]/a3

                    if j == i - 1:
                        a1 = z[i]*np.sin(psi[i])
                        
                        a2 = 2*r[i]*np.cos(psi[i])

                        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i]*Area[i-1]
                        
                        A[n,m] = -2*np.pi**2*( a1 + a2 )*r[i]/a3


                if N <= m < 2*N:
                    """---------- Calculate the nu values -------------------------------"""
                    if j == i:
                        a11 = (r[i+1] + r[i])**2 
                        a1 = a11*np.sin(psi[i])

                        a21 = z[i]*np.sin(psi[i])
                        a22 = 2*r[i]*np.cos(psi[i])
                        a2 = (a21 + a22)*z[i]

                        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i]**2
                        
                        A[n,m] = np.pi**2* (a1 + a2)/a3 

                    if j == i - 1:
                        a11 = (r[i+1] + r[i])*np.sin(psi[i])
                        a12 = r[i] + r[i-1]
                        a1 = a11*a12

                        a21 = z[i]*np.sin(psi[i]) 
                        a22 = 2*r[i]*np.cos(psi[i]) 
                        a23 = z[i] - z[i-1]
                        a2 = (a21 + a22)*a23

                        a3 = gamma(i=i,ds=ds,eta=eta)*Area[i]*Area[i-1]

                        A[n,m] = -np.pi**2*( a1 + a2 )/a3
        
            
            
            if N <= n < 2*N:
                """---------- For the Lagrangian derivative matrix part -------------------------------"""
                if 0 <= m < N:
                    if i == j:
                        A[n,m] = np.sin(psi[i])
                if N <= m < 2*N:
                    if i == j:
                        A[n,m] = -np.cos(psi[i])

        
        
        if 0 <= n < N:
            """---------- For the constraint eq's results part -------------------------------"""
            if 0 <= i < N-1:
                b11 = (z[i+1] - z[i] )*np.sin(psi[i]) 
                b12 = 2*r[i+1]*np.cos(psi[i]) 
                b1 = (b11 + b12 )/Area[i]

                Qp1 = Q_function(i=i+1,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,Area=Area,psi=psi,radi=r)

                b21 = 2*r[i]*np.cos(psi[i]) 
                b2 = (b11 - b21 )/Area[i]
                Q = Q_function(i=i,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,Area=Area,psi=psi,radi=r)

                b[n] = -np.pi*( b1*Qp1 + b2*Q )

            if i == N-1:
                b11 = z[i]*np.sin(psi[i]) 
                b21 = 2*r[i]*np.cos(psi[i]) 
                b2 = (b11 + b21 )/Area[i]

                Q = Q_function(i=i,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,Area=Area,psi=psi,radi=r)

                b[n] =  np.pi*b2*Q

        if N <= n < 2*N:
            """---------- For the Lagrangian derivative results part -------------------------------"""
            b[n] = dSdpsi_func(i=i,N=N,c0=c0,k=k,kG=kG,r=r,psi=psi,Area=Area)


    
    if print_matrix == True:
        np.set_printoptions(precision=1)
        print(f"A: {np.shape(A)[0]}x{np.shape(A)[1]}\n ",A)
        print("b:",b)
        x = np.linalg.solve(A,b)
        lamb_return = x[0:N]
        nu_return = x[N:2*N]
        print("len(lamb_return):",np.shape(lamb_return))
        print("len(nus_return):",np.shape(nu_return))
        print("lamb_return:",lamb_return)
        print("nus_return:",nu_return)
    else:
        x = np.linalg.solve(A,b)

        lamb_return = x[0:N]
        nu_return = x[N:2*N]

    return lamb_return,nu_return



if __name__ == "__main__":
    Ns = 3
    rs = [5*i for i in range(Ns+1)]
    zs = [i for i in range(Ns+1)]
    psis = [np.pi/3.5 for i in range(Ns)]
    As = [0.1 for i in range(Ns)]

    Lagrange_multipliers(N=Ns,k=1,c0=1,sigma=1,kG=1,tau=1,ds=1,eta=1,Area=As,r=rs,z=zs,psi=psis,print_matrix=True)