import numpy as np
from One_D_Constants import One_D_Constants,gamma


    
def Kroncker(i,j):
    if i == j :
        return 1
    else:
        return 0

def Lagran_multi(
        psi_list,t,k,c0,ds
        ,linalg_lstsq =True
        ,print_matrix = False
                 ):
    
    N = int(len(psi_list[0]))
    b = np.full(shape=(2*N),fill_value=10,dtype=float)
    A = np.full(shape=(2*N,2*N),fill_value=10,dtype=float)

    for i in range(2*N):
        for j in range(2*N):
            if i < N:
                b1 = 0
                if i == 0:
                    a1 =  (Kroncker(i,j) /gamma(i%N)- Kroncker(i+1,j) )/gamma((i+1)%N)*np.cos(psi_list[t][j%N])
                    a2 =(Kroncker(i+N,j) /gamma(i%N)- Kroncker(i+N+1,j) /gamma((i+1)))*np.sin(psi_list[t][j%N])
                
                else:
                    a11 = ( Kroncker(i,j) - Kroncker(i+1,j) )/gamma((i+1)%N)
                    a12 =( Kroncker(i-1,j) - Kroncker(i,j) )/gamma(i%N)
                    a1 = ( a11 - a12 )*np.cos(psi_list[t][j%N])
                    
                    a21 = ( Kroncker((i+N),j) - Kroncker(i+N+1,j) )/gamma((i+1)%N)            
                    a22 = ( Kroncker(i+N-1,j) - Kroncker(i+N,j) )/gamma(i%N)
                    a2 = (a21 - a22 )*np.sin(psi_list[t][j%N])                
            
            if i>=N:
                if i==N:
                    b1 = k*(psi_list[t][(i+1)%N]-psi_list[t][i%N]) -k*ds*c0
                    
                    a1 = np.sin(psi_list[t][j%N])*Kroncker(i,j)
                    a2 = - np.cos(psi_list[t][j%N])*Kroncker(i+N,j)

                if i > N:
                    b1 = k*(psi_list[t][(i+1)%N] + psi_list[t][(i-1)%N] - psi_list[t][i%N])

                    a1 = np.sin(psi_list[t][j%N])*Kroncker(i%N,j)
                    a2 = - np.cos(psi_list[t][j%N])*Kroncker(i%N+N,j)

            A[i][j] = a1+a2 #round(a1+a2,2)
            b[i] = b1 #round(b1,2)


    if print_matrix == True:
        print("A:",A)
        print("b:",b)
    if linalg_lstsq == True:
        x = np.linalg.lstsq(A,b,rcond=None)[0]
    else:
        x = np.linalg.solve(A,b)

    return x


def dPsidt(i,t,dt,multi,psi,deltaS):
    
    N = int(len(multi)/2)
    # from 0 -> N-1 is the lambda Lagrangian multipliers
    # from N -> 2N is the nu Lagrangian multipliers
    k = i + N
    if i == 0:
        a1 = (multi[i+1]/gamma(i+1) - multi[i]/gamma(i))*(np.sin(psi[t][i])/deltaS)
        a2 =(multi[k+1]/gamma(i+1) - multi[k]/gamma(i))*(np.cos(psi[t][i])/deltaS)

    else:
        a11 = (multi[i] - multi[i+1])/gamma(i+1) 
        a12 =- (multi[i-1] - multi[i])/gamma(i)
        a1 = (a11+a12)*(np.sin(psi[t][i])/deltaS)

        a21 =(multi[k] - multi[k+1])/gamma(i+1) 
        a22 = - (multi[k-1] - multi[k])/gamma(i)
        a2 = (a21 + a22)*(np.sin(psi[t][i])/deltaS)

    # updating the new psi.
    psi[t+1][i] = psi[t][i] + dt*(a1 +a2)
    



if __name__ == "__main__":
    args = One_D_Constants()
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
        )

    #print("x:",x)






    
    