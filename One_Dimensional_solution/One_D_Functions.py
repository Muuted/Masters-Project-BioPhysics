import numpy as np
from One_D_Constants import One_D_Constants,gamma


    
def Kroncker(i,j):
    if i == j :
        return 1
    else:
        return 0

def Lagran_multi(psi_list,t):
    
    N = int(len(psi_list[0]))

    b = np.full(shape=(2*N),fill_value=0)
    A = np.full(shape=(2*N,2*N),fill_value=10,dtype=float)

    row1,row2 = 0,N
    for i in range(2*N):
        for j in range(2*N):
            if i%N == 0:
                a1 =  (Kroncker(row1,j) /gamma(row1)- Kroncker(row1+1,j) /gamma(row1+1))*np.cos(psi_list[t][j%N])
                a2 =(Kroncker(row2,j) /gamma(row2)- Kroncker(row2+1,j) /gamma(row2+1))*np.sin(psi_list[t][j%N])

            else:
                a11 = ( Kroncker(row1,j) - Kroncker(row1+1,j) )/gamma(row1+1)
                a12 =( Kroncker(row1-1,j) - Kroncker(row1,j) )/gamma(row1)
                a1 = ( a11 - a12 )*np.cos(psi_list[t][j%N])
                
                a21 = ( Kroncker(row2,j) - Kroncker(row2+1,j) )/gamma(row2+1)            
                a22 = ( Kroncker(row2-1,j) - Kroncker(row2,j) )/gamma(row2)
                a2 = (a21 - a22 )*np.sin(psi_list[t][j%N])      


            A[i][j] =  round(a1 + a2,2) # a1 + a2
        
        row1 += 1
        row2 += 1
        if row1%N == 0 :
            row1 = 0
        if row2%N == 0 :
            row2 = N
        
    #print(A)
    print(np.shape(A))
    x = np.linalg.lstsq(A,b,rcond=None)
    #x = np.linalg.solve(A,b)

    return x



def dPsidt(i,multi,psi,deltaS):
    N = len(multi)/2
    # from 0 -> N-1 is the lambda Lagrangian multipliers
    # from N -> 2N is the nu Lagrangian multipliers
    k = i + N
    if i == 0:
        a1 = (multi[i+1]/gamma(i+1) - multi[i]/gamma(i))*(np.sin(psi[i])/deltaS)
        a2 =(multi[k+1]/gamma(i+1) - multi[k]/gamma(i))*(np.cos(psi[i])/deltaS)

    if 0<i<=N:
        a11 = (multi[i] - multi[i+1])/gamma(i+1) 
        a12 =- (multi[i-1] - multi[i])/gamma(i)
        a1 = (a11+a12)*(np.sin(psi[i])/deltaS)

        a21 =(multi[k] - multi[k+1])/gamma(i+1) 
        a22 = - (multi[k-1] - multi[k])/gamma(i)
        a2 = (a21 + a22)*(np.sin(psi[i])/deltaS)

    else:
        print("error")

    return a1 + a2



if __name__ == "__main__":
    args_list = One_D_Constants()
    L,r0,N,ds,T,dt = args_list[0:6]
    x_list ,z_list ,psi_list = args_list[6:9]
    lambda_list ,nu_list = args_list[9:11]
    
    print(
        "\n --------------------  \n"
    )

    x = Lagran_multi(
        psi_list=psi_list
        ,t=0
        )

    print("x:",x[0])






    
    