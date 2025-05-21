import  numpy as np
from One_D_Constants import One_D_Constants
import matplotlib.pyplot as plt

def plot_from_psi(psi,ds,r0):
    N = len(psi)
    x,z = np.zeros(N),np.zeros(N)
    L = N*ds + r0
    x[N-1] = L
    z[N-1] = 0
    for i in range(N-2,-1,-1):
        x[i] = x[i+1] - ds*np.cos(psi[i])
        z[i] = z[i+1] + ds*np.sin(psi[i])
    return x,z



if __name__ == "__main__":
    args = One_D_Constants()
    L,r0,N,ds,T,dt = args[0:6]
    x_list ,z_list ,psi_list =args[6:9]
    lambda_list ,nu_list = args[9:11]

    x,z = plot_from_psi(psi=psi_list[0],ds=ds,r0=r0)

    iter = [i for i in range(N-1)]
    val = [ (x[i+1]-x[i])**2 + (z[i+1]-z[i])**2 for i in range(0,N-1) ]

    plt.figure()
    plt.plot(iter,val,'-*')
    plt.title(f"ds^2 - x^2 - z^2 = 0, or it should be \n ds={ds**2}")

    plt.figure()
    plt.plot(x,z,'-*')
    plt.show()