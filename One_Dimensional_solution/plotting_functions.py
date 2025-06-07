import numpy as np
from One_D_Constants import One_D_Constants
import matplotlib.pyplot as plt
import progressbar
import pandas as pd



def plot_from_psi(psi:list,sim_steps: int,ds: float,r0:float, L:float):
    N = len(psi[0])
    x = np.full(shape=(sim_steps,N),dtype=float,fill_value=0)
    z = np.zeros(shape=(sim_steps,N),dtype=float) #Because we have N+1 points, but run dynamics on N
    print(np.shape(x))
    x[:,N-1] = L + r0
    print_list = []
    print_list2 = []
    print("\n Getting x,z progressbar")
    b2 = progressbar.ProgressBar(maxval=sim_steps)
    for t in range(sim_steps):
        b2.update(t)
        tolerance= 1e-10 #1e-4
        for i in range(N-2,-1,-1):
            x[t][i] = x[t][i+1] - ds*np.cos(psi[t][i])
            z[t][i] = z[t][i+1] + ds*np.sin(psi[t][i])
            
            a = np.sqrt((x[t][i+1]-x[t][i])**2 + (z[t][i+1]-z[t][i])**2)
            if ds*(1+tolerance) <= a  <= ds*(1-tolerance) :
                print(f"error on constant length , x[t][i]={x[t][i]} and z[t][i]={x[t][i]}")
                exit()

    return x,z
 


if __name__ == "__main__":
    args = One_D_Constants()
    
    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    df_name = args[15]
    df_sim = pd.read_pickle(data_path + df_name)
    psi_list = df_sim["Psi"][0]
    z = df_sim["z pos"][0]
    x = df_sim["x pos"][0]

    
    x,z = plot_from_psi(
        psi=psi_list
        ,sim_steps=sim_steps
        ,ds=ds ,r0=r0 ,L=L
        )


    xz_diff = [ (x[t][1]-x[t][0])**2 + (z[t][1]-z[t][0])**2 - ds**2 for t in range(sim_steps) ]
    t = [i for i in range(sim_steps)]
    plt.figure()
    plt.plot(t,xz_diff,'-.')

    plt.figure()
    plt.plot(z[:][0])
    plt.show()
    
    iter = [i for i in range(N-1)]
    val = [ (x[0][i+1]-x[0][i])**2 + (z[0][i+1]-z[0][i])**2  for i in range(0,N-1) ]

    plt.figure()
    plt.plot(iter,val,'-.')
    plt.title(f"ds^2 - x^2 - z^2 = 0, or it should be \n ds={ds**2}") 

    plt.figure()
    plt.plot(x[0],z[0],'-*')
    plt.title("plot of x,z positions")
    plt.show()