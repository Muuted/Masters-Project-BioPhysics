import numpy as np
import matplotlib.pyplot as plt
from One_D_Constants import One_D_Constants,gamma
from One_D_Functions import Lagran_multi, dPsidt_RungeKutta_4
from plotting_functions import plot_from_psi_V2 
import os
import pandas as pd
import progressbar
from Make_movie import Make_frames,Make_video


def cirle_fit(data_path, df_name):
    df_sim = pd.read_pickle(data_path + df_name)
    steps_tot = df_sim["sim_steps"][0]
    x = df_sim['x pos'][0][steps_tot-1]
    z = df_sim['z pos'][0][steps_tot-1]
    x_max ,x_min = max(x), min(x)
    z_max ,z_min = max(z), min(z)
    # coordinates of the barycenter
    x_m = np.mean(x)
    z_m = np.mean(z)

    # calculation of the reduced coordinates
    u = x - x_m
    v = z - z_m

    # linear system defining the center in reduced coordinates (uc, vc):
    #    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
    #    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
    Suv  = sum(u*v)
    Suu  = sum(u**2)
    Svv  = sum(v**2)
    Suuv = sum(u**2 * v)
    Suvv = sum(u * v**2)
    Suuu = sum(u**3)
    Svvv = sum(v**3)

    # Solving the linear system
    A = np.array([ [ Suu, Suv ], [Suv, Svv]])
    B = np.array([ Suuu + Suvv, Svvv + Suuv ])/2.0

    uc, vc = np.linalg.solve(A, B)

    xc_1 = x_m + uc
    zc_1 = z_m + vc

    # Calculation of all distances from the center (xc_1, yc_1)
    Ri_1      = np.sqrt((x-xc_1)**2 + (z-zc_1)**2)
    R_1       = np.mean(Ri_1)
    residu_1  = np.sum((Ri_1-R_1)**2)
    residu2_1 = np.sum((Ri_1**2-R_1**2)**2)
    
    x_circle ,z_circle = make_circle(
        xc=xc_1,zc=zc_1,R=R_1,ds=ds
        ,xlim=[x_max,x_min]
        ,zlim=[z_max, z_min]
    )

    df_sim["x circle points"] = [x_circle]
    df_sim["z circle points"] = [z_circle]
    df_sim["x circle center"] = [xc_1]
    df_sim["z circle center"] = [zc_1]
    df_sim["circle radius"] = [R_1]

    df_sim.to_pickle(data_path + df_name)

    return [xc_1 , zc_1 , R_1 , x_circle, z_circle]


def make_circle(xc,zc,R,ds,xlim,zlim):
    xmax,xmin = xlim
    zmax,zmin = zlim
    x_list,z_list = [],[]
    x,z = xmin, zmax

    steps_size = ds/30
    tol = steps_size #Tolerence for divation of radius
    on_circ = False
    while x <= xmax:
        on_circ = False
        r =  (x-xc)**2 + (z-zc)**2
        if R**2*(1-tol) < r < R**2*(1+tol) :
            x_list.append(x)
            z_list.append(z)
            x += steps_size
            on_circ = True
        
        if on_circ == False:
            z -= steps_size

    return [x_list,z_list]



if __name__ == "__main__":
    args = One_D_Constants()
    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    df_name= args[15]

    x_cen,z_cen,Radius,x_circle ,z_circle  = cirle_fit(
        data_path=data_path
        ,df_name=df_name
    )

    df_sim = pd.read_pickle(data_path + df_name)
    x = df_sim['x pos'][0][sim_steps-1]
    z = df_sim['z pos'][0][sim_steps-1]
    x_max ,x_min = max(x), min(x)
    z_max ,z_min = max(z), min(z)
    
    print(df_sim.info())
    plt.figure()

    plt.plot(x,z,label="data")
    plt.plot(x_circle,z_circle,label="circle")
    plt.xlim([x_min-ds,x_max])
    plt.ylim([-ds*10,ds*10])    
    plt.legend()
    plt.show()
