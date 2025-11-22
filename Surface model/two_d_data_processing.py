import numpy as np
#from Two_D_constants import Two_D_Constants, gamma, Two_D_paths, mass
#from Two_D_functions import Delta_s
import matplotlib.pyplot as plt
import os
import pandas as pd
import progressbar

def mass(i:int ,Area:list):
    rho = 1
    m = -1
    if i == 0:
        m = rho*Area[i]/2
    if i > 0:
        m = rho*( Area[i]+ Area[i-1])/2
    if m == -1:
        print("We have got the wrong correct mass, in the mass function")
        exit()
    return m

def check_area(
        t:int,N:int,r:list,z:list,Area:list
        ,tolerence:float=1e-10
        ):
    error = False
    #print(np.shape(Area))
    for i in range(N):
        area_change = np.pi*( r[i+1]+ r[i] )*np.sqrt( (r[i+1]- r[i])**2 + (z[i+1]- z[i])**2 )

        if Area[i] != area_change :
            print(
                f"we had a change in area \n"
                +f"Area[{i}]={Area[i]}  and area_change[{i}]={area_change} \n"
                +f"at time_step={t}  and position i={i} \n"
                f" Delta A = {Area[i] - area_change}"
                )
            error = True
            break

    return error

def check_area_from_data(
        df_name:str
        ,data_path:str
        ):
    
    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    Area = df_sim['area list'][0]
    dt = df_sim['dt'][0]
    L = df_sim["L"][0]
    tau = df_sim["tau"][0]
    sigma = df_sim["sigma"][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    gam2 = df_sim["gam(i>0)"][0]
    sim_steps = df_sim["sim_steps"][0]
    r0 = df_sim["r0"][0]
    #T_tot = df_sim["Total time [sec]"][0]

    area_change = np.zeros(shape=(sim_steps,N),dtype=float)
    tolerence = 1e-10

    error = False
    for t in range(1,sim_steps):
        for i in range(N):
            area_change[t][i] =np.pi*( r[t][i+1]**2 - r[t][i]**2 )

            if Area[i] != area_change[t][i] :
                print(
                    f"we had a change in area \n"
                    +f"Area[{i}]={Area[i]}  and area_change[{i}]={area_change[t][i]} \n"
                    +f"at time_step={t}  and position i={i} \n"
                    f" Delta A = {Area[i] - area_change[t][i]}"
                    )
                error = True
                break#exit()
        if error == True:
            break

    """print(
        f"\n check all the areas and no change in area was found. with dt={dt}"
    )"""
    return error
    
def rotate_coords(
         df_name:str
        ,data_path:str
    ):
    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    N = df_sim["N"][0]
    sim_steps = df_sim["sim_steps"][0]

    r = df_sim['r'][0][sim_steps-1] # use this as reference maybe?

    x = df_sim['r'][0][sim_steps-1] # calling it x, initially
    z = df_sim['z'][0][sim_steps-1]
    y = np.zeros(N+1)

    pos_vec = [np.zeros(3) for i in range(N+1)]# 3 pos for x,y,z and then the number of time steps in the other direction

    for i in range(N+1):
        pos_vec[i][0] = x[i]
        pos_vec[i][1] = y[i]
        pos_vec[i][2] = z[i]

    num_angles = 20
    phi = np.zeros(num_angles)
    dphi = 2*np.pi/num_angles
    for i in range(num_angles):
        phi[i] = i*dphi
    

def tot_area(
        N:int,r:list,z:list
        ):
    Area = 0
    for i in range(N):
        Area += np.pi*( r[i+1]+ r[i] )*np.sqrt( (r[i+1]- r[i])**2 + (z[i+1]- z[i])**2 )
    return Area


def E_pot_original(
        N:int, k:float, kG:float ,sigma:float 
        ,tau:float ,c0:float
        ,r:list,z:list,psi:list
        ,Area:list
        ):
    Atot = np.sum(Area)
    Epot = tau*r[0] - sigma*Atot/(2*np.pi)
    for i in range(N):
        if i < N -1:
            Epot += (
                (k*Area[i]/(2*np.pi*(r[i+1]+r[i])))*(
                    np.pi*(psi[i+1]-psi[i])*(r[i+1]+r[i])/Area[i] + np.sin(psi[i])/r[i] - c0
                    )**2*r[i]
                + sigma*Area[i]*r[i]/(np.pi*(r[i+1]+r[i]))
                + kG*(psi[i+1]-psi[i])*np.sin(psi[i])
            )
        if i == N -1 :
            Epot += (
                (k*Area[i]/(2*np.pi*(r[i+1]+r[i])))*(
                    np.pi*(-psi[i])*(r[i+1]+r[i])/Area[i] + np.sin(psi[i])/r[i] - c0
                    )**2*r[i]
                + sigma*Area[i]*r[i]/(np.pi*(r[i+1]+r[i]))
                + kG*(-psi[i])*np.sin(psi[i])
            )

    return Epot

def E_pot(
        N:int, k:float, kG:float
        ,tau:float ,c0:float
        ,r:list,psi:list
        ,Area:list
        )->float:
    Epot = tau*r[0]
    for i in range(N):
        if i < N - 1:
            Epot += (
                (k*Area[i]/(4*np.pi))*(
                    np.pi*(psi[i+1]-psi[i])*(r[i+1]+r[i])/Area[i] + np.sin(psi[i])/r[i] - c0
                    )**2
                + kG*(psi[i+1]-psi[i])*np.sin(psi[i])
            )
        if i == N - 1 :
            Epot += (
                (k*Area[i]/(4*np.pi))*(
                    -np.pi*psi[i]*(r[i+1]+r[i])/Area[i] + np.sin(psi[i])/r[i] - c0
                    )**2
                - kG*psi[i]*np.sin(psi[i])
            )

    return Epot

def E_kin(
        N:int, t:int, dt:float
        ,r:list ,z:list ,Area:list
        )->float:
    Ekin = 0
    for i in range(N):
        m = mass(i=i,Area=Area)

        dot_r = (r[t+1][i] - r[t][i])/dt
        dot_z = (z[t+1][i] - z[t][i])/dt

        Ekin += m*( dot_r**2 + dot_z**2 )

    return Ekin



def Xsqaured_test(N:int
                  ,r_init:list,z_init:list,psi_init:list
                  ,r:list,z:list,psi:list)->float:
    X = 0
    for i in range(N):
        X += (r_init[i] - r[i])**2 + (z_init[i] - z[i])**2 #+ (psi_init[i] - psi[i])**2
        #The angles are removed as they are between -pi and pi, which could completely dominate
        #the values for X given that r and z could potentially be much smaller.
    return X
    
def Excess_Area(rmax:float,Area_tot:float):
    return Area_tot - np.pi*rmax**2





def cirle_fit(
        data_path ,df_name 
        ,edge_point ,end_point
        ,save_data = False
        ):

    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    steps_tot = df_sim["sim_steps"][0]
    ds = df_sim["ds"][0]

    Radius_list = []
    x_center_list = []
    z_center_list = []
    for sim_step in range(1,steps_tot):
        #pass
        #x = df_sim['x pos'][0][sim_step][edge_point:end_point]
        #z = df_sim['z pos'][0][sim_step][edge_point:end_point]
        x = df_sim['r'][0][sim_step][edge_point:end_point]
        z = df_sim['z'][0][sim_step][edge_point:end_point]

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
    
        Radius_list.append(R_1)
        x_center_list.append(xc_1)
        z_center_list.append(zc_1)

    if save_data == True:    
        df_sim["x circle center"] = [x_center_list]
        df_sim["z circle center"] = [z_center_list]
        df_sim["circle radius"] = [Radius_list]

        print(df_sim.info())
        df_sim.to_pickle(data_path + df_name)


    return [xc_1 , zc_1 , R_1]


def make_circle(xc,zc,R,ds,xlim,zlim):
    xmin,xmax = xlim
    zmin,zmax = zlim
    x_list,z_list = [],[]
    x,z = xmin, zmax

    steps_size = ds/30
    tol = steps_size #Tolerence for divation of radius
    z_count,x_c = 0,0
    dz = steps_size
    while x <= xmax:
        on_circ = False
        r =  (x-xc)**2 + (z-zc)**2
        if R**2*(1-tol) < r < R**2*(1+tol) :
            x_list.append(x)
            z_list.append(z)
            x += steps_size
            on_circ = True
        
        if on_circ == False:
            z += dz #steps_size
            z_count += 1
            if z_count > zmax/steps_size:
                z = zmax
                z_count = 0
                dz *= -1
                

    return [x_list,z_list]


def make_circle_V2(rc,zc,R,rmax,rmin,zmax,zmin,step_size) ->list:

    steps = int((rmax-rmin)/step_size)
    d= rmax - rmin #diameter
    r_range = np.linspace(-d/2 -d/10,d/2 +d/10,steps)

    r_return, z_return = [], []
    for r in r_range:
        z = np.sqrt( R**2 - r**2 )

        r_return.append(r + rc)
        z_return.append(z + zc)

    for r in r_range:
        z = np.sqrt( R**2 - r**2 )
        r_return.append(r + rc)
        z_return.append(-z + zc)

    return r_return ,z_return


if __name__ == "__main__":
    exit()
