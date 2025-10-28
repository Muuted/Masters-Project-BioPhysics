#from Two_D_constants import Two_D_Constants
import numpy as np
import scipy
from scipy.special import kv
import matplotlib.pyplot as plt
import pandas as pd


def dSds(s,S,k,sigma,c0,tau,kG):#,k,sigma,c0):
    psi, r, z,n,lambs,nus,A = S
    #k,kG,sigma,c0 = p
    drds = np.cos(psi)
    dzds = np.sin(psi)
    dlambs_ds = (k/2)*( (n-c0)**2 - np.sin(psi)**2/r**2) + sigma
    dnu_ds = 0

    dpsids  = n 
    dnds = np.sin(psi)*np.cos(psi)/r**2 - n*np.cos(psi)/r + lambs*np.sin(psi)/(k*r) #- nus*np.cos(psi)/(k*r)

    dAds = 2*np.pi*r

    return [dpsids ,drds ,dzds ,dnds ,dlambs_ds ,dnu_ds ,dAds]


def descritize_sim_results(r,z,ds,max_num_points=""):
    j=len(r)-1
    index_list = [j]
    for i in range(len(r)-1,0,-1):
        dist = np.sqrt( (r[j]-r[i])**2 + (z[j]-z[i])**2)
        if dist >= ds:
            #r_return.append(r[i-1])
            index_list.append(i-1)
            j = i-1
            if max_num_points != "":
                if len(index_list) > max_num_points:
                    break
    return index_list

def Get_angle(x1,y1,x2,y2):
    """ x1,y1 is closer to the free edge along the membrane"""
    x1 = [x1 - x2]
    y1 = [y1 - y2]
    psi = np.pi - np.arctan2(y1,x1)[0]
    return -psi

def edge_tension(t,y,k,sigma,c0,tau,kG):
        return tau - y[4]
    
def edge_ratio(t,y,k,sigma,c0,tau,kG):
    dpsidt_1 = y[3]
    psi_1 = y[0]
    r_1 = y[1]
    alpha = kG/k
    val = (1-dpsidt_1)*r_1/np.sin(psi_1)-1 + alpha
    return val


def find_init_stationary_state(
        sigma,k,c0,tau,psi_L,r_L,z_L,s0,sN,ds,kG
        ,total_points = ""
        ,find_edge_condition = True
        ):

    #initial values
    lc = 1/np.sqrt(c0**2/2 + sigma/k) # characterisitic length in the aymptotic regime.
    n_L = (psi_L/kv(1,r_L/lc))*( -kv(0,r_L/lc) - kv(1,r_L/lc)/(r_L/lc))/lc
    lambs_L = (k*c0**2/2 + sigma )*r_L
    nus_L = 0 # nu(s_1) = nu(s_2) = 0 from that we know the outer value
    A = 0 #2*np.pi*( r_L**2 - r0**2  )

    #constants
    #args_list = (k ,sigma ,c0)
    args_list = (k,sigma,c0,tau,kG)
    #initial values list
    init_conditions = (psi_L ,r_L ,z_L ,n_L ,lambs_L ,nus_L ,A)
    #The integration part.
    
    edge_tension.terminal = True
    ans_odeint = scipy.integrate.solve_ivp(
        dSds
        ,t_span = [sN ,s0]
        ,t_eval = np.linspace(start=sN,stop=s0,num=10000)
        ,y0 = init_conditions
        ,args = args_list
        ,method="LSODA"
        ,rtol=1e-10
        ,atol=1e-20
        ,events=(edge_tension,edge_ratio)
    )
    psi = ans_odeint.y[0]
    r = ans_odeint.y[1]
    z = ans_odeint.y[2]
    dpsidt = ans_odeint.y[3]
    lambs = ans_odeint.y[4]
    nus = ans_odeint.y[5]
    
    dpsidt_1 = dpsidt[len(r)-1]
    psi_1 = psi[len(r)-1]
    r_1 = r[len(r)-1]
    alpha = kG/k
    
    index_list = descritize_sim_results(
        r = ans_odeint.y[1]#[0:m]
        ,z = ans_odeint.y[2]#[0:m]
        ,ds = ds
        ,max_num_points = total_points
        )

    r_discrete,z_discrete = [],[]
    psi_discrete,dpsidt_discrete = [],[]
    lambs_discrete, nus_discrete = [],[]
    for i in index_list:
        r_discrete.append(ans_odeint.y[1][i])
        z_discrete.append(ans_odeint.y[2][i])
        dpsidt_discrete.append(ans_odeint.y[3][i])
        lambs_discrete.append(ans_odeint.y[4][i])
        nus_discrete.append(ans_odeint.y[5][i])

    for i in range(len(r_discrete)-1):
        psi_discrete.append(
            Get_angle(
                x1=r_discrete[i] ,y1=z_discrete[i]
                ,x2=r_discrete[i+1] ,y2=z_discrete[i+1]
                )
        )
    
    r = ans_odeint.y[1]#[0:m]
    z = ans_odeint.y[2]#[0:m]
    return psi_discrete,r_discrete,z_discrete, r,z
    
if __name__ == "__main__":
    #integrate_solution()
    #Test_with_matlab_integrate_solution()
    pass