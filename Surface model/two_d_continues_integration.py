from Two_D_constants import Two_D_Constants
import numpy as np
import scipy
from scipy.special import kv
import matplotlib.pyplot as plt
import pandas as pd


def dSds(s,S,k,sigma,c0):
    psi, r, z,n,lambs,nus,A = S
    #k,kG,sigma,c0 = p
    drds = np.cos(psi)
    dzds = np.sin(psi)
    dlambs_ds = (k/2)*( (n-c0)**2 -np.sin(psi)**2/r**2) + sigma
    dnu_ds = 0

    dpsids  = n 
    dnds = np.sin(psi)*np.cos(psi)/r**2 - n*np.cos(psi)/r + lambs*np.sin(psi)/r #- nus*np.cos(psi)/(k*r)

    dAds = 2*np.pi*r

    return [dpsids ,drds ,dzds ,dnds ,dlambs_ds ,dnu_ds ,dAds]

def integrate_solution():
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    
    
    #args list
    sigma_c = k*c0**2
    k_c = k
    kG = -0.75*k
    tau_c = k*c0
    c0_c = 0.25
    lc = 1/np.sqrt(0.5 + sigma_c)

    # making the dimless variables.
    c0 = c0/c0_c
    tau = 1 #tau/tau_c
    sigma = 0.1 #sigma/sigma_c
    k = k/k_c

    args_list = (k ,kG ,sigma ,c0)
    
    #initial values 
    psi_L = -7.3648e-8
    r_L =  20 #(tau*lc**2/k)*1.01
    z_L = 0
    n_L = (psi_L/kv(1,r_L/lc))*( -kv(0,r_L/lc) - kv(1,r_L/lc)/(r_L/lc)    )
    lambs_L = (k*c0**2/2 + sigma)*r_L
    nus_L = 0 # nu(s_1) = nu(s_2) = 0 from that we know the outer value
    A = 0#2*np.pi*( r_L**2 - r0**2  )
    
    #initial values list
    init_conditions = (psi_L ,r_L ,z_L ,n_L ,lambs_L ,nus_L ,A)
    

    #The integration start and stop
    s0 ,sN = 0 ,r_L + 10

    #The integration part.
    ans_odeint = scipy.integrate.solve_ivp(
        dSds,
        t_span = [sN ,s0],
        t_eval = np.linspace(start=sN,stop=s0,num=100),
        y0 = init_conditions,
        args = args_list
    )

    print(f"y =[r ,z ,psi ,dpsids ,lambda ,nu ,A]")
    print(ans_odeint)

    r = ans_odeint.y[1]
    z = ans_odeint.y[2]

    m = 0
    for i in range(len(ans_odeint.y[5])-1):
        if ans_odeint.y[5][i] < tau:
            #print(f"==tau found and i={i}")
            m = i
            #break
    print(f"m={m}")
    r = ans_odeint.y[1][0:m]
    z = ans_odeint.y[2][0:m]
    fig, ax = plt.subplots()
    plt.plot(r,z,"-")
    plt.xlabel("r")
    plt.ylabel("z")
    plt.show()

def descritize_sim_results(r,z,ds,max_num_points=""):
    j=len(r)-1
    index_list = [j]
    for i in range(len(r)-2,0,-1):
        dist = np.sqrt( (r[j]-r[i])**2 + (z[j]-z[i])**2)
        if dist > ds:
            #r_return.append(r[i-1])
            index_list.append(i-1)
            j = i-1
            if max_num_points != "":
                if len(index_list) > max_num_points:
                    break
            
    #print(f"index_list={index_list}")
    return index_list

def Test_with_matlab_integrate_solution():

    const_args = Two_D_Constants(
        print_val=False
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    
    
    #args list
    sigma_c = k*c0**2
    k_c = k
    kG = -0.75*k
    tau_c = k*c0
    c0_c = c0

    # making the dimless variables.
    c0 = c0/c0_c
    tau = 1#tau/tau_c #1
    sigma = 0.1#sigma/sigma_c #0.1
    k = k/k_c
    #kG = -0.75*k

    args_list = (k ,sigma ,c0)
    
    #initial values
    lc = 1/np.sqrt(0.5 + sigma) # characterisitic length in the aymptotic regime.
    psi_L = -7.3648e-8 
    r_L =  20.0 #(tau*lc**2/k)*1.01
    z_L = 0.0
    n_L = (psi_L/kv(1,r_L/lc))*( -kv(0,r_L/lc) - kv(1,r_L/lc)/(r_L/lc))/lc# 5.8973e-08 #
    
    lambs_L = (k*c0**2/2 + sigma)*r_L
    print(lambs_L)
    nus_L = 0 # nu(s_1) = nu(s_2) = 0 from that we know the outer value
    A = 0#2*np.pi*( r_L**2 - r0**2  )
    """print(
        f"nL={n_L}  ,   lambs_L={lambs_L}  \n"
        +f"lc={lc}  ,   k={k}   ,   c0 = {c0}   ,   sigma={sigma} \n"
        +f" p2D ={r_L/lc}"
        )"""
    #initial values list
    init_conditions = (psi_L ,r_L ,z_L ,n_L ,lambs_L ,nus_L ,A)
    

    #The integration start and stop
    s0 ,sN = 0,(r_L + 10)

    #The integration part.
    ans_odeint = scipy.integrate.solve_ivp(
        dSds,
        t_span = [sN ,s0],
        t_eval = np.linspace(start=sN,stop=s0,num=5003),
        y0 = init_conditions,
        args = args_list,
        method="LSODA"#"RK45"
    )

    #print(f"y =[r ,z ,psi ,dpsids ,lambda ,nu ,A]")
    #print(ans_odeint)
    tol = 1e-2
    m = len(ans_odeint.y[0])
    for i in range(len(ans_odeint.y[0])-1):
        #if tau*(1-tol) <ans_odeint.y[4][i] < tau*(1.+tol):
        if ans_odeint.y[4][i] < tau:
            m = i
            break
    print(f"m=",m)
    
    r = ans_odeint.y[1][0:m]
    z = ans_odeint.y[2][0:m]

    """Then lets load the data from matlab"""
    matlab_data_path = r"C:\\Users\\AdamSkovbjergKnudsen\\Documents\\GitHub\\Masters-Project-BioPhysics\\Matlab masters project\\saved data\\"    
    matlab_file_name = "Compare integration results.txt"
    df_matlab_data = pd.read_csv(matlab_data_path + matlab_file_name)
    
    m = len(df_matlab_data.loc[0])
    
    r_matlab = df_matlab_data.loc[0][0:m]
    z_matlab = df_matlab_data.loc[4][0:m]
    lambs_matlab = df_matlab_data.loc[3]
    
    fig, ax = plt.subplots()
    plt.plot(r,z,".-"
             ,label="Python results"
             )
    plt.plot(r_matlab,z_matlab,".-",label="matlab results")
    index_list = descritize_sim_results(r=r,z=z,ds=0.5,max_num_points=10)
    r_descreet,z_descreet = [],[]
    for i in index_list:
        r_descreet.append(r[i])
        z_descreet.append(z[i])
    
    plt.plot(r_descreet,z_descreet,"*-",label="discreet version",color="k")
    plt.xlabel("r",fontsize=20)
    plt.ylabel("z",fontsize=20)
    plt.title(
        f"The constants values: \n"
        r"$\tau$"+f"={tau} ,"
        +r"$\sigma$"+f"={sigma} ,"
        +r"$k$"+f"={k} ,"
        +r"$c_0$"+f"={c0}"
        )
    plt.legend(fontsize=20)
    

    r = ans_odeint.y[1]
    plt.figure()
    plt.plot(r,ans_odeint.y[4],".-",label="Python")
    plt.plot(r_matlab,lambs_matlab,".-",label="Matlab")
    plt.hlines(y=tau,xmin=min(r),xmax=max(r),label=r"$\tau$="+f"{tau}")


    plt.title(r"$\lambda$ or tD Lagrange multiplier")
    plt.xlabel("r",fontsize=20)
    plt.ylabel(r"$\lambda$ or tD in matlab",fontsize=20)
    plt.legend(fontsize=20)
    

    plt.show()

def find_init_stationary_state(sigma,k,c0,tau,psi_L,r_L,z_L,s0,sN,total_points=""):

    #initial values
    lc = 1/np.sqrt(0.5 + sigma) # characterisitic length in the aymptotic regime.
    n_L = (psi_L/kv(1,r_L/lc))*( -kv(0,r_L/lc) - kv(1,r_L/lc)/(r_L/lc))/lc# 5.8973e-08 #
    lambs_L = (k*c0**2/2 + sigma)*r_L
    nus_L = 0 # nu(s_1) = nu(s_2) = 0 from that we know the outer value
    A = 0#2*np.pi*( r_L**2 - r0**2  )


    #constants
    args_list = (k ,sigma ,c0)
    #initial values list
    init_conditions = (psi_L ,r_L ,z_L ,n_L ,lambs_L ,nus_L ,A)
    #The integration part.
    ans_odeint = scipy.integrate.solve_ivp(
        dSds,
        t_span = [sN ,s0],
        t_eval = np.linspace(start=sN,stop=s0,num=10003),
        y0 = init_conditions,
        args = args_list,
        method="LSODA"
    )

    m = len(ans_odeint.y[0])-1
    for i in range(len(ans_odeint.y[0])-1):
        if ans_odeint.y[4][i] < tau: # lambda = ans_odeint.y[4]
            m = i
            break
    
    index_list = descritize_sim_results(
        r=ans_odeint.y[1][0:m],z=ans_odeint.y[2][0:m]
        ,ds=0.5,max_num_points=total_points
        )
    r_discrete,z_discrete = [],[]
    psi_discrete,dpsidt_discrete = [],[]
    lambs_discrete, nus_discrete = [],[]
    for i in index_list:
        psi_discrete.append(ans_odeint.y[0][i])
        r_discrete.append(ans_odeint.y[1][i])
        z_discrete.append(ans_odeint.y[2][i])
        dpsidt_discrete.append(ans_odeint.y[3][i])
        lambs_discrete.append(ans_odeint.y[4][i])
        nus_discrete.append(ans_odeint.y[5][i])
    
    
    return psi_discrete,r_discrete,z_discrete,dpsidt_discrete,lambs_discrete,nus_discrete
    
if __name__ == "__main__":
    #integrate_solution()
    #Test_with_matlab_integrate_solution()
    c0 = 0.25
    k = 1
    sigma_c = k*c0**2
    tau_c = k*c0
    lc = 1/c0
    sigma = 0.1*sigma_c
    lc_asymp = 1/np.sqrt(0.5 + sigma_c)
    
    tau = 1/(k*c0)
    rs2 = 20
    zs2 = 0
    s0, sN = 0, 30/lc
    psi_L = -7.6348e-8

    psi,r1,z1,dpsidt,lambs1,nus = find_init_stationary_state(
        sigma=sigma,k=k,c0=c0,tau=tau
        ,psi_L=psi_L,r_L=rs2,z_L=zs2,s0=s0,sN=sN
        ,total_points=""
    )
    psi,r,z,dpsidt,lambs,nus = find_init_stationary_state(
        sigma=0.1,k=1,c0=1,tau=1
        ,psi_L=-7.6348e-8,r_L=20,z_L=0,s0=0,sN=30
        ,total_points=""
    )
    plt.figure()
    plt.plot(r,z,".-",label="sim units")
    #plt.plot(r[0],z[0],"o")
    plt.plot(r1,z1,".-",label="dimless")
    #plt.plot(r1[0],z1[0],"o")
    plt.legend()

    plt.figure()
    plt.plot(r1,lambs1,".-",label="sim units")
    plt.hlines(y=tau,xmin=min(r1),xmax=max(r1),label="sim units")

    plt.plot(r,lambs,".-",label="dimless")
    plt.hlines(y=1,xmin=min(r),xmax=max(r),label="dimless")
    plt.legend()
    plt.show()
   