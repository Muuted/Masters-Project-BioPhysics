from Two_D_constants import Two_D_Constants, Two_D_paths
from Two_D_simulation_function import Two_D_simulation_V2, Two_D_simulation_V3
from Make_movie import Make_frames, Make_video
from two_d_data_processing import tot_area, E_pot, E_kin
from Two_D_functions import Langrange_multi, Epsilon_values
from two_d_continues_integration import dSds, descritize_sim_results,find_init_stationary_state
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import os
import progressbar
import scipy
from scipy.special import kv

def test_Lagrange_multi():
    const_args = Two_D_Constants(
        print_val=True
        ,init_rand_psi=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    make_movies = False
    """lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area
        ,psi=psi_list[0]
        ,radi=radi_list[0]
        ,z_list=z_list[0]
        ,print_matrix=True
    )"""
    t = 0
    lambs,nus = Langrange_multi(
        N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
        ,Area=Area
        ,psi=psi_list[t]
        ,radi=radi_list[t]
        ,z_list=z_list[t]
            )

def test_make_frames():
    const_args = Two_D_Constants(
        print_val=True
    )
    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]


    Make_frames(
        data_path=data_path
        ,figs_save_path=video_fig_path
        ,df_name=df_name
    )

def test_make_video():
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    Make_video(
        output_path = video_save_path
        ,input_path = video_fig_path
        ,video_name = df_name
        ,fps=12
    )

def test_epsilon_value():
    const_args = Two_D_Constants(
        print_val=False#True
        #,init_rand_psi=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    sim_steps = 3

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    print(f"r={radi_list[0]}")
    print(f"psi={psi_list[0]}")
    print(f"Area={Area}")
    ef,eg = Epsilon_values(
        N=N,r=radi_list[0],z=z_list[0],psi=psi_list[0],Area=Area
        ,print_matrix=True
    )

    print(f"ef={ef}")
    print(f"eg={eg}")
    
def test_tot_area():

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    dt = df_sim["dt"][0]
    sim_steps = df_sim["sim_steps"][0]

    Area_change= np.zeros(sim_steps)
    time = np.zeros(sim_steps)
    for t in range(sim_steps):
        Area_change[t] += tot_area(N=N,r=r[t],z=z[t])
        time[t] = t*dt

    Amin, Amax = min(Area_change) ,max(Area_change)
    Aratio = Amax/Amin 
    fig, ax=plt.subplots()
    plt.plot(time,Area_change,'.')
    plt.xlabel("time [s]")
    plt.ylabel("total area")
    plt.title(
        f"Ratio of Amax=Amin={Aratio} \n "
        +f"Amax - AMin={Amax-Amin}"
        )
    ax.ticklabel_format(useOffset=False)
    
    dA = np.zeros(sim_steps-1)
    for t in range(sim_steps-1):
        dA[t] = Area_change[t+1] - Area_change[t]

    fig, ax = plt.subplots()
    plt.plot(dA[0:sim_steps-1],'.')
    ax.ticklabel_format(useOffset=False)
    plt.title("Change in Area")
    plt.show()
    
def testing_Epot_Ekin():
    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,video_fig_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    const_args = Two_D_Constants(
        print_val=True
    )
    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    #exit()
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    psi = df_sim['psi'][0]
    Area = df_sim['area list'][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    dt = df_sim["dt"][0]
    sim_steps = df_sim["sim_steps"][0]

    T = []
    S = []

    for t in range(sim_steps-1):
        T.append(
            E_kin(N=N,t=t,dt=dt,r=r,z=z,Area=Area)
            )
        S.append(
            E_pot(N=N,k=k,kG=kG,sigma=sigma,tau=tau,c0=c0
                  ,r=r[t],z=z[t],psi=psi[t],Area=Area)
        )


    plt.figure()
    plt.plot(T)
    plt.title("E kin")

    plt.figure()
    plt.plot(S)
    plt.title("E pot")

    plt.show()

def test_Area_diff_dt():

    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]


    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,figs_for_video_path = path_args[2:4]
    df_name_ref, fps_movie ,num_frames = path_args[4:7]

    dt_list = [dt/2 ,dt ,dt*2]
    Area_data = np.zeros(shape=(3,sim_steps))
    time_vec = np.zeros(shape=(3,sim_steps))
    name_list = [df_name_ref + f" dt={dt}" for dt in dt_list]
    #sim_steps = int(5e2)

    """ To check if all the data files exits as to not
      redo the simulations if the files allready exsits"""
    all_data_sims = False
    redo = False
    file_count = 0
    for i in range(len(name_list)):
        if os.path.isfile(data_path + name_list[i]):
            file_count += 1
    if file_count == len(name_list):
        all_data_sims = True
    

    """Do the simulations to gain the files"""
    if all_data_sims == False or redo == True:
        for i in range(len(dt_list)):
            dt = dt_list[i]
            df_name = name_list[i]
            Two_D_simulation_V3(
                N=N ,k=k ,c0=c0 ,sigma=sigma ,dt=dt ,ds=ds
                ,kG=kG ,tau=tau ,sim_steps=sim_steps
                ,L=L, r0=r0
                ,Area=Area_list
                ,psi=psi_list
                ,radi=radi_list
                ,z_list=z_list
                ,df_name = df_name
                ,num_frames = num_frames
                ,data_path = data_path
                ,Tolerence=1e-4
            )

    """ To plot the data from the simulations"""
    for i in range(len(dt_list)):
        dt = dt_list[i]
        df_name = name_list[i]
        df_sim = pd.read_pickle(data_path + df_name)

        r = df_sim['r'][0]
        z = df_sim['z'][0]
        N = df_sim["N"][0]
        dt = df_sim["dt"][0]
        sim_steps = df_sim["sim_steps"][0]

        
        for t in range(sim_steps):
            Area_data[i][t] = tot_area(N=N,r=r[t],z=z[t])
            time_vec[i][t] = t*dt

    fig, ax = plt.subplots()
    for i in range(len(dt_list)):
        plt.plot(
            time_vec[i],Area_data[i]
            ,".",label=f"dt={dt_list[i]}"
            )
    ax.ticklabel_format(useOffset=False)
    plt.xlabel("time [s]")
    plt.ylabel("Area [arbirary unit]")
    plt.title(label="The total membrane area evolution for different time steps")
    plt.legend(fontsize=15)
    plt.show()

def test_area_correction_difference():
    const_args = Two_D_Constants(
        print_val=True
    )
    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]


    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,figs_for_video_path = path_args[2:4]
    df_name_ref, fps_movie ,num_frames = path_args[4:7]

    df_name_nocorr = df_name_ref + f"No correction"
    df_name_only_corr = df_name_ref + f"only correction"
    name_list = [df_name_nocorr,df_name_only_corr]

    all_data_sims = False
    redo = True
    file_count = 0
    for i in range(len(name_list)):
        if os.path.isfile(data_path + name_list[i]):
            file_count += 1
            print(f"file count={file_count}")
    if file_count == len(name_list):
        all_data_sims = True


    if all_data_sims == False or redo == True:
        Two_D_simulation_V3(
            N=N ,k=k ,c0=c0 ,sigma=sigma ,dt=dt ,ds=ds
            ,kG=kG ,tau=tau ,sim_steps=sim_steps
            ,L=L, r0=r0
            ,Area=Area_list
            ,psi=psi_list
            ,radi=radi_list
            ,z_list=z_list
            ,df_name = df_name_nocorr
            ,num_frames = num_frames
            ,data_path = data_path
            ,Tolerence=1e100
        )

        Two_D_simulation_V3(
            N=N ,k=k ,c0=c0 ,sigma=sigma ,dt=dt ,ds=ds
            ,kG=kG ,tau=tau ,sim_steps=sim_steps
            ,L=L, r0=r0
            ,Area=Area_list
            ,psi=psi_list
            ,radi=radi_list
            ,z_list=z_list
            ,df_name = df_name_only_corr
            ,num_frames = num_frames
            ,data_path = data_path
            ,Tolerence=1e-5
        )


    """ To plot the data from the simulations"""
    Area_data = np.zeros(shape=(2,sim_steps))
    time_vec = np.zeros(shape=(2,sim_steps))
    
    for i in range(len(name_list)):
        df_name = name_list[i]
        df_sim = pd.read_pickle(data_path + df_name)

        r = df_sim['r'][0]
        z = df_sim['z'][0]
        N = df_sim["N"][0]
        dt = df_sim["dt"][0]
        sim_steps = df_sim["sim_steps"][0]

        
        for t in range(sim_steps):
            Area_data[i][t] = tot_area(N=N,r=r[t],z=z[t])
            time_vec[i][t] = t*dt

    fig, ax = plt.subplots()
    for i in range(len(name_list)):
        plt.plot(
            time_vec[i],Area_data[i]
            ,".",label=f"{name_list[i]}"
            )
    
    ax.ticklabel_format(useOffset=False)
    plt.xlabel("time [s]")
    plt.ylabel("Area [arbirary unit]")
    plt.title(
        label="The total membrane area evolution \n"
        +"with and without correction made."
        )
    plt.legend(fontsize=15)
    #plt.show()



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
    print(
        f"nL={n_L}  ,   lambs_L={lambs_L}  \n"
        +f"lc={lc}  ,   k={k}   ,   c0 = {c0}   ,   sigma={sigma} \n"
        +f" p2D ={r_L/lc}"
        )
    #initial values list
    init_conditions = (psi_L ,r_L ,z_L ,n_L ,lambs_L ,nus_L ,A)
    

    #The integration start and stop
    s0 ,sN = 0,(r_L + 10)

    #The integration part.
    ans_odeint = scipy.integrate.solve_ivp(
        dSds,
        t_span = [sN ,s0],
        #t_eval = np.linspace(start=sN,stop=s0,num=5003),
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
    
    r = ans_odeint.y[1]#[0:m]
    z = ans_odeint.y[2]#[0:m]
    psi = ans_odeint.y[0]#[0:m]

    #Then lets load the data from matlab
    matlab_data_path = r"C:\\Users\\AdamSkovbjergKnudsen\\Documents\\GitHub\\Masters-Project-BioPhysics\\Matlab masters project\\saved data\\"    
    matlab_file_name = "Compare integration results.txt"
    df_matlab_data = pd.read_csv(matlab_data_path + matlab_file_name)
    
    m = len(df_matlab_data.loc[0])
    
    r_matlab = df_matlab_data.loc[0][0:m]
    z_matlab = df_matlab_data.loc[4][0:m]
    psi_matlab = df_matlab_data.loc[1]
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
    
    #plt.plot(r_descreet,z_descreet,"*-",label="discreet version",color="k")
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
    

    #r = ans_odeint.y[1]
    plt.figure()
    plt.plot(r,ans_odeint.y[4],".-",label="Python")
    plt.plot(r_matlab,lambs_matlab,".-",label="Matlab")
    plt.hlines(y=tau,xmin=min(r),xmax=max(r),label=r"$\tau$="+f"{tau}")


    plt.title(r"$\lambda$ or tD Lagrange multiplier")
    plt.xlabel("r",fontsize=20)
    plt.ylabel(r"$\lambda$ or tD in matlab",fontsize=20)
    plt.legend(fontsize=20)
    

    plt.figure()
    plt.plot(r,psi,label="Python")
    plt.plot(r_matlab,psi_matlab,label="Matlab")
    plt.title("compare psi in Python and matlab",fontsize=15)
    plt.xlabel("r",fontsize=15)
    plt.ylabel(r"$\psi$",fontsize=15)
    plt.legend()
    plt.show()
    plt.draw()


def test_of_sim_variables_in_stationary_configuration():
    
    c0 = 0.25
    k = 1
    sigma_c = k*c0**2
    tau_c = k*c0
    lc = 1/c0

    sigma = 0.1*sigma_c
    tau = 1*tau_c
    rs2 = 20*lc
    zs2 = 0
    s0, sN = 0, 30*lc
    psi_L = -7.3648e-8
    
    psi,r1,z1,dpsidt,lambs1,nus = find_init_stationary_state(
        sigma=sigma,k=k,c0=c0,tau=tau, ds=0.5
        ,psi_L=psi_L,r_L=rs2,z_L=zs2,s0=s0,sN=sN
        ,total_points=""
    )
    psi,r,z,dpsidt,lambs,nus = find_init_stationary_state(
        sigma=0.1,k=1,c0=1,tau=1, ds = 0.5/4
        ,psi_L=-7.3648e-8,r_L=20,z_L=0,s0=0,sN=30
        ,total_points=""
    )

    plt.figure()
    #plt.plot(r/max(r),z/max(z),"*-",label="dimless")
    plt.plot(r,z,"*-",label="dimless")
    #plt.plot(r[0],z[0],"o")
    #plt.plot(r1/max(r1),z1/max(z1),".-",label="sim units")
    plt.plot(r1,z1,".-",label="sim units")
    #plt.plot(r1[0],z1[0],"o")
    plt.legend()

    plt.figure()
    plt.plot(r,lambs,".-",label="dimless")
    plt.hlines(y=1,xmin=min(r),xmax=max(r),label="dimless")

    plt.plot(r1,lambs1,".-",label="sim units")
    plt.hlines(y=tau,xmin=min(r1),xmax=max(r1),label="sim units")



    plt.legend()
    plt.show()
   

def test_if_constraint_diff_is_correct():
    from Two_D_constants import Two_D_Constants_stationary_state
    from Two_D_functions import Epsilon_values,constraint_f,constraint_g, c_diff_f,c_diff_g

    const_args = Two_D_Constants_stationary_state(
        print_val=True
        ,show_stationary_state=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    
    

    grad_f_r,grad_f_z,grad_f_psi = [[] for i in range(len(psi_list[0]))],[[] for i in range(len(psi_list[0]))],[[] for i in range(len(psi_list[0]))]
    grad_g_r,grad_g_z,grad_g_psi = [[] for i in range(len(psi_list[0]))],[[] for i in range(len(psi_list[0]))],[[] for i in range(len(psi_list[0]))]
    dfdr,dfdz,dfdpsi = [[] for i in range(len(psi_list[0]))],[[] for i in range(len(psi_list[0]))],[[] for i in range(len(psi_list[0]))]
    dgdr,dgdz,dgdpsi = [[] for i in range(len(psi_list[0]))],[[] for i in range(len(psi_list[0]))],[[] for i in range(len(psi_list[0]))]

    """
    As h -> 0  the tolerance for the eror can also decrease and still not get any errors.
    So i think that is evidence that the gradient are correct.
    """
    h = 1e-6
    error_tolerence = 1e-4

    def increase_by_h(dh,x,j):
        var_list = []
        for i in range(len(x)):
            if i == j:
                var_list.append(x[i]+dh)
            else:
                var_list.append(x[i])
        return var_list

    for i in range(len(psi_list[0])):
        for j in range(len(psi_list[0])):
            diff_f_r =(
                    constraint_f(i=i,N=N,psi=psi_list[0],Area=Area_list
                                ,r= increase_by_h(dh=h,x=radi_list[0],j=j)
                                )
                -constraint_f(i=i,N=N,r=radi_list[0],psi=psi_list[0],Area=Area_list)
                )
            grad_f_r[i].append(
                diff_f_r/h
            )

            dfdr[i].append(
                c_diff_f(
                    i=i,j=j,N=N,r=radi_list[0],psi=psi_list[0],Area=Area_list
                    ,diff_var="r")
                    )
            
            diff_f_z =(
                    constraint_f(i=i,N=N,r= radi_list[0]
                    ,psi=psi_list[0],Area=Area_list)
                -constraint_f(i=i,N=N,r=radi_list[0],psi=psi_list[0],Area=Area_list)
                )
            grad_f_z[i].append(
                diff_f_z/h
            )

            dfdz[i].append(
                c_diff_f(
                    i=i,j=j,N=N,r=radi_list[0],psi=psi_list[0],Area=Area_list
                    ,diff_var="z"
                ))
            
            diff_f_psi =(
                    constraint_f(i=i,N=N,r= radi_list[0],Area=Area_list
                    ,psi=increase_by_h(dh=h,x=psi_list[0],j=j)
                    )
                -constraint_f(i=i,N=N,r=radi_list[0],psi=psi_list[0],Area=Area_list)
                )
            grad_f_psi[i].append(
                diff_f_psi/h
            )

            dfdpsi[i].append(
                c_diff_f(
                    i=i,j=j,N=N,r=radi_list[0],psi=psi_list[0],Area=Area_list
                    ,diff_var="psi"
                ))

    for i in range(len(psi_list[0])):
        for j in range(len(psi_list[0])):
            diff_g_r =(
                    constraint_g(i=i,N=N,z=z_list[0],psi=psi_list[0],Area=Area_list
                                ,r=increase_by_h(dh=h,x=radi_list[0],j=j)
                                )
                -constraint_g(i=i,N=N,r=radi_list[0],z=z_list[0],psi=psi_list[0],Area=Area_list)
                )
            grad_g_r[i].append(
                diff_g_r/h
            )
            
            dgdr[i].append(
                c_diff_g(
                    i=i,j=j,N=N,r=radi_list[0],psi=psi_list[0],Area=Area_list,z=z_list[0]
                    ,diff_var="r"
                ))
            
            
            diff_g_z =(
                    constraint_g(i=i,N=N,r=radi_list[0],psi=psi_list[0],Area=Area_list
                                ,z=increase_by_h(dh=h,x=z_list[0],j=j)
                                )
                -constraint_g(i=i,N=N,r=radi_list[0],z=z_list[0],psi=psi_list[0],Area=Area_list)
                )
            grad_g_z[i].append(
                diff_g_z/h
            )

            dgdz[i].append(
                c_diff_g(
                    i=i,j=j,N=N,r=radi_list[0],psi=psi_list[0],Area=Area_list,z=z_list[0]
                    ,diff_var="z"
                ))
            
            
            diff_g_psi =(
                    constraint_g(i=i,N=N,r=radi_list[0],Area=Area_list,z=z_list[0]
                                ,psi=increase_by_h(dh=h,x=psi_list[0],j=j)
                                )
                -constraint_g(i=i,N=N,r=radi_list[0],z=z_list[0],psi=psi_list[0],Area=Area_list)
                )
            
            grad_g_psi[i].append(
                diff_g_psi/h
            )

            dgdpsi[i].append(
                c_diff_g(
                    i=i,j=j,N=N,r=radi_list[0],psi=psi_list[0],Area=Area_list,z=z_list[0]
                    ,diff_var="psi"
                ))


    """ Findig where the gradient are different"""
    error_index_f_r,error_index_f_z,error_index_f_psi = [[],[]],[[],[]],[[],[]]
    error_index_g_r,error_index_g_z,error_index_g_psi = [[],[]],[[],[]],[[],[]]
    
    def find_error_index(
            i,j,error_tol,
            Delta_constraint
            ,diff_constraint
            ):
        null_val =""
        return_i ,return_j = null_val,null_val
        if Delta_constraint[i][j] == 0:
                if  -error_tol < diff_constraint[i][j] and diff_constraint[i][j] < error_tol:
                    return_i = null_val
                    return_j = null_val
                else:
                    return_i = i
                    return_j = j
        else:
            if abs(diff_constraint[i][j]/Delta_constraint[i][j]) < 1 - error_tolerence or abs(diff_constraint[i][j]/Delta_constraint[i][j]) > 1 + error_tolerence:
                return_i = i
                return_j = j
        
        return return_i,return_j, null_val

    for i in range(len(psi_list[0])):
        for j in range(len(psi_list[0])):
            """The f constraint functions"""
            error_i, error_j,null = find_error_index(
                i=i,j=j,error_tol=error_tolerence
                ,Delta_constraint=grad_f_r
                ,diff_constraint=dfdr)
            if error_i != null and error_j != null:
                error_index_f_r[0].append(error_i)
                error_index_f_r[1].append(error_j)


            error_i, error_j,null = find_error_index(
                i=i,j=j,error_tol=error_tolerence
                ,Delta_constraint=grad_f_z
                ,diff_constraint=dfdz)
            if error_i != null and error_j != null:
                error_index_f_z[0].append(error_i)
                error_index_f_z[1].append(error_j)


            error_i, error_j,null = find_error_index(
                i=i,j=j,error_tol=error_tolerence
                ,Delta_constraint=grad_f_psi
                ,diff_constraint=dfdpsi)
            if error_i != null and error_j != null:
                error_index_f_psi[0].append(error_i)
                error_index_f_psi[1].append(error_j)

            """ The g  constraint functions"""
            error_i, error_j,null = find_error_index(
                i=i,j=j,error_tol=error_tolerence
                ,Delta_constraint=grad_g_r
                ,diff_constraint=dgdr)
            if error_i != null and error_j != null:
                error_index_g_r[0].append(error_i)
                error_index_g_r[1].append(error_j)


            error_i, error_j,null = find_error_index(
                i=i,j=j,error_tol=error_tolerence
                ,Delta_constraint=grad_g_z
                ,diff_constraint=dgdz)
            if error_i != null and error_j != null:
                error_index_g_z[0].append(error_i)
                error_index_g_z[1].append(error_j)


            error_i, error_j,null = find_error_index(
                i=i,j=j,error_tol=error_tolerence
                ,Delta_constraint=grad_g_psi
                ,diff_constraint=dgdpsi)
            if error_i != null and error_j != null:
                error_index_g_psi[0].append(error_i)
                error_index_g_psi[1].append(error_j)
    
    print(
        f"error index f_r ={error_index_f_r} \n"
        +f"error index f_z ={error_index_f_z} \n"
        +f"error index f_psi ={error_index_f_psi} \n"
        +f"\n"
        f"error index g_r ={error_index_g_r} \n"
        +f"error index g_z ={error_index_g_z} \n"
        +f"error index g_psi ={error_index_g_psi} "
    )
    exit()
    for i in range(N):
        fig, ax = plt.subplots(nrows=3,ncols=2)
        
        #manager = plt.get_current_fig_manager()
        #manager.resize(*manager.window.maxsize())

        plt.get_current_fig_manager().full_screen_toggle()

        #fig.set_size_inches(w=10,h=10)
        font_size= 15
        ax[0,0].set_title(
            r"dfdr =$\frac{\partial f_{i} }{\partial q_{i}}$ and diff_$f_{i}$ =$\frac{ f_{i}(q_{i} +h) - f_{i}(q_{i})}{h}$"
            +f" \n i={i}"
            ,fontsize=font_size
            )
        ax[0,0].plot(grad_f_r[i],"o-",label="diff_f")
        ax[0,0].plot(dfdr[i],".-",label="dfdr")
        ax[0,0].legend()
        
        ax[1,0].plot(grad_f_z[i],"o-",label="diff_f")
        ax[1,0].plot(dfdz[i],".-",label="dfdz")
        ax[1,0].legend()

        ax[2,0].plot(grad_f_psi[i],"o-",label="diff_f")
        ax[2,0].plot(dfdpsi[i],".-",label=r"dfd$\psi$")
        ax[2,0].set_xlabel(r"$i$",fontsize=font_size)
        ax[2,0].legend()


        ax[0,1].set_title(
            r"dgdr=$\frac{\partial g_{i}}{\partial q_{j}}$ and diff_$g_{i}=\frac{ g_{i} (q_{j}+h) - g_{i} (q_{j})}{h}$"
            +f" \n i={i}"
            ,fontsize=font_size
            )
        ax[0,1].plot(grad_g_r[i],"o-",label="diff_g")
        ax[0,1].plot(dgdr[i],".-",label="dgdr")
        ax[0,1].legend()

        ax[1,1].plot(grad_g_z[i],"o-",label="diff_g")
        ax[1,1].plot(dgdz[i],".-",label="dgdz")
        ax[1,1].legend()

        ax[2,1].plot(grad_g_psi[i],"o-",label="diff_g")
        ax[2,1].plot(dgdpsi[i],".-",label=r"dgd$\psi$")
        ax[2,1].set_xlabel(r"$i$",fontsize=font_size)
        ax[2,1].legend()

        #plt.show()
        plt.draw()
        plt.pause(5)
        fig.clear()


def testing_new_epsilon_matrix():
    from Two_D_constants import Two_D_Constants_stationary_state
    from Two_D_functions import Epsilon_values,constraint_f,constraint_g, c_diff_f,c_diff_g
    from Testing_ideas import Epsilon_v2
    const_args = Two_D_Constants_stationary_state(
        print_val=True
        ,show_stationary_state=True
    )
    print("\n \n \n")
    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]
    
    t = 0
    print("old epsilon function")
    ebf_old, ebg_old, A_original,b_original = Epsilon_values(
        N=N, r=radi_list[t], z=z_list[t] ,psi=psi_list[t] ,Area=Area_list
        ,print_matrix=False
        ,testing= True
                )
    
    #print(f"ebf={ebf} \n ebg={ebg}")
    print("new epsilon function")
    ebf_new,ebg_new, A_new, b_new = Epsilon_v2(
        N=N, r=radi_list[t], z=z_list[t] ,psi=psi_list[t] ,Area=Area_list
        ,print_matrix=False
        ,testing= True
    )
    #print(f"ebf={ebf} \n ebg={ebg}")

    Diff_A = A_original - A_new
    print(f"diff = {ebf_old - ebf_new}")
    print(f"diff = {ebg_old - ebg_new}")
    print(f"diff ={A_original-A_new}")

    print(f"sum of differenc of A = {np.sum(A_original-A_new)}")


    print("Using direct calculations to verify the materix")
    A11 = 4*np.pi**2*( radi_list[0][1]**2 + radi_list[0][0]**2)/Area_list[0] + (np.sin(psi_list[0][0]))**2

    print("A_{11} ="+f"{A11}" + f" and for the matrix value we have A[0,0]={A_new[0,0]}")

    print(
        f"radi={radi_list[0]} \n"
        +f"z={z_list[0]} \n"
        +f"psi={psi_list[0]}"
    )


def Testing_total_area_function():
    from two_d_data_processing import tot_area
    from Two_D_constants import Two_D_Constants_stationary_state
    flat = False

    if flat  == True:
        const_args = Two_D_Constants(
            print_val=True
        )

        L,r0,N,ds,T,dt = const_args[0:6]
        k,c0,sim_steps = const_args[6:9]
        sigma, tau, kG = const_args[9:12]
        Area, psi_list = const_args[12:14]
        radi_list,z_list = const_args[14:16]

        """ first lets test if the Area is the same"""
        

        Area_from_list = np.sum(Area)
        Area_from_function = tot_area(N=N,r=radi_list[0],z=z_list[0])

        print(
            f"Area_from_list ={Area_from_list} \n"
            +f"Area_from_function ={Area_from_function} \n"
            +f"difference = {Area_from_list-Area_from_function}"
        )
    
    if flat == False:
        const_args = Two_D_Constants_stationary_state(
        print_val=False
        ,show_stationary_state=False
        )

        L,r0,N,ds,T,dt = const_args[0:6]
        k,c0,sim_steps = const_args[6:9]
        sigma, tau, kG = const_args[9:12]
        Area_list, psi_list = const_args[12:14]
        radi_list,z_list = const_args[14:16]

        """ first lets test if the Area is the same"""
        Area_from_list = np.sum(Area_list)
        Area_from_function = tot_area(N=N,r=radi_list[0],z=z_list[0])

        print(
            f"Area_from_list ={Area_from_list} \n"
            +f"Area_from_function ={Area_from_function} \n"
            +f"difference = {Area_from_list-Area_from_function}"
        )



def testing_perturbation_function():
    from Two_D_functions import Perturbation_of_inital_state
    from Two_D_constants import Two_D_Constants_stationary_state

    const_args = Two_D_Constants_stationary_state(
        print_val=True
        ,show_stationary_state=False
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]


    r_init,z_init,psi_init =[i for i in radi_list[0]],[i for i in z_list[0]],[i for i in psi_list[0]]
    
    Perturbation_of_inital_state(
        points_perturbed=3, ds=ds
        ,r=r_init,z=z_init,psi=psi_init
        ,delta_psi= -1e-6
    )

    plt.figure()
    font_size=15
    plt.plot(radi_list[0],z_list[0],"o-",label="Unperturbed")
    plt.plot(r_init,z_init,"o-",label="Perturbed")
    plt.plot(radi_list[0][0],z_list[0][0],"o",label="start point i=0",color="k")
    plt.xlabel("r",fontsize=font_size)
    plt.ylabel("z",fontsize=font_size)
    plt.title("Comparing the Unperturbed state and the Perturbed state",fontsize=font_size)
    plt.legend()



    plt.figure()
    plt.plot(radi_list[0][0:N],psi_list[0],"o-",label="Unperturbed")
    #plt.plot(r_init[0:N],psi_init,".-",label="r,psi")
    plt.plot(radi_list[0][0:N],psi_init,".-",label="perturbed")
    plt.title(r"Pertubation of $\psi$ at the end",fontsize=font_size)
    plt.xlabel("r",fontsize=font_size)
    plt.ylabel(r"$\psi$",fontsize=font_size)
    plt.legend()
    plt.show()


def testing_values_in_epsilon():
    from Two_D_constants import Two_D_Constants_stationary_state
    from Two_D_functions import Epsilon_values,constraint_f,constraint_g, c_diff_f,c_diff_g, Epsilon_v2
    const_args = Two_D_Constants_stationary_state(
        print_val=False
        ,show_stationary_state=False
    )
    print("\n \n \n")
    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]
    
    t = 0    
    #print(f"ebf={ebf} \n ebg={ebg}")
    print("new epsilon function")
    ebf_new,ebg_new, A_new, b_new = Epsilon_v2(
        N=N, r=radi_list[t], z=z_list[t] ,psi=psi_list[t] ,Area=Area_list
        ,print_matrix=False
        ,testing= True
    )


    print("Using direct calculations to verify the materix")
    A11 = 4*np.pi**2*( radi_list[0][1]**2 + radi_list[0][0]**2)/Area_list[0]**2 + (np.sin(psi_list[0][0]))**2
    A12 = -4*np.pi**2*(radi_list[0][1]**2)/(Area_list[0]*Area_list[1])
    A13 = 2*np.pi**2*(z_list[0][1]-z_list[0][0])*(radi_list[0][1]-radi_list[0][0])/Area_list[0]**2 - np.cos(psi_list[0][0])*np.sin(psi_list[0][0])
    zz = (z_list[0][2] - z_list[0][1])**2
    rr = (radi_list[0][2] + radi_list[0][1])**2 
    A44 = np.pi**2*( zz + rr)/Area_list[1]**2 + np.cos(psi_list[0][1])**2
    print(
        "\n \n A_{11} ="+f"{A11}" + f" and for the matrix value we have A[0,0]={A_new[0,0]} \n"
        +f" A[0,0]-A_11 ={A_new[0,0]-A11} \n \n"
        +f"A_12={A12}  and A[0,1]={A_new[0,1]} \n A[0,1]-A_12= {A_new[0,1]-A12} \n \n"
        +f"A_12={A13}  and A[0,2]={A_new[0,2]} \n A[0,1]-A_12= {A_new[0,2]-A13} \n \n"
        +f"A_44={A44}  and A[3,3]={A_new[3,3]} \n A[3,3]-A_44= {A_new[3,3]-A44} \n \n"
        )
    


def testing_initial_angles():
    from Two_D_constants import Two_D_Constants_stationary_state
    print("\n \n \n")
    const_args = Two_D_Constants_stationary_state(
        print_val=False
        ,show_stationary_state=True
    )
    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    c = ds*1.5

    plt.figure()
    plt.plot(radi_list[0],z_list[0],"o-",label="initial state")

    x = [
        radi_list[0][1]
        , radi_list[0][1] + c*np.cos(psi_list[0][0])
        ]
    y = [
        z_list[0][1]
        , z_list[0][1] + c*np.sin(psi_list[0][0])
        ]

    plt.plot(x,y,color="k",label="angle test")
    plt.plot(x[0],y[0],color="g",marker=".")

    plt.legend()
    for i in range(1,N):
        j = i 
        x = [radi_list[0][j+1], radi_list[0][j+1] + c*np.cos(psi_list[0][j])]
        y = [z_list[0][j+1], z_list[0][j+1] + c*np.sin(psi_list[0][j])]
        plt.plot(x,y,color="k",label="angle test")
        plt.plot(x[0],y[0],color="g",marker=".")
    
    
    #plt.figure()
    #plt.plot(radi_list[0][0:N],psi_list[0])
    plt.show()


def testing_if_constraints_are_true():
    from Two_D_constants import Two_D_Constants_stationary_state
    from Two_D_functions import Epsilon_values,constraint_f,constraint_g, c_diff_f,c_diff_g, Epsilon_v2
    print("\n \n \n")
    const_args = Two_D_Constants_stationary_state(
        print_val=False
        ,show_stationary_state=True
    )
    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    #psi_list[0] = [i + np.pi for i in psi_list[0]]
    Tol= 1e-5

    f_cons, i_count_f = [],[]
    g_cons ,i_count_g = [], []
    for i in range(N):
        f = constraint_f(i=i,N=N,r=radi_list[0],psi=psi_list[0],Area=Area_list)
        if f < - Tol or f > Tol :
            i_count_f.append(i)
            f_cons.append(f)

        g = constraint_g(i=i,N=N,r=radi_list[0],z=z_list[0],psi=psi_list[0],Area=Area_list)
        if g < - Tol or g > Tol :
            i_count_g.append(i)
            g_cons.append(g)
    

    plt.figure()
    plt.plot(i_count_f,f_cons,".-",label="f")
    plt.plot(i_count_g,g_cons,".-",label="g")
    plt.title(
        f"N={N}, len(f)={len(f_cons)} and len(g)={len(g_cons)}"
    )
    plt.legend()
    plt.show()

def testing_arctan2_function():
    from two_d_continues_integration import Get_angle
    theta = [0, round(np.pi/4,2), round(np.pi/2,2) , round(3*np.pi/4,2), round(np.pi,2)]
    print(f"Possible  \n angles = [0, np.pi/4, np.pi/2 ,3*np.pi/4, np.pi] \n angles = {theta}")

    def draw_circle(r,a=0,b=0):
        x = []
        dr = 0.01
        y_pos, y_neg = [],[]
        for i in np.arange(-r,r+dr,dr):
            x.append(i)
            #print(f"np.sqrt( r**2 - (i-a)**2) and i = {i}")
            if i == -r or i == r:
                y_pos.append(0)
                y_neg.append(0)
            elif i < -r or i > r:
                y_pos.append(0)
                y_neg.append(0)
            else:
                y = abs(np.sqrt( r**2 - (i-a)**2) + b)
                y_pos.append(y)
                y_neg.append(-y)

        return x, y_pos,y_neg
    
    c = 1
    ang = 3*np.pi/4
    p1 = [1,1]
    p1 = [0,0]
    p2 = [p1[0] + c*np.cos(ang) ,   p1[1] + np.sin(ang)]
    

    #fig,ax = plt.subplots(1,2)
    plt.figure()

    plt.plot([p1[0], p2[0]], [p1[1], p2[1]], "o-",label="start points")

    x = [p1[0] - p1[0]  , p2[0] - p1[0]  ]
    y = [p1[1] - p1[0]  ,   p2[1] - p1[1]]
    plt.plot(x,y, "o-",label="moved points")

    x_circ, y_pos_circ,y_neg_circ = draw_circle(r=0.5)

    plt.plot(x_circ,y_pos_circ,"k")
    plt.plot(x_circ,y_neg_circ,"k")
    theta = np.arctan2(y[1],x[1])
    psi = Get_angle(
        x1=p1[0],y1=p1[1]
        ,x2=p2[0],y2=p2[1]
    )
    print(f"theta = {round(theta,2)}")
    print(f"psi = {round(psi,2)}")
    print(f"that-pi = {round(theta-np.pi,2)}")

    plt.legend()
    plt.grid()
    plt.xlim(-1.2,1.2)
    plt.ylim(-1.2,1.2)
    plt.show()

if __name__ == "__main__":
    #test_Lagrange_multi()
    #test_make_frames()
    #test_make_video()
    #test_check_area()
    #test_epsilon_value()
    #test_tot_area()
    #testing_Epot_Ekin()
    #test_Area_diff_dt()
    #test_area_correction_difference()
    
    #Test_with_matlab_integrate_solution()
    #test_of_sim_variables_in_stationary_configuration()

    """from Two_D_constants import Two_D_Constants_stationary_state
    Two_D_Constants_stationary_state(
        show_stationary_state=False
        ,print_val=True
    )"""

    #test_if_constraint_diff_is_correct()
    #testing_new_epsilon_matrix()
    #testing_perturbation_function()
    #Testing_total_area_function()
    #testing_values_in_epsilon()
    testing_initial_angles()
    testing_if_constraints_are_true()
    #testing_arctan2_function()
