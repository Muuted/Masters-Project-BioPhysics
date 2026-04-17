#from Two_D_constants import Two_D_Constants#, Two_D_paths
#from Two_D_simulation_function import Two_D_simulation_V2, Two_D_simulation_V3
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
np.set_printoptions(legacy='1.25')

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
    from Two_D_constants import Two_D_Constants_stationary_state

    const_args = Two_D_Constants_stationary_state(
        print_val=True
        ,show_stationary_state=True
        ,start_flat=False
        ,perturb=False
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]
    r_unperturbed, z_unperturbed = const_args[16:18]
    eta = const_args[18]
    
    
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

    args_list = (k,sigma,c0,tau,kG)#(k ,sigma ,c0)
    
    #initial values
    lc = 1/np.sqrt(0.5 + sigma) # characterisitic length in the aymptotic regime.
    psi_L = -7.3648e-8 
    r_L =  20.0 #(tau*lc**2/k)*1.01
    z_L = 0.0
    n_L = (psi_L/kv(1,r_L/lc))*( -kv(0,r_L/lc) - kv(1,r_L/lc)/(r_L/lc))/lc# 5.8973e-08 #
    
    lambs_L = (k*c0**2/2 + sigma)*r_L
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
        dSds
        ,t_span = [sN ,s0]
        ,t_eval = np.linspace(start=sN,stop=s0,num=5003)
        ,y0 = init_conditions
        ,args = args_list
        ,method="LSODA"#"RK45",
        ,atol=1e-20
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
    
    r1 = ans_odeint.y[1][0:m]
    z1 = ans_odeint.y[2][0:m]
    psi1 = ans_odeint.y[0][0:m]

    plt.figure()
    plt.plot(r1,z1,label="membrane")
    plt.title(r"with stopping condition $\lambda( s_1 ) = \tau $ with atol")
    plt.xlabel("r",fontsize=15)
    plt.ylabel("z",fontsize=15)
    plt.legend()

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


    save_path = "C:\\Users\\AdamSkovbjergKnudsen\\Desktop\\skole\\1 Tidligere semestre\\Kandidat speciale\\Figures\\"

    fig, ax = plt.subplots()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
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
    #plt.savefig(save_path + "")
    

    #r = ans_odeint.y[1]
    plt.figure()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.plot(r,ans_odeint.y[4],".-",label="Python")
    plt.plot(r_matlab,lambs_matlab,".-",label="Matlab")
    plt.hlines(y=tau,xmin=min(r),xmax=max(r),label=r"$\tau$="+f"{tau}")

    plt.title(r"$\lambda$ or tD Lagrange multiplier with atol")
    plt.xlabel("r",fontsize=20)
    plt.ylabel(r"$\lambda$ or tD in matlab",fontsize=20)
    plt.legend(fontsize=20)
    plt.savefig(save_path + "constraint_multi_lambda_comparision_tau=1.png")
    plt.savefig(save_path + "constraint_multi_lambda_comparision_tau=1.svg")

    
    plt.figure()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.plot(r,ans_odeint.y[4],".-",label="Python")
    plt.plot(r_matlab,lambs_matlab,".-",label="Matlab")
    plt.hlines(y=tau,xmin=min(r),xmax=max(r),label=r"$\tau$="+f"{tau}")

    plt.title(r"$\lambda$ or tD Lagrange multiplier with atol")
    plt.xlabel("r",fontsize=20)
    plt.ylabel(r"$\lambda$ or tD in matlab",fontsize=20)
    plt.xlim(0.1,2.6)
    plt.ylim(-1,1.6)
    plt.legend(fontsize=15,loc="lower right")
    plt.savefig(save_path + "constraint_multi_lambda_comparision_tau=1_zoomed.png")
    plt.savefig(save_path + "constraint_multi_lambda_comparision_tau=1_zoomed.svg")


    plt.figure()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.plot(r,psi,".-",label="Python")
    plt.plot(r_matlab,psi_matlab,label="Matlab")
    plt.title("compare psi in Python and matlab",fontsize=15)
    plt.xlabel("r",fontsize=15)
    plt.ylabel(r"$\psi$",fontsize=15)
    plt.legend()
    #plt.savefig(save_path + "")
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
        points_perturbed=4 #N-1
        ,ds=ds
        ,r=r_init,z=z_init,psi=psi_init
        ,delta_psi= 5e-1
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

    c = -ds*1.5

    plt.figure()
    plt.plot(radi_list[0],z_list[0],"o-",label="initial state")

    #psi_list[0] = [i + np.pi for i in psi_list[0]]
    
    print(f"init angle = {psi_list[0][0]}")
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
    plt.xlabel("r",fontsize=15)
    plt.ylabel("z",fontsize=15)
    plt.title(
        r"Line segment from the $\psi_{i}$ for a given (r,z) point"
        +f"\n to test if the angles are correct"
        ,fontsize=15
        )
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
    #Area_list, psi_list = const_args[12:14]
    #radi_list,z_list = const_args[14:16]

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,figs_for_video_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    df_name += f" dt={dt} and N={N}"

    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    #exit()
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    psi = df_sim['psi'][0]
    Area = df_sim['area list'][0]
    sim_steps = df_sim["sim_steps"][0]
    Tol= 1e-8

    f_cons, i_count_f = [[] for t in range(sim_steps)]  ,[[] for t in range(sim_steps)]
    g_cons ,i_count_g = [[] for t in range(sim_steps)]  ,[[] for t in range(sim_steps)]
    for t in range(sim_steps-1):
        for i in range(N):
            f = constraint_f(i=i,N=N,r=r[t],psi=psi[t],Area=Area)
            if f < -Tol or f > Tol :
                i_count_f[t].append(i)
                f_cons[t].append(f)

            g = constraint_g(i=i,N=N,r=r[t],z=z[t],psi=psi[t],Area=Area)
            if g < -Tol or g > Tol :
                i_count_g[t].append(i)
                g_cons[t].append(g)
    
    print(f"len(f[0])={len(f_cons[0])}")
    error_count = 0
    for t in range(sim_steps):
        if len(i_count_f[t]) >0:
            error_count += 1
            print(f"i_f={i_count_f[t]}")
    print(f"error count ={error_count} of {sim_steps} time steps")
    print(f"for a Tolerance of {Tol}")
    """plt.figure()
    plt.plot(i_count_f,f_cons,".-",label="f")
    plt.plot(i_count_g,g_cons,".-",label="g")
    plt.title(
        f"values of f,g that are non-zero \n"
        +f"N={N}, len(f)={len(f_cons)} and len(g)={len(g_cons)}"
    )
    plt.xlabel("index i")
    plt.ylabel("")
    plt.legend()
    plt.show()"""

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

def testing_integration_with_events():
    from Two_D_constants import Two_D_Constants_stationary_state
    from Two_D_functions import Epsilon_values,constraint_f,constraint_g, c_diff_f,c_diff_g, Epsilon_v2
    from two_d_continues_integration import dSds, descritize_sim_results,find_init_stationary_state
    print("\n \n \n")
    const_args = Two_D_Constants_stationary_state(
        print_val=False
        ,show_stationary_state=True
        ,pause_timer=0.1
    )
    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

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

    #args_list = (k ,sigma ,c0)
    args_list = (k ,sigma ,c0 ,tau ,kG)
    
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

    def dSds_test(s,S,k,sigma,c0,tau,kG):
        psi, r, z,n,lambs,nus,A = S
        #k,kG,sigma,c0 = p
        a,b =tau, kG
        drds = np.cos(psi)
        dzds = np.sin(psi)
        dlambs_ds = (k/2)*( (n-c0)**2 -np.sin(psi)**2/r**2) + sigma
        dnu_ds = 0

        dpsids  = n 
        dnds = np.sin(psi)*np.cos(psi)/r**2 - n*np.cos(psi)/r + lambs*np.sin(psi)/r #- nus*np.cos(psi)/(k*r)

        dAds = 2*np.pi*r

        return [dpsids ,drds ,dzds ,dnds ,dlambs_ds ,dnu_ds ,dAds]
    def edge_tension(t,y,k,sigma,c0,tau,kG):
        #tau = 1
        return tau - y[4]
    
    def edge_ratio(t,y,k,sigma,c0,tau,kG):
        dpsidt_1 = y[3]
        psi_1 = y[0]
        r_1 = y[1]
        alpha = kG/k
        val = (1-dpsidt_1)*r_1/np.sin(psi_1)-1 + alpha
        return val
    
    #edge_ratio.terminal = True
    edge_tension.terminal = True
    #The integration part.
    ans_odeint = scipy.integrate.solve_ivp(
        dSds_test
        ,t_span = [sN ,s0]
        ,t_eval = np.linspace(start=sN,stop=s0,num=10003)
        ,y0 = init_conditions
        ,args = args_list
        ,method="LSODA"#"RK45"
        ,atol=1e-10
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
    
    print(
    f"test ={(1-dpsidt_1)*r_1/np.sin(psi_1)-1 + alpha}"
    )

    print(ans_odeint.t_events)
    plt.figure()
    plt.plot(r[len(r)-1],z[len(r)-1],marker="o")
    plt.plot(r,z)
    plt.show()



def testing_for_no_correction_on_initial_state():
    from Two_D_constants import Two_D_Constants_stationary_state
    from Two_D_functions import Epsilon_v2,c_diff,Epsilon_values
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
    
    t=0
    epsilon = Epsilon_v2(
            N=N, r=radi_list[t], z=z_list[t] ,psi=psi_list[t] ,Area=Area_list
                    )
    
    efb = epsilon[0:N]
    egb = epsilon[N:2*N]
    i_list = [i for i in range(N)]
    plt.figure()
    plt.plot(i_list,efb,".-",label=r"$\epsilon_{f,i}$")
    plt.plot(i_list,egb,".-",label=r"$\epsilon_{g,i}$")
    plt.xlabel("i")
    plt.ylabel(r"$\epsilon_{i}$")
    plt.title(
        r"Values of the $\epsilon_\beta$ multipliers found for area correction "
        +f"\n on the initial state"
        )
    plt.legend()
    plt.show()


    print(f"epsilon vec = {epsilon}")


def testing_gradient_of_constraints():
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
    fig, ax = plt.subplots(nrows=3,ncols=2)
    plt.get_current_fig_manager().full_screen_toggle()
    for i in range(N):
        
        #manager = plt.get_current_fig_manager()
        #manager.resize(*manager.window.maxsize())

       

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
    

def testing_gradient_for_S():
    from Two_D_functions import Q_function,B_function, dSdpsi_func
    from two_d_data_processing import E_pot
    from Two_D_constants import Two_D_Constants_stationary_state

    const_args = Two_D_Constants_stationary_state(
        print_val=False
        ,show_stationary_state=False
        ,start_flat=False
        ,perturb=False
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]
    r_unperturbed, z_unperturbed = const_args[16:18]
    eta = const_args[18]

    grad_S_r,grad_S_psi = [],[]
    dSdr,dSdpsi = [],[]
    
    h = 1e-8
    error_tolerence = 1e-4

    def increase_by_h(dh,x,j):
        var_list = []
        for i in range(len(x)):
            if i == j:
                var_list.append(x[i] + dh)
            else:
                var_list.append(x[i])
        return var_list
    
    def find_error_index(
            j,error_tol,
            Delta_constraint
            ,diff_constraint
            ):
        null_val =""
        return_i ,return_j = null_val,null_val
        if Delta_constraint[j] == 0:
                if  -error_tol < diff_constraint[j] and diff_constraint[j] < error_tol:
                    return_i = null_val
                    return_j = null_val
                else:
                    return_i = i
                    return_j = j
        else:
            if abs(diff_constraint[j]/Delta_constraint[j]) < 1 - error_tolerence or abs(diff_constraint[j]/Delta_constraint[j]) > 1 + error_tolerence:
                return_i = i
                return_j = j
        
        return return_i,return_j, null_val
    
    def determine_error(i,dfdq,grad_f_q,tol):
        val =  abs( dfdq - grad_f_q )
        #print(f"val={val}")
        if tol < val:
            return i
        else:
            return ""

    
    # Finding the two gradient for r
    
    for j in range(N):
        # Finding the dS_i_dr_j parts
        dSdr.append(
            -Q_function(i=j,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau
                        ,Area=Area_list,psi=psi_list[0],radi=radi_list[0])
            )

        Si_rj_h = E_pot(N=N,k=k,kG=kG,sigma=sigma,tau=tau,c0=c0
            ,Area=Area_list,z=z_list[0],psi=psi_list[0]
            ,r = increase_by_h(dh=h,x=radi_list[0],j=j)
            )

        Si_rj = E_pot(N=N,k=k,kG=kG,sigma=sigma,tau=tau,c0=c0,z=z_list[0],psi=psi_list[0],Area=Area_list
            ,r=radi_list[0]
            )

        grad_S_r.append(
            (Si_rj_h - Si_rj)/h
            )

        # Fidning the dS_i_dpsi_j parts.
        dSdpsi.append(
            dSdpsi_func(i=j,N=N,c0=c0,k=k,kG=kG
                    ,r=radi_list[0] ,psi=psi_list[0],Area=Area_list)
                )

        Si_psij_h = E_pot(N=N,k=k,kG=kG,sigma=sigma,tau=tau,c0=c0
            ,z=z_list[0]
            ,psi=increase_by_h(dh=h,x=psi_list[0],j=j)
            ,Area=Area_list
            ,r=radi_list[0]
        )

        Si_psij = E_pot(N=N,k=k,kG=kG,sigma=sigma,tau=tau,c0=c0
                        ,z=z_list[0],psi=psi_list[0],r=radi_list[0]
                        ,Area=Area_list
                        )

        grad_S_psi.append(
            (Si_psij_h - Si_psij)/h
            )



    # Finding errors
    error_index_S_r,error_index_S_psi = [],[]
    for j in range(N):

        err_r = determine_error(i=j,dfdq=dSdr[j] ,grad_f_q=grad_S_r[j],tol=error_tolerence)
        if err_r != "":
            error_index_S_r.append(err_r)
        
        err_psi = determine_error(i=j,dfdq=dSdpsi[j],grad_f_q=grad_S_psi[j],tol=error_tolerence)
        if err_psi != "":
            error_index_S_psi.append(err_psi)

        #print(f"err_r:{err_r} and err_psi:{err_psi}")
    print(f"dSdr err ={error_index_S_r}")
    print(f"dSdpsi err ={error_index_S_psi}")
    


def find_overflow_error():
    from Two_D_functions import Q_function,B_function
    from two_d_data_processing import E_pot
    from Two_D_constants import Two_D_Constants_stationary_state

    const_args = Two_D_Constants_stationary_state(
        print_val=False
        ,show_stationary_state=False
        ,start_flat=False
        ,perturb=False
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]
    r_unperturbed, z_unperturbed = const_args[16:18]
    eta = const_args[18]


    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,figs_for_video_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    df_name += f" N,ds,dt,T={N,ds,dt,T}" # c0={c0} tau={tau}"#f" ds={dt}and N={N} and ds={ds} c0={c0} kG={kG}"

    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())

    print(radi_list[1])    

    exit()
    corr_count = df_sim["correction count"][0]
    time = [i*dt for i in range(sim_steps-1)]
    float_count = 0
    for t in range(sim_steps-1):
        for i in range(N):
            if isinstance(radi_list[t][i],float) != True:
                print(f"radi_list[t][i]={radi_list[t][i]}")
            if isinstance(radi_list[t][i],float):
                float_count += 1
    
        if radi_list[t].all() == 0:
            print(f"t={t*dt}")
            print(f"radi={radi_list[t]}")
            exit()




def test_convergence_of_alpha():
    from Two_D_constants import Two_D_Constants_stationary_state


    tau_list = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    psi2_L_list = [
        -0.3264e-4,-0.1634e-4,-0.0522e-4,-0.0322e-4,-0.0207e-4,-0.0139e-4,-0.0096e-4,-0.0050e-4
        ,-0.0037e-4,-0.0028e-4,-0.0021e-4,-0.0017e-4,-0.0013e-4,-0.0010e-4,-0.0008e-4,-0.0007e-4
        ,-0.0005e-4,-0.0004e-4,-0.0003e-4,-0.0003e-4,-0.0002e-4,-0.0002e-4,-0.0002e-4,-0.0001e-4
    ]

    sigma_list = [
        -0.4000,-0.3696,-0.3089,-0.2785,-0.2481,-0.2177,-0.1873,-0.1266,-0.0962
        ,-0.0658,-0.0354,-0.0051,0.0253,0.0557,0.0861,0.1165
        ,0.1468,0.1772,0.2076,0.2380,0.2684,0.2987,0.3291,0.3595
    ]

    alpha_discrete_list = [[] for i in range(len(tau_list))]
    alpha_list = []
    ds_min = 1e-10
    ds_max = 1.5e-2*1.5
    ds_vec = np.linspace(ds_max,ds_min, 100)
    N0 = 20
    ds0 = 1.5e-2
    L0 = N0*ds0
    alpha_ref,N_ref,ds_ref , i_ref= -1e3 ,-1e3 ,-1e3 , 0
    for j in range(len(tau_list)):
        print(f"{int(j*100/len(tau_list))}%",end="\r")
        for i in range(len(ds_vec)):
            ds = ds_vec[i]
            n = 3#int(L0/ds)
            #print(f"n={n}, ds={ds:e} , sigma={sigma_ref} , psi2={psi2_ref} ,tau={tau_ref}")

            const_args = Two_D_Constants_stationary_state(
                print_val=False
                ,show_stationary_state=False
                ,start_flat=False
                ,perturb=False
                ,N=n
                ,ds=ds
                ,tilde_sigma=sigma_list[j], psi_L=psi2_L_list[j] ,tilde_tau=tau_list[j]
            )
            L,r0,N,ds,T,dt = const_args[0:6]
            k,c0,sim_steps = const_args[6:9]
            #sigma, tau, kG = const_args[9:12]
            Area_list, psi_list = const_args[12:14]
            radi_list,z_list = const_args[14:16]
            r_unperturbed, z_unperturbed = const_args[16:18]
            eta = const_args[18]
            dpsi_perturb , psi_unperturb ,psi_L, alpha = const_args[19:23]

            dpsids = (psi_list[0][1] - psi_list[0][0])/ds
            alpha_discrete = (c0 - dpsids)*radi_list[0][0]/np.sin(psi_list[0][0]) - 1

            alpha_discrete_list[j].append(alpha_discrete/alpha)
            if i == 0:
                alpha_list.append(alpha)
        
    
    bad_sigma = 10 #sigma_list[len(sigma_list)-2]
    fig, ax = plt.subplots()
    for j in range(len(tau_list)):
        if sigma_list[j] != bad_sigma:
            plt.plot(ds_vec,alpha_discrete_list[j],".-",label=r"$\sigma$="+f"{sigma_list[j]}")
            plt.plot(ds_vec[0],alpha_discrete_list[j][0],"*")
            
    
    plt.hlines(y=1,xmin=-1,xmax=max(ds_vec)*1.3,linestyles="dashed")
    plt.xlabel("ds")
    #plt.legend()
    plt.xlim(0,max(ds_vec)*1.1)
    plt.show()


def Testing_RungeKutta():
    from Two_D_constants import Two_D_Constants_stationary_state
    from Runge_Kutta import RungeKutta45
    const_args = Two_D_Constants_stationary_state(
        print_val=True
        ,show_stationary_state=True
        ,start_flat=False
        ,perturb=False
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]
    r_unperturbed, z_unperturbed = const_args[16:18]
    eta,dpsi_perturb_val,psi_unperturbed = const_args[18:21]
    psi2_init ,alpha = const_args[21:23]

    t=0
    lambs,nus = Langrange_multi(
                N=N,k=k,c0=c0,sigma=sigma
                ,kG=kG,tau=tau,ds=ds,eta=eta
                ,Area=Area_list
                ,psi=psi_list[t]
                ,radi=radi_list[t]
                ,z_list=z_list[t]
            )
    kr,kz,kpsi = RungeKutta45(
        N=N,dt=dt,k=k,c0=c0, sigma=sigma
        ,kG=kG ,tau=tau, ds=ds,eta=eta
        ,Area=Area_list 
        ,psi_init=psi_list[t],r_init=radi_list[t], z_init=z_list[t]
        ,lamb=lambs , nu=nus
    )
    #print(kr)
    print(np.shape(kr))
    print("kr=",kr)
    print(f"kr[:][0]={kr[:,0]}")
    print(f"sum(kr[:][0])={np.sum(kr[:,0]):.1e}")



def test_flat_model_object():
    from Surface_class import Surface_membrane
    from two_d_data_processing import get_files
    N = 30
    T= 4e-8 #1e-7 #5e-9'
    dt = 1e-11
    tau = 0#1.0e3 #max = 1.0
    var_corr = [
        1e-5, 
        1e-4,
        2.5e-4,
        5e-4,
        6.25e-4,
        7.5e-4,
        9.75e-4,
        1e-3,
        1e-2
        ]
    for i in range(len(var_corr)):
        save_path = f"2D sim results\\rolling test\\tau={tau:0.1e} dt={dt}\\"
        membrane = Surface_membrane(
            N=N, T=T , dt = dt
            ,const_index = 0
            ,save_path= save_path +f"varcorr={var_corr[i]:0.1e}\\"
        )

        membrane.start_flat = True
        membrane.tau = tau
        membrane.const_length_diff_N_density = False
        membrane.var_corr_tol = var_corr[i]

        #membrane.make_plots = True
        #membrane.make_movie = True
        #membrane.close_final_plots = True
        #membrane.init_config_show_time = 40
        #membrane.setup_simulation()
        #membrane.run_sim()
        #membrane.plotting_n_movie_data()
    
    files = get_files(save_path)

    for i in range(len(files)-1):
        for _ in range(30):
            df1 = pd.read_pickle(files[i])
            df2 = pd.read_pickle(files[i+1])

            tol1 = df1["Tolerance"][0]
            tol2 = df2["Tolerance"][0]
            if tol1 < tol2 :
                files[i] ,files[i+1] = files[i+1], files[i]


    fig,ax = plt.subplots()
    d = 0
    for data in files:
        df = pd.read_pickle(data)
        #print(df.info())
        #print(data)
        Epot = df["Epot"][0]
        Epot_before = df["Epot before correction"][0]
        corr_count = df["correction count"][0]
        sim_steps = df["sim_steps"][0]
        dt = df["dt"][0]
        tol_varr = df["Tolerance"][0]
        reduce = 0
        time_vec = np.linspace(0,(sim_steps-reduce)*dt , sim_steps-reduce )

        diff_Epot = []
        count = 0
        for i in range(sim_steps-1):
            dE = Epot[i+1] - Epot[i]
            if dE > 0:
                count += 1
                #diff_Epot.append(1)
            elif dE <= 0:
                pass #diff_Epot.append(0)
            diff_Epot.append(count)
        #print(len(diff_Epot))

        #if any(diff_Epot) != 0:
        ax.plot(
            time_vec[0:sim_steps-1],diff_Epot
            ,label=f"tol={tol_varr:.1e} , num violations={np.sum(diff_Epot)}"
            ,marker="."
            )
    
    plt.xlabel("time [s]")
    plt.ylabel("num times Epot[t+1]-Epot[t] > 0")
    plt.title("Checking the amount of times that the \n potential energy increases during the simulation.")
    plt.legend()
    plt.grid()
    plt.draw()
    plt.pause(0.5)
    plt.savefig(save_path + "compare when dEmorethan0.png")


    fig,ax = plt.subplots()
    for data in files:
        df = pd.read_pickle(data)
        Epot = df["Epot"][0]
        sim_steps = df["sim_steps"][0]
        dt = df["dt"][0]
        tol_varr = df["Tolerance"][0]
        reduce = 0#1#2
        time_vec = np.linspace(0,(sim_steps-reduce)*dt , sim_steps-reduce )

        ax.plot(
            time_vec,Epot
            ,label=f"tol={tol_varr:.1e}"
        )


    plt.title("all Epot")
    plt.xlabel("time [s]")
    plt.ylabel("Epot [zJ]")
    plt.legend()

    df = pd.read_pickle(files[0])
    print(df.info())
    fig,ax = plt.subplots()
    for data in files:
        df = pd.read_pickle(data)
        Epot = df["Epot"][0]
        count_corr = df["correction count"][0]
        sim_steps = df["sim_steps"][0]
        tol_varr = df["Tolerance"][0]
        sum_corrs = []
        k = 0
        diff_Epot = []
        count = 0
        for i in range(sim_steps-2):
            k += count_corr[i]
            sum_corrs.append(k)

            dE = Epot[i+1] - Epot[i]
            if dE > 0:
                count += 1
            diff_Epot.append(count)
        
        ax.plot(
            sum_corrs
            ,diff_Epot
            ,label=f"tol={tol_varr:.1e}"
            ,marker="."
            )

    plt.xlabel("sum corrections")
    plt.ylabel("")
    plt.legend()


    fig,ax = plt.subplots()
    d = 0
    for data in files:
        df = pd.read_pickle(data)
        #print(df.info())
        #print(data)
        S_pot = df["Epot"][0]
        S_pot_before = df["Epot before correction"][0]
        corr_count = df["correction count"][0]
        sim_steps = df["sim_steps"][0]
        dt = df["dt"][0]
        tol_varr = df["Tolerance"][0]
        reduce = 0
        time_vec = np.linspace(0,(sim_steps-reduce)*dt , sim_steps-reduce )

        corr_place_y = []
        corr_place_x = []

        Epot_x, Epot_y = [], []
        difference_Epot = np.zeros(sim_steps)
        for t in range(sim_steps):
            if corr_count[t] != 0:
                corr_place_y.append( (S_pot[t] + S_pot_before[t])/2)
                corr_place_x.append( time_vec[t] )
            
            Q = S_pot[t] - S_pot_before[t]
            if Q > 0.005:
                Epot_x.append(time_vec[t])
                Epot_y.append(Q)
        
        if len(Epot_y) > 0 :
            ax.plot(
                Epot_x,Epot_y
                ,label=f"tol={tol_varr:0.1e}"
                ,marker="."
                #,markersize=10
            )
        
    plt.title(r"$ E_{pot} - E_{before,pot}$")
    plt.xlabel("dt")
    plt.ylabel(r"$ \Delta E_{pot} $")
    plt.grid()

    plt.legend()


    plt.show()
    


def test_gradients_again(
    data_path:str = "2D sim results\\obj\\plus larger ds\\N=40\\(N,T,dt)=(40,1.0e-07,1.0e-11)\\",
    compare_df_name:str = "compare_df.pkl",
    h:float = 1e-9,
    hpsi:float  = 1e-9,
    make_new_data:bool = True
):
    from Two_D_functions import Q_function, dzdt_func,dpsidt_func,gamma, c_diff_f,c_diff_g, dSdpsi_func
    from Two_D_functions import drdt_func, constraint_f,constraint_g, Langrange_multi, Q_function
    from two_d_data_processing import E_pot, get_files
    

    def constraints_multi(
            N,Area
            ,nus,lambs
            ,r,psi,z
        ):
        grad_constraint = 0
        for i in range(N):
            grad_constraint += (
                lambs[i]*constraint_f(i=i,N=N,r=r,psi=psi,Area=Area)
                + nus[i]*constraint_g(i=i,N=N,r=r,psi=psi,z=z,Area=Area)
            )
        return grad_constraint

    def theory_constraint(
            j,N,r,z,psi,Area,nus,lambs,diff_var
        ):

        dc = 0
        for i in range(N):
            dc += (
                lambs[i]*c_diff_f(i=i,j=j,N=N,r=r,psi=psi,Area=Area,diff_var=diff_var)
                + nus[i]*c_diff_g(i=i,j=j,N=N,r=r,z=z,psi=psi,Area=Area,diff_var=diff_var)
            )
        
        return dc

     
    file = get_files(data_path)
    df = pd.read_pickle(file[0])
    #print(df.info())
    sim_steps = df["sim_steps"][0]

    dt = df["dt"][0]
    N = df["N"][0]
    k = df["k"][0]
    kG = df["kG"][0]
    c0 = df["c0"][0]
    sigma = df["sigma"][0]
    tau = df["tau"][0]
    eta =  1 #df["eta"][0] 
    ds = df["ds"][0]
    Area = df["area list"][0]
    S_pot = df["Epot"][0]
    S_pot_before = df["Epot before correction"][0]
    corr_count = df["correction count"][0]

    r = df["r"][0]
    z = df["z"][0]
    psi = df["psi"][0]


    time_vec = np.linspace(0,(sim_steps)*dt,sim_steps)

    grad_test_dLdr = np.zeros(shape=(sim_steps,N+1))
    grad_test_dLdz = np.zeros(shape=(sim_steps,N+1))
    grad_test_dLdpsi = np.zeros(shape=(sim_steps,N))

    grad_test_r = np.zeros(shape=(sim_steps,N+1))
    grad_test_z = np.zeros(shape=(sim_steps,N+1))
    grad_test_psi = np.zeros(shape=(sim_steps,N))

    grad_test_constraint_r = np.zeros(shape=(sim_steps,N+1))
    grad_test_constraint_z = np.zeros(shape=(sim_steps,N+1))
    grad_test_constraint_psi = np.zeros(shape=(sim_steps,N))

    further_data = []
    
    sim_steps = int(sim_steps/100)
    tol = 1e-2 # for the tolerance of deviation
    print("\n \n")
    if make_new_data == True:
        print_scale = sim_steps/1000
        start_time = time.time()
        for t in range(sim_steps):
            if int(t%print_scale) == 0 :
                time1 = time.time()-start_time
                time_left = (time1/(t+1))*(sim_steps-t)
                time_h_start, time_m_start, time_s_start = int((time1/60**2)%24), int((time1/60)%60), int(time1%60)
                time_h_end, time_m_end, time_s_end = int((time_left/60**2)%24), int((time_left/60)%60), int(time_left%60)
                print(
                    f"completion : {round(t/(print_scale*10),1)}%       " 
                    +f"Time since start = {time_h_start}h {time_m_start}m {time_s_start}s        "
                    +f"Estimated time left = {time_h_end}h {time_m_end}m {time_s_end}s"
                    , end="\r"
                )

            lambs, nus = Langrange_multi(
                N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,ds=ds,eta=eta,Area=Area
                ,psi=psi[t],radi=r[t],z_list=z[t]
                )
            
            S = - E_pot(N=N,k=k,kG=kG,tau=tau,c0=c0,r=r[t],psi=psi[t],Area=Area)
            constraint_eq = constraints_multi(N=N,r=r[t],psi=psi[t],z=z[t],Area=Area
                                              ,nus=nus ,lambs=lambs)

            grad_L_ref = S + constraint_eq            

            for i in range(N):
                rh = [ r[t][n] + h if n==i else r[t][n] for n in range(N+1)]
                zh = [ z[t][n]+ h  if n==i else z[t][n] for n in range(N+1)]
                psih = [ psi[t][n]+ hpsi  if n==i else psi[t][n] for n in range(N)]

                """------------- Gradient test for r ---------------------------------------------------------------------"""
                lambs_rh ,nus_rh= Langrange_multi(
                    N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,ds=ds,eta=eta,Area=Area
                    ,psi=psi[t] ,z_list=z[t]
                    ,radi=rh
                )                
                                
                S_rh = - E_pot(N=N,k=k,kG=kG,tau=tau,c0=c0,psi=psi[t],Area=Area
                               ,r=rh
                            )

                constraint_eq_rh = constraints_multi(
                    N=N,Area=Area ,psi=psi[t], z=z[t]
                    ,nus=nus_rh ,lambs=lambs_rh
                    ,r=rh
                    )

                """grad_L_rh = S_rh + constraint_eq_rh 

                dLdrh = (grad_L_rh - grad_L_ref)/h  #the Newton derivative 

                dLdr = gamma(i=i,ds=ds,eta=eta)*drdt_func(
                            i=i,N=N,k=k,c0=c0,sigma=sigma,kG=kG
                            ,tau=tau,ds=ds,eta=eta
                            ,Area=Area
                            ,lamb=lambs,nu=nus
                            ,psi=psi[t] ,radi=r[t] ,z_list=z[t]
                            )"""
                

                diff_constraints_rh = (constraint_eq_rh - constraint_eq)/h                
                diff_constraints_r = theory_constraint(j=i,N=N,r=r[t],z=z[t],psi=psi[t],Area=Area,nus=nus,lambs=lambs,diff_var="r")

                Q_r = Q_function(i=i,N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,Area=Area,psi=psi[t],radi=r[t])

                grad_S_r = (S_rh - S)/h

                #grad_test_dLdr[t][i] = (dLdr - dLdrh)/dLdr #gets a % deviation of the theoretical result
                grad_test_constraint_r[t][i] = (diff_constraints_r - diff_constraints_rh)/diff_constraints_r
                grad_test_r[t][i] = (Q_r - grad_S_r)/Q_r

                if abs(grad_test_constraint_r[t][i]) > tol:
                    further_data.append([
                        grad_test_constraint_r[t][i],
                        f"grad_test_constraint_r[t][i] = (diff_constraints_r - diff_constraints_rh)/diff_constraints_r = ({diff_constraints_r} - {diff_constraints_rh})/{diff_constraints_r} = ",
                        f"grad test constraint r (t*dt,i)=({t*dt:.3e},{i}):  = {grad_test_constraint_r[t][i]} ",
                        f"lambs rh={lambs_rh[i]}" ,
                        f"lambs   ={lambs[i]} ",
                        f"nus rh  ={nus_rh[i]}" ,
                        f"nus     ={nus[i]}",
                        f"constraint eq rh={constraint_eq_rh}   ",
                        f"constraint eq   ={constraint_eq}   ",
                        f"S rh={S_rh}" ,
                        f"S   ={S}",
                        f"grad S r={grad_S_r}  " ,
                        f"Q r     ={Q_r}   ",
                        f"diff_constraints_rh={diff_constraints_rh} " ,
                        f"diff_constraints_r ={diff_constraints_r}    ",
                    ])
                                        

                if abs(grad_test_r[t][i]) > tol:
                    further_data.append([
                        grad_test_r[t][i] ,
                        f"grad_test_r[t][i] = (Q_r - grad_S_r)/Q_r =({Q_r} - {grad_S_r})/{Q_r}",
                        f"grad test r: (t,i)={t,i}:  = {grad_test_r[t][i]} ",
                        f"lambs rh={lambs_rh[i]} ",
                        f"lambs  ={lambs[i]} ",
                        f"nus rh={nus_rh[i]}",
                        f"nus   ={nus[i]}",
                        f"constraint eq rh={constraint_eq_rh}",
                        f"constraint eq   ={constraint_eq}",
                        f"S rh ={S_rh} ",
                        f"S    ={S}",
                        f"grad S r ={grad_S_r}",
                        f"Q r      ={Q_r}",
                        f"diff_constraints_rh ={diff_constraints_rh}",
                        f"diff_constraints_r  ={diff_constraints_r}"
                    ])
                    

                """------------- Gradient test for z --------------------------------------------------------------------------"""
                lambs_zh, nus_zh = Langrange_multi(
                    N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,ds=ds,eta=eta,Area=Area
                    ,psi=psi[t] ,radi=r[t]
                    ,z_list=zh
                )                
                                
                #S_zh = - E_pot(N=N,k=k,kG=kG,tau=tau,c0=c0,r=r[t],psi=psi[t],Area=Area)

                constraint_eq_zh = constraints_multi(
                    N=N,Area=Area ,r=r[t] ,psi=psi[t]
                    ,nus=nus_zh ,lambs=lambs_zh
                    ,z=zh
                    )

                """grad_L_zh = S_zh + constraint_eq_zh 

                dLdzh = (grad_L_zh - grad_L_ref)/h  #the Newton derivative 

                dLdz = gamma(i=i,ds=ds,eta=eta)*dzdt_func(
                    i=i,ds=ds,eta=eta,Area=Area,radi=r[t],nu=nus
                )"""
                

                diff_constraints_zh = (constraint_eq_zh - constraint_eq)/h                
                diff_constraints_z = theory_constraint(j=i,N=N,r=r[t],z=z[t],psi=psi[t],Area=Area,nus=nus,lambs=lambs,diff_var="z")

                
                #grad_S_z = (S_zh - S)/h

                #grad_test_dLdz[t][i] = (dLdz - dLdzh)/dLdz  #gets a % deviation of the theoretical result
                grad_test_constraint_z[t][i] = (diff_constraints_z - diff_constraints_zh)/diff_constraints_z
                #grad_test_z[t][i] = (diff_constraints_z - diff_constraints_zh)/diff_constraints_z 

                if abs(grad_test_constraint_z[t][i]) > tol:
                    further_data.append([
                        grad_test_constraint_z[t][i],
                        f"grad_test_constraint_z[t][i] = (diff_constraints_z - diff_constraints_zh)/diff_constraints_z = ({diff_constraints_z} - {diff_constraints_zh})/{diff_constraints_z}",
                        f"grad test constraint z: (t*dt,i)=({t*dt:.3e},{i}):  = {grad_test_constraint_z[t][i]} ",
                        f"lambs zh={lambs_zh[i]}",
                        f"lambs   ={lambs[i]}",
                        f"nus zh  ={nus_zh[i]}",
                        f"nus     ={nus[i]}",
                        f"constraint eq zh={constraint_eq_zh}",
                        f"constraint eq   ={constraint_eq}",
                        f"diff_constraints_zh={diff_constraints_zh}",
                        f"diff_constraints_z ={diff_constraints_z}"
                    ])


                """------------- Gradient test for psi --------------------------------------------------------------"""

                lambs_psih, nus_psih = Langrange_multi(
                    N=N,k=k,c0=c0,sigma=sigma,kG=kG,tau=tau,ds=ds,eta=eta,Area=Area
                    ,z_list=z[t],radi=r[t]
                    ,psi=psih 
                )

                #dpsids = d
                dSdpsi = -dSdpsi_func(i=i,N=N,c0=c0,k=k,kG=kG,r=r[t],psi=psi[t],Area=Area)

                S_psih = - E_pot(N=N,k=k,kG=kG,tau=tau,c0=c0,r=r[t],Area=Area
                                 ,psi=psih
                                 )

                diff_constraints_psi = theory_constraint(j=i,N=N,r=r[t],z=z[t],psi=psi[t],Area=Area,nus=nus,lambs=lambs,diff_var="psi")
                
                constraint_eq_psih = constraints_multi(
                    N=N,Area=Area ,r=r[t], z=z[t]
                    ,nus=nus_psih ,lambs=lambs_psih
                    ,psi=psih
                    )
                
                grad_S_psi = (S_psih - S)/hpsi
                
                """grad_L_psih = S_psih + constraint_eq_psih

                dLdpsih = (grad_L_psih - grad_L_ref)/h
                
                dLdpsi = (dSdpsi + diff_constraints_psi)"""
                
                diff_constraints_psih = (constraint_eq_psih - constraint_eq)/hpsi

                grad_test_psi[t,i] = (dSdpsi - grad_S_psi)/dSdpsi
                #grad_test_dLdpsi[t,i] = ( dLdpsi - dLdpsih )/dLdpsi 
                grad_test_constraint_psi[t,i] = (diff_constraints_psi - diff_constraints_psih)/diff_constraints_psi
                
                if abs(grad_test_constraint_psi[t][i]) > tol:
                    further_data.append([
                        grad_test_constraint_psi[t,i],
                        f"grad_test_constraint_psi[t,i] = (diff_constraints_psi - diff_constraints_psih)/diff_constraints_psi = ({diff_constraints_psi} - {diff_constraints_psih})/{diff_constraints_psi}",
                        f"grad test constraint psi:(t*dt,i)=({t*dt:.3e},{i}):  = {grad_test_constraint_psi[t][i]} ",
                        f"lambs psih={lambs_psih[i]}",
                        f"lambs     ={lambs[i]}",
                        f"nus psih  ={nus_psih[i]}",
                        f"nus       ={nus[i]}",
                        f"constraint eq psih={constraint_eq_psih}",
                        f"constraint eq     ={constraint_eq:}",
                        f"S psih={S_psih}",
                        f"S     ={S}",
                        f"grad S psi={grad_S_psi}",
                        f"dSdpsi={dSdpsi}",
                        f"diff_constraints_psih={diff_constraints_psih}",
                        f"diff_constraints_psi ={diff_constraints_psi}    "
                    ])

                if abs(grad_test_psi[t][i]) > tol:
                    further_data.append([
                        grad_test_psi[t,i],
                        f"grad_test_psi[t,i] = (dSdpsi - grad_S_psi)/dSdpsi = ({dSdpsi} - {grad_S_psi})/{dSdpsi}",
                        f"grad test psi: (t*dt,i)=({t*dt:.3e},{i}):  = {grad_test_psi[t][i]:.3e} ",
                        f"lambs psih={lambs_psih[i]}",
                        f"lambs     ={lambs[i]}",
                        f"nus psih  ={nus_psih[i]}",
                        f"nus       ={nus[i]}",
                        f"constraint eq psih={constraint_eq_psih}",
                        f"constraint eq     ={constraint_eq:}",
                        f"S psih={S_psih}",
                        f"S     ={S}",
                        f"grad S psi={grad_S_psi}",
                        f"dSdpsi={dSdpsi}",
                        f"diff_constraints_psih={diff_constraints_psih}",
                        f"diff_constraints_psi ={diff_constraints_psi}    "
                    ])

                """------------- Saving the data --------------------------------------------------------------------"""

                df2 = pd.DataFrame({
                    'grad test r': [grad_test_r],
                    'grad test z': [grad_test_z],
                    'grad test psi': [grad_test_psi],
                    'grad test dLdr': [grad_test_dLdr],
                    'grad test dLdz': [grad_test_dLdz],
                    'grad test dLdpsi': [grad_test_dLdpsi],
                    'grad test constraint r':[grad_test_constraint_r],
                    'grad test constraint z':[grad_test_constraint_z],
                    'grad test constraint psi':[grad_test_constraint_psi],
                    'further data':[further_data],
                    'h':h,
                    'hpsi':hpsi
                })

                df2.to_pickle(data_path + compare_df_name)

    elif make_new_data == False:
        df_compare_grad = pd.read_pickle(data_path + compare_df_name)
        grad_test_r = df_compare_grad["grad test r"][0]
        #grad_test_z = df_compare_grad["grad test z"][0]
        grad_test_psi = df_compare_grad["grad test psi"][0]

        grad_test_dLdr = df_compare_grad["grad test dLdr"][0]
        grad_test_dLdz = df_compare_grad["grad test dLdz"][0]
        grad_test_dLdpsi = df_compare_grad["grad test dLdpsi"][0]

        grad_test_constraint_r =  df_compare_grad["grad test constraint r"][0]
        grad_test_constraint_z =  df_compare_grad["grad test constraint z"][0]
        grad_test_constraint_psi =  df_compare_grad["grad test constraint psi"][0]
        further_data = df_compare_grad["further data"][0]
        h = df_compare_grad["h"][0]
        hpsi = df_compare_grad["hpsi"][0]
        
        
    
    print("\n \n")

    for d in further_data[0]:
        print(d)

    print("\n \n")
    d_max,imax = 0,0
    for i in range(len(further_data)):
        d = abs(further_data[i][0])
        if d > d_max:
            d_max = d
            imax = i

    print("The maximum deviation is: \n")
    for d in further_data[imax]:
        print(d)

    data = [
        grad_test_r,
        #grad_test_z,
        grad_test_psi,
        #grad_test_dLdr,
        #grad_test_dLdz,
       # grad_test_dLdpsi,
        grad_test_constraint_r,
        grad_test_constraint_z,
        grad_test_constraint_psi
    ]

    data_names = [
        r"Grad test r : $ \dfrac{  \dfrac{ \partial S_{theory} }{ \partial r_{i} }  - ( ~ S(r_i + h) - S(r_i) ~ ) }{  \dfrac{ \partial S_{theory} }{ \partial r_{i} } } $ ",
        #r"Grad test z : $ \dfrac{  \dfrac{ \partial S_{theory} }{ \partial z_{i} }  - ( ~ S(z_i + h) - S(z_i) ~ ) }{  \dfrac{ \partial S_{theory} }{ \partial z_{i} } } $ ",
        r"Grad test $\psi$ : $ \dfrac{  \dfrac{ \partial S_{theory} }{ \partial \psi_{i} }  - ( ~ S(\psi_i + h) - S(\psi_i) ~ ) }{  \dfrac{ \partial S_{theory} }{ \partial \psi_{i} } } $ ",
        #f"Grad test dLdr : " + r"$ \dfrac{dLdr_{theory} - dLdr_{r_i+h} }{dLdr_{theory}} $ ",
        #"Grad test dLdz : " + r"$ \dfrac{dLdz_{theory} - dLdz_{z_i+h} }{dLdz_{theory}} $ ",
        #"Grad test dLdz : " + r"$ \dfrac{dLd \psi_{theory} - dLd \psi _{\psi_i+h} }{dLd \psi _{theory}} $ ",
        "Grad test constraint r : " + r"$ \dfrac{constriant_{theory} - constriant_{r_i+h} }{constriant_{theory}} $ ",
        "Grad test constraint z : " + r"$ \dfrac{constriant_{theory} - constriant_{z_i+h} }{constriant_{theory}} $ ",
        "Grad test constraint psi : " + r"$ \dfrac{constriant_{theory} - constriant_{\psi_i+h} }{constriant_{theory}} $ ",
    ]

    save_name_list = [
        "dSdr",
        "dSdpsi",
        "grad constriant r",
        "grad constriant z",
        "grad constriant psi"
    ]
    
    h_list = [
        h,hpsi,h,h,hpsi
    ]
    
    t_switch = [i for i in S_pot].index(min(S_pot))

    fig, ax = plt.subplots()
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    ax.plot(time_vec,S_pot,label="Epot")
    ax.vlines(x=time_vec[t_switch],ymin=min(S_pot),ymax=max(S_pot),label="dSdt > 0",color="k")

    ax.set_xlabel("time [s]")
    ax.set_ylabel("Epot")
    plt.legend()
    plt.grid()
    plt.draw()
    plt.pause(0.5)
    plt.savefig(data_path + "potential energy" + ".png")
    plt.savefig(data_path + "potential energy" + ".svg")



    for l in range(len(data)):
        i_problem_list = []
        for t in range(sim_steps):
            for i in range(N):
                if abs(data[l][t][i]) > tol and i not in i_problem_list:
                    # 1% deviation from the theory, is the value when 1e-2 = 1/100
                    i_problem_list.append(i)
                    
        i_problem_list.sort()
        fig, ax = plt.subplots(2,2)
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        N_start = 0
        dN_plot = int(N/4)
        for n in range(2):
            for m in range(2):
                for i in range(N_start,N_start + dN_plot):
                    if i == N_start:
                        ax[n,m].vlines(x=time_vec[t_switch],ymin=min(data[l][:,i]),ymax=max(data[l][:,i]),color="k",label="dSdt>0")
                    else:
                        ax[n,m].vlines(x=time_vec[t_switch],ymin=min(data[l][:,i]),ymax=max(data[l][:,i]),color="k")

                    ax[n,m].plot(
                        time_vec,data[l][:,i]
                        ,label=f"i={i}"
                    )

                    ax[n,m].grid("on")
                    ax[n,m].legend()
                    ax[n,m].ticklabel_format(style='sci', scilimits=(0,0))
                    
                    if n == 0:
                        ax[n,m].set_ylabel("")
                    if n == 1:
                        ax[n,m].set_xlabel("time [s]")

                N_start += dN_plot

        plt.suptitle(
           data_names[l] + f"   deviations greater than 1e-2:  i" + r"$\in$" + f"{i_problem_list}"
           +f"\n for h = {h_list[l]}"
        )

        plt.draw()
        plt.pause(0.5)
        plt.savefig(data_path + save_name_list[l] + ".png")
        plt.savefig(data_path + save_name_list[l] + ".svg")

    plt.show()
   

def compare_potentential_energy():
    from two_d_data_processing import get_files

    data_path = "2D sim results\\obj\\plus\\N=40\\"
    compare_df_name = "compare_df.pkl"
    file = get_files(data_path)
    df = pd.read_pickle(file[0])
    print(df.info())
    sim_steps = df["sim_steps"][0]

    dt = df["dt"][0]
    N = df["N"][0]
    k = df["k"][0]
    kG = df["kG"][0]
    c0 = df["c0"][0]
    sigma = df["sigma"][0]
    tau = df["tau"][0]
    eta =  1 #df["eta"][0] 
    ds = df["ds"][0]
    Area = df["area list"][0]
    S_pot = df["Epot"][0]
    S_pot_before = df["Epot before correction"][0]
    T_kin = df["Ekin"][0]
    corr_count = df["correction count"][0]

    time_vec = np.linspace(0,(sim_steps)*dt,sim_steps)

    corr_place_x,corr_place_y = [], []
    Epot_x, Epot_y = [], []
    Epot_path_x,Epot_path_y = [],[]
    for t in range(sim_steps):
        if corr_count[t] != 0:
            corr_place_y.append( (S_pot[t] + S_pot_before[t])/2)
            corr_place_x.append( time_vec[t] )
        
        Q = S_pot[t] - S_pot_before[t]
        if Q !=0:
            Epot_x.append(time_vec[t])
            Epot_y.append(Q)
        #difference_Epot[t] = S_pot[t] - S_pot_before[t]
        
        Epot_path_x.append(time_vec[t])
        Epot_path_y.append(S_pot_before[t])

        Epot_path_x.append(time_vec[t])
        Epot_path_y.append(S_pot[t])

    fig, ax = plt.subplots()
    ax.plot(
        time_vec,S_pot_before
        ,label="Epot before"
        ,marker="."
        ,linewidth=4
    )
    ax.plot(
        time_vec,S_pot
        ,label="Epot"
        ,marker="."
        ,linestyle="dashed"
        ,linewidth=4
    )
    """ax.plot(
        corr_place_x,corr_place_y
        ,label="time of corr"
        ,marker="|"
        ,markersize = 20
        #,linestyle=False
        #,linewidth=3
    )
    """
    ax.plot(
        Epot_path_x,Epot_path_y
        ,label="path taken"
        ,linestyle="-"#"dashed"
        ,marker="."
        ,color="k"
        )
    
    plt.title("potential energy")
    plt.legend()
    plt.grid()


    fig,ax = plt.subplots()
    ax.plot(
        #time_vec,difference_Epot
        Epot_x,Epot_y
        ,label="difference Epot"
        )
    
    
    plt.legend()
    plt.title(r"$ E_{pot} - E_{before,pot} $")
    plt.grid()


    fig,ax = plt.subplots()
    ax.plot(
        time_vec[1::],T_kin[1::]
        ,marker="."
    )
    #plt.legend()
    plt.title(r"$ E_{kin} $")
    plt.grid()

    plt.show()




def test_if_variable_correction_causes_Epot_increase():
    from two_d_data_processing import get_files
    from Two_D_functions import dpsidt_func,drdt_func,dzdt_func,Make_variable_corrections
    from Runge_Kutta import RungeKutta45
    data_path = "2D sim results\\object results\\RK4\\T=3e-07\\plus\\(N,T,dt)=(20,3.0e-07,2.5e-11)\\"
    #data_path = "2D sim results\\object results T=1e-06\\plus\\(N,T,dt)=(20,1.0e-06,2.5e-11)\\"
    #compare_df_name = "compare_df.pkl"
    file = get_files(data_path)
    

    df = pd.read_pickle(file[0])
    print(df.info())
    r_list = df["r"][0]
    z_list = df["z"][0]
    psi_list = df["psi"][0]

    r_ref = r_list.copy()
    z_ref = z_list.copy()
    psi_ref = psi_list.copy()

    r_var_corr = r_list.copy()
    z_var_corr = z_list.copy()
    psi_var_corr = psi_list.copy()

    Area_list = df["area list"][0]

    sim_steps = df["sim_steps"][0]
    dt = df["dt"][0]
    N = df["N"][0]
    k = df["k"][0]
    kG = df["kG"][0]
    c0 = df["c0"][0]
    sigma = df["sigma"][0]
    tau = df["tau"][0]
    eta =  1 #df["eta"][0] 
    ds = df["ds"][0]
    Area = df["area list"][0]
    S_pot = df["Epot"][0]
    S_pot_ref = S_pot.copy()
    #S_pot_before = df["Epot before correction"][0]
    T_kin = df["Ekin"][0]
    corr_count = df["correction count"][0]
    var_corr_tol = df["Tolerance"][0]
    

    overflow_err = False
    integration_method = "RK4"
    time_vec = np.linspace(0,sim_steps*dt,sim_steps)

    dEpotdt_before ,dEpotdt = np.zeros(sim_steps-1),np.zeros(sim_steps-1)
    compare_dEdt = np.zeros(sim_steps-1)
    t_points = np.zeros(9)
    m = 0
    for t in range(sim_steps-1):
        #dEpotdt_before[t] = (S_pot_before[t+1] - S_pot_before[t])/dt
        #dEpotdt[t] = (S_pot[t+1] - S_pot[t])/dt
        #compare_dEdt[t] = dEpotdt[t] - dEpotdt_before[t]

        if t%(int(sim_steps/9)) == 0:
            t_points[m] = int(t)
            m += 1
            if m == 9:
                break

    
    num_of_free_steps = 40
    do_correction = False#True
    for t_start in t_points:
        t_stop = int(t_start + num_of_free_steps)

        for t in range(int(t_start),t_stop):      
            lambs,nus = Langrange_multi(
                    N=N,k=k,c0=c0,sigma=sigma
                    ,kG=kG,tau=tau,ds=ds,eta=eta
                    ,Area=Area_list
                    ,psi=psi_list[t]
                    ,radi=r_list[t]
                    ,z_list=z_list[t]
                )            
                
            if integration_method == "RK4":
                kr,kz,kpsi = RungeKutta45(
                    N=N,dt=dt,k=k,c0=c0, sigma=sigma
                    ,kG=kG ,tau=tau, ds=ds,eta=eta
                    ,Area=Area_list,lamb=lambs,nu=nus
                    ,psi_init=psi_list[t],r_init=r_list[t], z_init=z_list[t]
                    )
                
                for i in range(N+1):
                    if i == N:
                        z_list[t+1][i] = z_list[t][i]
                        r_list[t+1][i] = r_list[t][i]
                        #psi[t+1][i] = psi[t][i]

                    if i < N:
                        r_list[t+1][i] = r_list[t][i] + (dt/6)*(kr[1][i] + 2*kr[2][i] + 2*kr[3][i] + kr[4][i])
                        z_list[t+1][i] = z_list[t][i] + (dt/6)*(kz[1][i] + 2*kz[2][i] +2* kz[3][i] + kz[4][i])
                        psi_list[t+1][i] = psi_list[t][i] + (dt/6)*(kpsi[1][i] + 2*kpsi[2][i] + 2*kpsi[3][i] + kpsi[4][i])

                    
                    S_pot[t] = E_pot(
                    N=N,k=k,kG=kG,tau=tau,c0=c0
                    ,r=r_list[t],psi=psi_list[t],Area=Area_list
                    )


                if do_correction == True:
                    correction_count = Make_variable_corrections(
                        N = N
                        ,r = r_list[t+1] ,z = z_list[t+1], psi = psi_list[t+1] 
                        ,Area = Area_list
                        ,Area_init = np.sum(Area_list)
                        ,Tolerence = var_corr_tol
                        ,corr_max = 10
                        ,t = t
                    )
                    #correct_count_list[t+1] = correction_count       
            else:
                print("No integration method was choosen, program terminates")
                exit()




    fig,ax = plt.subplots()
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    ax.plot(
        time_vec[0:sim_steps],S_pot
        ,label="Epot"
        #,marker="."
        #,markersize=2
    )

    ax.plot(
        time_vec[0:sim_steps],S_pot_ref
        ,label="Epot ref"
        #,marker="."
        #,markersize=2
    )
    for m in range(len(t_points)):
        t_start = int(t_points[m])
        t_stop = int(t_start + num_of_free_steps + 20)
        if t_start != 0:
            t_start = int(t_points[m] - 20)
            t_stop = int(t_start + num_of_free_steps + 40)


        ymax = max(S_pot[int(t_start):t_stop])
        ymin = min(S_pot[int(t_start):t_stop])
        scaleing = ymax-ymin
        ymax = ymax + scaleing*0.3
        ymin = ymin - scaleing*0.3

        xleft = time_vec[t_start] - time_vec[sim_steps-1]/40
        xright = time_vec[t_start] + time_vec[sim_steps-1]/40

        line_width = 1
        ax.hlines(xmin=xleft,xmax=xright ,y=ymin ,color="k",linewidth=line_width)
        ax.hlines(xmin=xleft,xmax=xright,y=ymax ,color="k",linewidth=line_width)

        ax.vlines(x=xleft, ymin=ymin , ymax=ymax ,color="k",linewidth=line_width)
        ax.vlines(x=xright, ymin=ymin , ymax=ymax ,color="k",linewidth=line_width)

        ax.text(x=xright,y=ymax
                ,s=f"box:{m}")
    plt.title("Comparing the potential energy of the reference simulation \n "
              +f"and the loaded version w/o variable corrections"
              )
    
    plt.legend()
    plt.grid()

    """
    for t_start in t_points:
        t_stop = int(t_start + num_of_free_steps)
        
        fig,ax = plt.subplots()
        ax.plot(
            time_vec[0:sim_steps-1],S_pot
            ,label=f"Epot, tstart={t_start}"
            ,marker="."
        )

        plt.xlim(time_vec[int(t_start)],time_vec[t_stop])
        plt.ylim(
            min(S_pot[int(t_start):t_stop]) , max(S_pot[int(t_start):t_stop])
        )
        plt.legend()
        plt.grid()


    plt.show()"""
    #exit()
    fig,ax = plt.subplots(3,3)
    #fig.canvas.manager.full_screen_toggle()
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    #wm = plt.get_current_fig_manager()
    #wm.window.state('zoomed')
    m = 0
    for i in range(3):
        for j in range(3):
            t_start = int(t_points[m])
            t_stop = int(t_start + num_of_free_steps + 20)
            if t_start != 0:
                t_start = int(t_points[m] - 20)
                t_stop = int(t_start + num_of_free_steps + 40)
            

            ax[i,j].plot(
                    time_vec[int(t_start):t_stop],S_pot[int(t_start):t_stop]
                    ,label=f"Epot, box: {m}"
                    ,marker="."
                )
            
            ax[i,j].plot(
                    time_vec[int(t_start):t_stop],S_pot_ref[int(t_start):t_stop]
                    ,label="Epot reference"
                )
            ax[i,j].legend()
            ax[i,j].set_xlim([time_vec[t_start],time_vec[t_stop]])
            ax[i,j].grid()
            ymax = max(S_pot[int(t_start):t_stop])
            ymin = min(S_pot[int(t_start):t_stop])
            scaleing = ymax-ymin
            ax[i,j].set_ylim(
                ymin - scaleing*0.1 ,ymax + scaleing*0.1
            )
            if i == 2:
                ax[i,j].set_xlabel("time [s]")
            if j == 0:
                ax[i,j].set_ylabel("Epot [zJ]")

            m += 1

    plt.suptitle("Zoomed in parts of the total potential energy plot to get a close look")
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

    #test_if_constraint_diff_is_correct()
    #testing_new_epsilon_matrix()
    #testing_perturbation_function()
    #Testing_total_area_function()
    #testing_values_in_epsilon()
    #testing_initial_angles()
    #testing_if_constraints_are_true()
    #testing_arctan2_function()
    #testing_integration_with_events()
    #testing_for_no_correction_on_initial_state()

    #testing_gradient_of_constraints()
    #testing_gradient_for_S()

    #find_overflow_error()

    #test_convergence_of_alpha()

    #Testing_RungeKutta()


    #test_flat_model_object()
    test_gradients_again()
    #compare_potentential_energy()
    #test_if_variable_correction_causes_Epot_increase()