from Two_D_constants import Two_D_Constants, Two_D_paths
from Two_D_simulation_function import Two_D_simulation
from Make_movie import Make_frames, Make_video
from two_d_data_processing import tot_area, E_pot, E_kin, Xsqaured_test, Plot_3D_mesh
from Two_D_functions import Langrange_multi, Epsilon_values
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

    
def plot_tot_area(
        data_path:str,df_name:str,output_path:str
        ):

    if data_path == "" and df_name== "":
        print(f" No paths were given in the plot_tot_area function")
        exit()
    save_name = df_name.replace(".avi","")
    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    
    #print(df_sim.info())
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    dt = df_sim["dt"][0]
    sim_steps = df_sim["sim_steps"][0]
    corr_count = df_sim["correction count"][0]
    Xsqrt = df_sim["Chi squared test"][0]

    Area_change= np.zeros(sim_steps)
    time = np.zeros(sim_steps)
    for t in range(sim_steps):
        Area_change[t] = tot_area(N=N,r=r[t],z=z[t])
        time[t] = t*dt

    Amin, Amax = min(Area_change) ,max(Area_change)
    Aratio = Amax/Amin 

    """-------------------------------------- Total Area ----------------------------------------------------------"""
    fig, ax=plt.subplots()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.plot(time,Area_change,'-')
    plt.xlabel("time [s]",fontsize=15)
    plt.ylabel("total area",fontsize=15)
    plt.title(
        f"Ratio of Amax/Amin={Aratio} \n "
        +f"Amax - AMin={Amax-Amin}"
        ,fontsize=15)
    plt.grid()
    ax.ticklabel_format(useOffset=False)
    
    save_name_1 = save_name + "Atot"# Total area over time"
    save_name_1 = "Atot"
    #ax.set_aspect("equal",adjustable="box")
    plt.draw()
    plt.pause(2)
    plt.savefig(output_path + save_name_1 + ".png")
    plt.savefig(output_path + save_name_1 + ".svg")
    


    """-------------------------------------- Variable corrections pr time ----------------------------------------------------------"""
    fig,ax = plt.subplots()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.plot([i*dt for i in range(len(corr_count))],corr_count,linestyle="-")
    plt.title(
        f"correction counts pr time step \n "
        +r"$\frac{corrections}{time step} \approx$"
        +f"{round(np.sum(corr_count)/sim_steps ,3)}"
        ,fontsize=15
        )
    plt.xlabel("t[s]",fontsize=15)
    plt.ylabel("number of variables corrections",fontsize=15)
    plt.grid()
    
    save_name_1 = save_name + " var corr"# Total area over time"
    save_name_1 =  " var corr"
    #ax.set_aspect("equal",adjustable="box")
    plt.draw()
    plt.pause(2)
    
    plt.savefig(output_path + save_name_1 + ".png")
    plt.savefig(output_path + save_name_1 + ".svg")
    



    dA = np.zeros(sim_steps-1)
    for t in range(sim_steps-1):
        dA[t] = Area_change[t+1] - Area_change[t]

    """--------------------------------------Change in Area----------------------------------------------------------"""
    fig, ax = plt.subplots()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    #wm = plt.get_current_fig_manager()
    #wm.window.state('zoomed')
    plt.plot([i*dt for i in range(len(dA))],dA[0:sim_steps-1],'-')
    #ax.ticklabel_format(useOffset=False)
    plt.title("Change in Area",fontsize=15)
    plt.xlabel("t [s]",fontsize=15)
    plt.ylabel("Area",fontsize=15)
    plt.grid()
    save_name_2 = save_name + " dA"
    save_name_2 = " dA"
    #ax.set_aspect("equal",adjustable="box")
    plt.draw()
    plt.pause(2)
    plt.savefig(output_path + save_name_2 + ".png")
    plt.savefig(output_path + save_name_2 + ".svg")



    """--------------------------------------Chi----------------------------------------------------------"""
    fig,ax = plt.subplots()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.plot(time[0:sim_steps-1],Xsqrt,label=r"$\chi^2$ test")
    plt.title(r"$\chi^2$ test for deviation from the unperturbed state, so $\sigma_i$=1",fontsize=15)
    plt.xlabel("t [s]",fontsize=15)
    plt.ylabel(r"$\chi^2$ [$\mu m^2$]", fontsize=15)
    plt.grid()
    save_name_3 = save_name + "chisqrt"
    save_name_3 = "chisqrt"
    #ax.set_aspect("equal",adjustable="box")
    plt.draw()
    plt.pause(2)
    plt.savefig(output_path + save_name_3 + ".png")
    plt.savefig(output_path + save_name_3 + ".svg")




def plot_Epot_Ekin(
        data_path:str,df_name:str,output_path:str
        ):

    if data_path == "" or df_name== "" or output_path=="":
        print(f" No paths were given, in the plot_Epot_Ekin function")
        exit()
    save_name = df_name.replace(".avi","")

    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    #exit()
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    r_unperturbed = df_sim['r unperturbed'][0]
    z_unperturbed = df_sim['z unperturbed'][0]
    psi = df_sim['psi'][0]
    Area = df_sim['area list'][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    k = df_sim['k'][0]
    kG = df_sim['kG'][0]
    ds = df_sim['ds'][0]
    sigma = df_sim['sigma'][0]
    tau = df_sim['tau'][0]
    dt = df_sim["dt"][0]
    sim_steps = df_sim["sim_steps"][0]
    S = df_sim['Epot'][0]
    T = df_sim['Ekin'][0]
    font_size= 15
    

    t_vec = [dt*i for i in range(sim_steps-1)]
    fontsize = 10
    
    """--------------------------------------Kinetic & Potential energy plot ----------------------------------------------------------"""
    fig, ax = plt.subplots(2,1)
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')

    ax[0].plot(t_vec,T,"-",label="Kinetic energy")
    ax[0].set_xlabel("time [s]",fontsize=fontsize)
    ax[0].set_ylabel(r"$E_{kin} [zJ]$",fontsize=fontsize)
    ax[0].set_title("Kinetic energy, not scale properly \n" +f"min(E_kin)={min(T)}",fontsize=fontsize)
    ax[0].legend(fontsize=15)
    ax[0].ticklabel_format(useOffset=False)
    ax[0].grid()
    #plt.figure()
    ax[1].plot(t_vec,S,"-",label="Potential energy")
    ax[1].set_xlabel("time [s]",fontsize=fontsize)
    ax[1].set_ylabel(r"$E_{pot}$ [zJ]",fontsize=fontsize)
    ax[1].set_title("Potential energy  " +r"$min(E_{pot}) \approx$"+f"{round(min(S),3)}  and " +r"$\Delta E_{pot} \approx$"+f"{max(S)-min(S):0.1e}",fontsize=fontsize)
    ax[1].legend(fontsize=15)
    ax[1].ticklabel_format(useOffset=False)
    ax[1].grid()
    plt.draw()
    plt.pause(0.1)
    save_name_1 = "potential and kinetic energy" + df_name 
    save_name_1 = df_name +" S&T"
    save_name_1 = "S&T"
    #ax.set_aspect("equal",adjustable="box")
    plt.draw()
    plt.pause(2)
    plt.savefig(output_path + save_name_1 + ".png")
    plt.savefig(output_path + save_name_1 + ".svg")
    
    """--------------------------------------Kinetic energy plot stand alone----------------------------------------------------------"""
    fig, ax = plt.subplots()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.plot(t_vec,T,"-",label="Kinetic energy")
    plt.xlabel("time [s]",fontsize=15)
    plt.ylabel(r"$E_{kin}$ [zJ]",fontsize=15)
    plt.title("Kinetic Energy (Not correct scale)",fontsize=15)
    plt.ticklabel_format(useOffset=False)
    plt.legend(fontsize=15)
    plt.grid()
    save_name_4 = "Ekin"
    #ax.set_aspect("equal",adjustable="box")
    plt.draw()
    plt.pause(2)
    plt.savefig(output_path + save_name_4 + ".png")
    plt.savefig(output_path + save_name_4 + ".svg")
    
    """---------------------------------------Potential energy plot stand alone---------------------------------------------------------"""
    fig, ax = plt.subplots()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.plot(t_vec,S,"-",label="Potential energy")
    plt.xlabel("time [s]",fontsize=15)
    plt.ylabel(r"$E_{pot}$ [zJ]",fontsize=15)
    plt.title("Potential Energy",fontsize=15)
    plt.legend(fontsize=15)
    plt.ticklabel_format(useOffset=False)
    plt.grid()
    save_name_5 = "Epot"
    #ax.set_aspect("equal",adjustable="box")
    plt.draw()
    plt.pause(2)
    plt.savefig(output_path + save_name_5 + ".png")
    plt.savefig(output_path + save_name_5 + ".svg")


    """---------------------------------Inital positon and end positon plot---------------------------------------------------------------"""
    fig, ax = plt.subplots()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    if r_unperturbed[0] - r[0][0] == 0:
        pass
    else:
        plt.plot(
            r_unperturbed,z_unperturbed
            ,marker="v",linestyle="dotted"
            ,color="b"
            ,label="Unperturbed initial pos"
            )
    plt.plot(
        r[0],z[0]
        ,marker="o",linestyle="dashed"
        ,color="k"
        ,label="intial positions"
        )
    if r[sim_steps-1].all() == 0:
        plt.plot(
        r[sim_steps-1],z[sim_steps-1]
        ,marker="s",linestyle="solid"
        ,color="g"
        ,label="Error occured, and end isnt there"
        )
    else:
        plt.plot(
        r[sim_steps-1],z[sim_steps-1]
        ,marker="s",linestyle="solid"
        ,color="g"
        ,label="end positions"
        )
    plt.xlabel(r"r [$\mu m$]",fontsize=15)
    plt.ylabel(r"z [$\mu m$]",fontsize=15)
    plt.title("show difference from start and end positions",fontsize=15)
    plt.legend()#fontsize=font_size)
    plt.grid()
    save_name_2 = df_name + "init&end"
    save_name_2 =  "init&end"
    #ax.set_aspect("equal",adjustable="box")
    plt.draw()
    plt.pause(2)
    plt.savefig(output_path + save_name_2 +".png")
    plt.savefig(output_path + save_name_2 +".svg")


    """---------------------------------Inital positon and end positon scaled correcly plot---------------------------------------------------------------"""
    fig, ax = plt.subplots()
    font_size= 15
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    
    if r_unperturbed[0] - r[0][0]==0:
       pass
    else:
        plt.plot(
            r_unperturbed,z_unperturbed
            ,marker="v",linestyle="dotted"
            ,color="b"
            ,label="Unperturbed initial pos"
            )
    plt.plot(
        r[0],z[0]
        ,marker="o",linestyle="dashed"
        ,color="k"
        ,label="intial positions"
        )
    if r[sim_steps-1].all() == 0:
        plt.plot(
        r[sim_steps-1],z[sim_steps-1]
        ,marker="s",linestyle="solid"
        ,color="g"
        ,label="Error occured, and end isnt there"
        )
    else:
        plt.plot(
        r[sim_steps-1],z[sim_steps-1]
        ,marker="s",linestyle="solid"
        ,color="g"
        ,label="end positions"
        )

    plt.xlabel(r"r [$\mu m$]",fontsize=font_size)
    plt.ylabel(r"z [$\mu m$]",fontsize=font_size)
    xmin = min([min(r[t]) for t in range(sim_steps)]) - ds
    xmax  = r[0][N] + ds
    ymin = min([min(z[t]) for t in range(sim_steps)]) - ds
    ymax = ymin + (xmax - xmin)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.title("show difference from start and end positions",fontsize=15)
    plt.legend(fontsize=font_size)
    plt.grid()

    save_name_3 = df_name + "init&end scaled"
    save_name_3 =  "init&end scaled"
    ax.set_aspect("equal",adjustable="box")
    plt.draw()
    plt.pause(2)
    plt.savefig(output_path + save_name_3 +".png")
    plt.savefig(output_path + save_name_3 +".svg")

    """---------------------------------Show different positons a 5 different t---------------------------------------------------------------"""
    fig, ax = plt.subplots()
    font_size= 15
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')

    t_ref_list = [0, int(sim_steps/4), int(2*sim_steps/4), int(3*sim_steps/4), int(sim_steps)-1 ]
    style_list =["solid","dotted","dashed","dashdot","solid"]
    color_list = ["k","g","r","b","m","c"]
    marker_list= ["o","v","s","^","*",]
    for i in range(len(t_ref_list)):
        plt.plot(r[t_ref_list[i]],z[t_ref_list[i]]
        ,marker=marker_list[i]
        ,linestyle=style_list[i]
        ,color=color_list[i]
        ,label=r"t$\approx$"+f"{t_ref_list[i]*dt:0.1e}"
        )

    plt.xlabel(r"r [$\mu m$]",fontsize=font_size)
    plt.ylabel(r"z [$\mu m$]",fontsize=font_size)
    xmin = min([min(r[t]) for t in range(sim_steps)]) - ds
    xmax  = r[0][N] + ds
    ymin = min([min(z[t]) for t in range(sim_steps)]) - ds
    ymax = ymin + (xmax - xmin)


    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    ax.set_aspect("equal",adjustable="box")
    plt.title("show difference from start and end positions",fontsize=15)
    plt.legend(fontsize=font_size)
    plt.grid()
    save_name_3 = df_name + "init&end scaled"
    save_name_3 =  "4 different postions in time"
    plt.draw()
    plt.pause(2)
    plt.savefig(output_path + save_name_3 +".png")
    plt.savefig(output_path + save_name_3 +".svg")

    
    
    





def plot_reference_fig_for_finding_what_to_simulate():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib import cycler
    import os

    from two_d_continues_integration import find_init_stationary_state
    from Two_D_constants import Two_D_Constants_stationary_state


    np.set_printoptions(legacy='1.25')

    df_name = "Matlab data"
    data_path2 = "2D sim results\\"
    figure_save_path = data_path2 + "Data for thesis\\"


    df = pd.read_pickle(data_path2 + df_name)


    print(df.info())
    #print(len(df["tau_list_Flat"]))
    #print(len(df["tau_list_Curved"]))
    #print(df["tau_list_Curved"])


    max_tau = -1e20
    min_tau = 1e20
    diff_taus = []
    for n in range(len(df["tau_list_Curved"])):
        for i in range(len(df["tau_list_Curved"][n])):
            if df["tau_list_Curved"][n][i] < min_tau:
                min_tau = df["tau_list_Curved"][n][i]
            
            if df["tau_list_Curved"][n][i] > max_tau:
                max_tau = df["tau_list_Curved"][n][i]

            if df["tau_list_Curved"][n][i] not in diff_taus:
                diff_taus.append(df["tau_list_Curved"][n][i])



    for n in range(len(df["tau_list_Flat"])):
        for i in range(len(df["tau_list_Flat"][n])):
            if df["tau_list_Flat"][n][i] < min_tau:
                min_tau = df["tau_list_Flat"][n][i]
            
            if df["tau_list_Flat"][n][i] > max_tau:
                max_tau = df["tau_list_Flat"][n][i]
            
            if df["tau_list_Flat"][n][i] not in diff_taus:
                diff_taus.append(df["tau_list_Flat"][n][i])


    #pos_A_min_x1 ,pos_A_min_x2 ,neg_A_min_x3 = 36 ,3.4 ,-1.08#-12
    #pos_A_max_x1 ,pos_A_max_x2 ,neg_A_max_x3 = 37 ,3.6 ,-1.06#-11
    #pos_r1_min_y1 ,pos_r1_min_y2 ,neg_r1_min_y3 = 3 ,7.35 , 0.64#2.725
    #pos_r1_max_y1 ,pos_r1_max_y2 ,neg_r1_max_y3 = 3.05 ,7.45 ,0.66#2.8
    
    """ Triangle """
    pos_A_min_x1,pos_A_max_x1, pos_r1_min_y1, pos_r1_max_y1 = 23.6 ,23.8 ,2.83 ,2.84#36 ,37 ,3 ,3.05 # 12.0 ,12.05,1.81,1.815 #
    """ plus """
    pos_A_min_x2,pos_A_max_x2, pos_r1_min_y2, pos_r1_max_y2 = 3.4 ,3.6 ,7.35 ,7.45
    """ cross """
    neg_A_min_x3,neg_A_max_x3,neg_r1_min_y3,neg_r1_max_y3 = -1.08 ,-1.06 ,0.64 ,0.66
    n_neg_A ,n_pos_A = [] ,[]
    i_neg ,i_pos = [],[]

    for n in range(len(df["ExcessAreaCurved"])):
        for i in range(len(df["ExcessAreaCurved"][n])):
            if len(df["ExcessAreaCurved"]) > 0:
                A_curve = df["ExcessAreaCurved"][n][i]
                r1_curve = df["r1Curved"][n][i]
                if pos_A_min_x1 < A_curve < pos_A_max_x1 and pos_r1_min_y1 < r1_curve < pos_r1_max_y1:
                    n_pos_A.append(n)
                    i_pos.append(i)
                    print("first")
                if pos_A_min_x2 < A_curve < pos_A_max_x2 and pos_r1_min_y2 < r1_curve < pos_r1_max_y2:
                    n_pos_A.append(n)
                    i_pos.append(i)
                    print("2nd")

    for n in range(len(df["ExcessAreaFlat"])):
        for i in range(len(df["ExcessAreaFlat"][n])):
            if len(df["ExcessAreaFlat"]) > 0:
                A_Flat = df["ExcessAreaFlat"][n][i]
                r1_Flat = df["r1Flat"][n][i]
                if neg_A_min_x3 < A_Flat < neg_A_max_x3 and neg_r1_min_y3 < r1_Flat < neg_r1_max_y3:
                    n_neg_A.append(n)
                    i_neg.append(i)

    print(f"pos n={n_pos_A} i pos={i_pos}")
    print(f"pos n={n_neg_A} i pos={i_neg}")

    #print(diff_taus)
    #print(len(diff_taus))
    print(f"mintau={min_tau} and max tau = {max_tau}")

    """----------------------------------- Matlab data ----------------------------------------------"""
    const_args = Two_D_Constants_stationary_state(
            print_val = False
            ,show_stationary_state = False
            ,start_flat = False
            ,perturb = False
        )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0= const_args[6:8]
    #sigma, tau = const_args[9:11]

    lc = 1/c0
    sigma_c = k*c0**2
    tau_c = k*c0

    rs2 = 20*lc 
    zs2 = 0
    s0, sN = 0, 50*lc

    fig, ax = plt.subplots(1,2)
    cmap = plt.cm.coolwarm
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')

    index = [i for i in range(len(diff_taus))]
    i = 0
    plot_tau_ref = []

    for n in range(len(df["ExcessAreaCurved"])):
        if len(df["tau_list_Curved"][n]) > 0 :
            i = (df["tau_list_Curved"][n][0]-1)/(max_tau - 1)
            if df["tau_list_Curved"][n][0] not in plot_tau_ref:
                #plt.plot(df["ExcessAreaCurved"][n],df["NeckRadiusCurved"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Curved"][n][0]:0.1f}")
                ax[0].plot(
                    [i/lc**2 for i in df["ExcessAreaCurved"][n]]
                    ,[i/lc for i in df["r1Curved"][n]],".-",color=cmap(i)
                    ,label=r"$\tau \approx$"+f"{(df["tau_list_Curved"][n][0]/lc**2)/1000:0.1f} nN"
                    )
                plot_tau_ref.append(df["tau_list_Curved"][n][0])
            else:
                #plt.plot(df["ExcessAreaCurved"][n],df["NeckRadiusCurved"][n],".-")#,color=cmap(i))
                ax[0].plot(
                    [i/lc**2 for i in df["ExcessAreaCurved"][n]]
                    , [i/lc for i in df["r1Curved"][n]]
                    ,".-",color=cmap(i)
                    )#,color=cmap(i))


    for n in range(len(df["ExcessAreaFlat"])):
        if len(df["tau_list_Flat"][n]) > 0:
            i = (df["tau_list_Flat"][n][0]-1)/(max_tau -1)
            if df["tau_list_Flat"][n][0] not in plot_tau_ref:
                #plt.plot(df["ExcessAreaFlat"][n],df["NeckRadiusFlat"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Flat"][n][0]:0.1f}")
                ax[0].plot(
                    [i/lc**2 for i in df["ExcessAreaFlat"][n]]
                    ,[i/lc for i in df["r1Flat"][n]]
                    ,".-",color=cmap(i),label=r"$\tau \approx$"+f"{(df["tau_list_Flat"][n][0]/lc**2)/1000:0.1f} nN"
                    )
                plot_tau_ref.append(df["tau_list_Flat"][n][0])
            else:
                #plt.plot(df["ExcessAreaFlat"][n],df["NeckRadiusFlat"][n],".-")
                ax[0].plot(
                    [i/lc**2 for i in df["ExcessAreaFlat"][n]]
                    ,[i/lc for i in df["r1Flat"][n]]
                    ,".-",color=cmap(i))
        

    markers = ["^","+","x"]
    markers_latex = [r"$\triangle$",r"$\plus$" ,r"$\times$"]
    ax[0].plot(
        df["ExcessAreaCurved"][n_pos_A[0]][i_pos[0]]/lc**2
        ,df["r1Curved"][n_pos_A[0]][i_pos[0]]/lc
        ,marker=markers[0] 
        ,color="k"
        ,markersize = 15
        ,mfc = "none"
        )

    ax[0].plot(
        df["ExcessAreaCurved"][n_pos_A[1]][i_pos[1]]/lc**2
        ,df["r1Curved"][n_pos_A[1]][i_pos[1]]/lc
        ,marker=markers[1] 
        ,color="k"
        ,markersize = 12
        )
    ax[0].plot(
        df["ExcessAreaFlat"][n_neg_A[0]][i_neg[0]]/lc**2
        ,df["r1Flat"][n_neg_A[0]][i_neg[0]]/lc
        ,marker=markers[2] 
        ,color="k"
        ,markersize = 12
        ,mfc = "none"
        )

    plt.subplots_adjust(
        #left=0.05
        #,right=0.5
        wspace=0.6
    )
    #plt.legend(fontsize=13)
    ax[0].legend(
        bbox_to_anchor=(1.02, 1)
        , loc='upper left'
        , borderaxespad=0
        ,fontsize=15
                )
    ax[0].vlines(x=0,ymin=-100,ymax=300,colors="k",linestyles="--")
    #plt.xlim(-60, 50)
    ax[0].set_xlim(-1.2e4, 2.9e4)
    ax[0].set_ylim(0 ,230)
    ax[0].set_title(
        "Edge radius (r1) vs Excess Area"
        #+r"[$\tau]=\frac{\mu g\cdot \mu m }{s^2}$"
        ,fontsize=15)
    ax[0].set_xlabel(r"$\Delta A_{Excess} = A_{membrane} - A_{disc} $  [$\mu m^2$]",fontsize=15)
    ax[0].set_ylabel(r"Edge radius (r1) [$\mu m$]",fontsize=15)
    ax[0].grid()
    #ax[0].set_aspect("equal",adjustable="box")


    """---------------------------------------------- Initial configurations ---------------------------------------------------"""
    #fig,ax = plt.subplots()
    """lc = 1/c0
    sigma_c = k*c0**2
    tau_c = k*c0

    rs2 = 20*lc 
    zs2 = 0
    s0, sN = 0, 50*lc"""
    print(f"n pos A :{n_pos_A} and i pos = {i_pos}")
    print(f"n neg A :{n_neg_A} and i neg = {i_neg}")

    sigma_list = [
        df["sigma_list_Curved"][n_pos_A[0]][i_pos[0]]
        ,df["sigma_list_Curved"][n_pos_A[1]][i_pos[1]]
        ,df["sigma_list_Flat"][n_neg_A[0]][i_neg[0]]
        ]
    tau_list = [
        df["tau_list_Curved"][n_pos_A[0]][i_pos[0]]
        ,df["tau_list_Curved"][n_pos_A[1]][i_pos[1]]
        ,df["tau_list_Flat"][n_neg_A[0]][i_neg[0]]
        ]
    psi_L_list = [
        df["psi_L_list_Curved"][n_pos_A[0]][i_pos[0]]
        ,df["psi_L_list_Curved"][n_pos_A[1]][i_pos[1]]
        ,df["psi_L_list_Flat"][n_neg_A[0]][i_neg[0]]
        ]

    
    print("first integration")
    psi,r,z, r_contin_0, z_contin_0, alpha = find_init_stationary_state(
                k=k ,c0=c0 ,ds=ds
                ,r_L=rs2 ,z_L=zs2 ,s0=s0 ,sN=sN
                ,total_points = N
                ,tau=tau_list[0]*tau_c
                ,sigma=sigma_list[0]*sigma_c
                ,psi_L=psi_L_list[0]
            )
    print("second integration")
    psi,r,z, r_contin_1, z_contin_1, alpha = find_init_stationary_state(
                k=k ,c0=c0 ,ds=ds
                ,r_L=rs2 ,z_L=zs2 ,s0=s0 ,sN=sN
                ,total_points = N
                ,tau=tau_list[1]*tau_c
                ,sigma=sigma_list[1]*sigma_c
                ,psi_L=psi_L_list[1]
            )
    print("third integration")
    psi,r,z, r_contin_2, z_contin_2, alpha = find_init_stationary_state(
                k=k ,c0=c0 ,ds=ds
                ,r_L=rs2 ,z_L=zs2 ,s0=s0 ,sN=sN
                ,total_points = N
                ,tau=tau_list[2]*tau_c
                ,sigma=sigma_list[2]*sigma_c
                ,psi_L=psi_L_list[2]
            )

    spacing = 0.15 #max(z_contin_2) + max(z_contin_1) + max(z_contin_0)

    ax[1].plot(
        r_contin_0
        ,z_contin_0 + spacing
        ,label = markers_latex[0] 
        + r": $ r1 \approx $" + f"{df["r1Curved"][n_pos_A[0]][i_pos[0]]/lc:.1f}" +r" $\mu m$"
        + r", $\tau \approx $" + f"{(df["tau_list_Curved"][n_pos_A[0]][i_pos[0]]/lc**2)/1000:0.1f} nN"
        ,linestyle="-"
        ,linewidth=2
            )

    ax[1].plot(
        r_contin_1
        ,z_contin_1 
        ,label = markers_latex[1] 
        + r": $ r1 \approx $" + f"{df["r1Curved"][n_pos_A[1]][i_pos[1]]/lc:.1f}" +r" $\mu m$"
        + r", $\tau \approx $" + f"{(df["tau_list_Curved"][n_pos_A[1]][i_pos[1]]/lc**2)/1000:0.1f} nN"
        ,linestyle="dashed"
        ,linewidth=2
        )

    ax[1].plot(
        r_contin_2
        ,z_contin_2 - spacing
        ,label = markers_latex[2]  
        + r": $ r1  \approx $" + f" {df["r1Flat"][n_neg_A[0]][i_neg[0]]/lc:.1f}" +r" $\mu m$"
        + r", $\tau \approx $" + f"{(df["tau_list_Flat"][n_neg_A[0]][i_neg[0]]/lc**2)/1000:0.1f} nN"
        ,linestyle="-."
        ,linewidth=2
        )

    rmin = 0#min([min(r_contin_0) , min(r_contin_1) , min(r_contin_2)]) - ds
    rmax = max([max(r_contin_0) , max(r_contin_1) , max(r_contin_2)]) #+ ds

    ax[1].set_yticks([i for i in np.arange(-3,3,0.05)])
    ax[1].set_xticks([i for i in np.arange(-3,3,0.05)])
    ax[1].set_xlim(rmin,rmax +0.03*(rmax-rmin))
    lim_spacing = -1.5*spacing
    ax[1].set_ylim(lim_spacing,(rmax-rmin) + lim_spacing)
    #ax[1].set_ylim(-(rmax-rmin)/2,(rmax-rmin)/2)
    #ax[1].set_ylim(-0.15,0.2)
    ax[1].legend(fontsize=20)
    ax[1].set_xlabel(r"r [$\mu m$]",fontsize=15)
    ax[1].set_ylabel(r"z [$\mu m$]",fontsize=15)
    ax[1].set_title(
        "Initial configurations of the three chosen \n points in the phase space diagram on the left"
        ,fontsize=15
        )

    ax[1].set_aspect("equal",adjustable="box")
    ax[1].grid()
    plt.draw()
    plt.pause(1)
    plt.savefig(figure_save_path + "edge_radius_ExcessArea.png")
    plt.savefig(figure_save_path + "edge_radius_ExcessArea.svg")


    for i in range(len(sigma_list)):
        print(
            f"   symbol:  "+ markers[i] + "\n"
            ,f"sigma ={sigma_list[i]}    "
            ,f"tau ={tau_list[i]}   "
            ,f"psi2 ={psi_L_list[i]}    \n"
            ,f"sigma ={sigma_list[i]*sigma_c}    "
            ,f"tau ={tau_list[i]*tau_c}   "
        )
        print("---------------------------------")

    plt.show()





def Find_the_circle_radius_of_rolling_test():
    from two_d_data_processing import cirle_fit, make_circle, make_circle_V2

    data_path = "2D sim results\\Data for thesis\\Verification\\c0=c0 tau=0\\"
    #data_path = "C:\\Users\\adams\\Documents\\GitHub\\Masters-Project-BioPhysics\\Surface model\\2D sim results\\Data for thesis\\Verification\\c0=c0 tau=0\\"
    data_name = "2D surface N,ds,dt,T,tau,c0=(40, 0.015, 1e-11, 1e-06, 0, 25).pkl"

    df_sim = pd.read_pickle(data_path+data_name)

    print("I have read the file the first time around")
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    r_unperturbed = df_sim['r unperturbed'][0]
    z_unperturbed = df_sim['z unperturbed'][0]
    psi = df_sim['psi'][0]
    Area = df_sim['area list'][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    k = df_sim['k'][0]
    kG = df_sim['kG'][0]
    ds = df_sim['ds'][0]
    sigma = df_sim['sigma'][0]
    tau = df_sim['tau'][0]
    dt = df_sim["dt"][0]
    sim_steps = df_sim["sim_steps"][0]
    S = df_sim['Epot'][0]
    T = df_sim['Ekin'][0]


    
    start_point = 0
    end_point = 22
    rc1, zc1, R = cirle_fit(
        data_path=data_path
        ,df_name=data_name
        ,edge_point = start_point
        ,end_point= end_point
    )
    print("I have fitted a circle")

    """r_circ, z_circ = make_circle(
        xc=xc1,zc=zc1,R=R,ds=1e-5
        ,xlim=[min(r[sim_steps-1][start_point:end_point]) , max(r[sim_steps-1][start_point:end_point])]
        ,zlim=[min(z[sim_steps-1][start_point:end_point]) , max(z[sim_steps-1][start_point:end_point])]
        )"""
    r_circ_up,r_circ_down, z_circ_up, z_circ_down = make_circle_V2(
        rc=rc1,zc=zc1,R=R
        ,rmax= max(r[sim_steps-1][start_point:end_point])
        ,rmin=min(r[sim_steps-1][start_point:end_point])
        ,zmax=max(z[sim_steps-1][start_point:end_point])
        ,zmin=min(z[sim_steps-1][start_point:end_point])
        ,step_size=1e-7
    )
    print("I have made the circle")
    """------------------------------------------------------------------------------------------------"""
    fig, ax = plt.subplots()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')

    plt.plot(r_circ_up,z_circ_up
             ,linestyle="-"
             ,color="k"
             ,label=f"circle fit"#+r"r$\approx$"+f"{R:0.3f}"
             ,linewidth = 3
             )
    
    plt.plot(r_circ_down,z_circ_down
             ,linestyle="-"
             ,color="k"
             ,linewidth = 3
             #,label=f"circle fit "+r"r$\approx$"+f"{R:0.3f}"
             )
    
    plt.plot(
        [rc1    ,r_circ_up[int( len(z_circ_down)/1.2)] ]
        ,[zc1  ,z_circ_up[int( len(z_circ_down)/1.2 )]]
        ,label=f"radius"+r"$\approx$"+f"{R:0.3f}" +r"$\mu m$"
        ,color="k"
        ,linestyle="-"
        ,marker=""
        )
    
    plt.plot(r[sim_steps-1][end_point],z[sim_steps-1][end_point]
             ,marker="v"
             ,linestyle=""
             ,color="b"
             #,label="stop point"
             ,markersize = 15
             )

    plt.plot(r[sim_steps-1],z[sim_steps-1]
             ,marker="*"
             ,linestyle="-"
             ,color="m"
             ,label=f"membrane t={(sim_steps-1)*dt:0.1e}s"
             #,linewidth= 5
             )
    
    
  
    
    rmin = min(r[sim_steps-1][start_point:end_point]) - ds
    rmax =  max(r[sim_steps-1][start_point:end_point]) + ds
    deltar = rmax - rmin
    zmin = 0
    zmax = deltar
    plt.xlim(rmin,rmax)
    plt.ylim(zmin,zmax)
    ax.set_aspect("equal",adjustable="box")
    plt.legend(fontsize=15)
    plt.title("")
    plt.grid()
    plt.xlabel(r"r [$\mu m$]",fontsize=15)
    plt.ylabel(r"z [$\mu m$]",fontsize=15)
    plt.draw()
    plt.pause(3)
    plt.savefig(data_path + "Rolling circle fit test.png")
    plt.savefig(data_path + "Rolling circle fit test.svg")
    plt.show()


def Investigating_chosen_configuration_1():

    data_path = "2D sim results\\Data for thesis\\Data simulation\\"
    folder_names = " sigma,tau,psi2=(1.8e+042.0e+03-1.5e-08)\\"
    save_name_list = ["-perturbed","+perturbed","unperturbed"]
    data_name_list = [
        "2D  N,ds,dt,T,tau,c0=(20, 0.015, 1.25e-13, 1e-07, 2000, 25)"
        ,"2D  N,ds,dt,T,tau,c0=(20, 0.015, 1.25e-13, 1e-07, 2000, 25)"
        ,"2D  N,ds,dt,T,tau,c0=(20, 0.015, 1.25e-13, 1e-08, 2000, 25)"
        ]

    df_sim_minus_perturb = pd.read_pickle(data_path + save_name_list[0] + folder_names +  data_name_list[0])
    df_sim_plus_perturb = pd.read_pickle(data_path + save_name_list[1] + folder_names +  data_name_list[1])
    df_sim_unperturb = pd.read_pickle(data_path + save_name_list[2] + folder_names +  data_name_list[2])


    #print(df_sim_minus_perturb.info())
    

    """----------------------------------------- Minus perturbed data-----------------------------------------------------"""
    X2_minus_perturb = df_sim_minus_perturb["Chi squared test"][0]
    sim_steps_minus_perturb = df_sim_minus_perturb["sim_steps"][0]
    dt_minus_perturb = df_sim_minus_perturb["dt"][0]
    r_minus_perturb = df_sim_minus_perturb["r"][0]
    z_minus_perturb = df_sim_minus_perturb["z"][0]
    time_minus_perturb = np.linspace(0,sim_steps_minus_perturb*dt_minus_perturb,sim_steps_minus_perturb-1)


    """----------------------------------------- Plus perturbed data-----------------------------------------------------"""
    X2_plus_perturb = df_sim_plus_perturb["Chi squared test"][0]
    sim_steps_plus_perturb = df_sim_plus_perturb["sim_steps"][0]
    dt_plus_perturb = df_sim_plus_perturb["dt"][0]
    r_plus_perturb = df_sim_plus_perturb["r"][0]
    z_plus_perturb = df_sim_plus_perturb["z"][0]
    time_plus_perturb = np.linspace(0,sim_steps_plus_perturb*dt_plus_perturb,sim_steps_plus_perturb-1)


    """----------------------------------------- Unperturbed data-----------------------------------------------------"""
    X2_unperturb = df_sim_unperturb["Chi squared test"][0]
    sim_steps_unperturb = df_sim_unperturb["sim_steps"][0]
    dt_unperturb = df_sim_unperturb["dt"][0]
    r_unperturb = df_sim_unperturb["r"][0]
    z_unperturb = df_sim_unperturb["z"][0]
    r_unperturb_init = df_sim_unperturb["r unperturbed"][0]
    z_unperturb_init = df_sim_unperturb["z unperturbed"][0]
    ds = df_sim_unperturb["ds"][0]
    time_unperturb = np.linspace(0,sim_steps_unperturb*dt_unperturb,sim_steps_unperturb-1)

    

    """----------------------------------------- Compare all X^2 results -----------------------------------------------------"""
    fig, ax = plt.subplots()
    #wm = plt.get_current_fig_manager()
    #wm.window.state('zoomed')
    
    line_width = 2
    plt.plot(
        time_minus_perturb,X2_minus_perturb
        ,label="-perturb"
        ,linestyle="-."
        ,linewidth=line_width
        )
    plt.plot(
        time_plus_perturb,X2_plus_perturb
        ,label="+perturb"
        ,linestyle="--"
        ,linewidth=line_width
        )
    plt.plot(
        time_unperturb,X2_unperturb
        ,label="unperturb"
        ,linestyle="-"
        ,linewidth=line_width
        )

    ymax = max([
         max(X2_minus_perturb[0:int(1e-8/dt_minus_perturb)])
        ,max(X2_plus_perturb[0:int(1e-8/dt_plus_perturb)])
        ,max(X2_unperturb[0:int(1e-8/dt_unperturb)])
        ])
    xmax = sim_steps_unperturb*dt_unperturb*1.01 
    xmin = -1e-10
    ymin = -0.3e-4
    plt.xlim(xmin, xmax)
    plt.ylim(ymin,ymax)
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.legend(fontsize=12)
    plt.xlabel("t [s]",fontsize=15)
    plt.ylabel(r"$\chi^2$ [$\mu m^2$]",fontsize=15)
    plt.grid()
    plt.draw()
    plt.pause(2)

    save_name_1 = "chisqrt for all 3"
    plt.savefig(data_path + save_name_1 +".png")
    plt.savefig(data_path + save_name_1 +".svg")


    """----------------------------------------- Compare all X^2 results -----------------------------------------------------"""    
    fig, ax = plt.subplots()
    ax.set_aspect("equal",adjustable="box")
    #wm = plt.get_current_fig_manager()
    #wm.window.state('zoomed')
    line_width = 2.5
    marker_size = 10
    plt.plot(
        r_unperturb_init,z_unperturb_init
        ,marker="o"
        ,linestyle="-"
        ,label="unperturbed initial pos"
        ,linewidth=line_width
        ,markersize=marker_size
        )
    
    plt.plot(
        r_minus_perturb[sim_steps_minus_perturb-1],z_minus_perturb[sim_steps_minus_perturb-1]
        ,marker="s"
        ,linestyle="--"
        ,label="-perturb"
        ,linewidth=line_width
        ,markersize=marker_size
        )
    
    plt.plot(
        r_plus_perturb[sim_steps_unperturb-1] ,z_plus_perturb[sim_steps_unperturb-1]
        ,marker="^"
        ,linestyle="-."
        ,label="+perturb"
        ,linewidth=line_width
        ,markersize=marker_size
        )
    
    plt.plot(
        r_unperturb[sim_steps_unperturb-1],z_unperturb[sim_steps_unperturb-1]
        ,marker="*"
        ,linestyle="dotted"
        ,label="unperturb"
        ,linewidth=line_width
        ,markersize=marker_size
        )
    
    
    xmax = max(r_unperturb_init) #+ ds
    xmin =  min([
        min(r_unperturb_init)
        ,min(r_minus_perturb[sim_steps_minus_perturb-1])
        ,min(r_plus_perturb[sim_steps_unperturb-1])
        ,min(r_unperturb[sim_steps_unperturb-1])
    ])
    ymin =  min([
        min(z_unperturb_init)
        ,min(z_minus_perturb[sim_steps_minus_perturb-1])
        ,min(z_plus_perturb[sim_steps_unperturb-1])
        ,min(z_unperturb[sim_steps_unperturb-1])
    ])
    ymax = xmax - xmin 
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.xlabel(r"r [$\mu m$]",fontsize=15)
    plt.ylabel(r"z [$\mu m$]",fontsize=15)
    plt.legend(fontsize=15)
    plt.grid()
    plt.show()



def figure_3_potential_energy_landscape_cases():
    import numpy as np
    import matplotlib.pyplot as plt
    x_span = 2.3
    x_vec = np.linspace(start=-x_span,stop=x_span,num=100)
    a = -1
    b = 3
    c = 0#-6
    d = 0#2
    x_shift = -1
    Y3 = [a*(x - x_shift)**3 + b*(x - x_shift)**2 + c*(x - x_shift) + d for x in x_vec]

    fig, ax = plt.subplots()

    """------------------- 3rd depgree poly plot------------------------------------"""

    plt.plot(x_vec,Y3,label="3rd")

    plt.hlines(y=0,xmin=-1.5,xmax=4,linestyles="--",colors="k")
    plt.hlines(y=4,xmin=-1.5,xmax=4,linestyles="--",colors="k")

    xycoords = 'figure fraction'
    ax.annotate("",xytext=(3,0),xy=(3,4),arrowprops=dict(arrowstyle="<->"))
    ax.text(x=3.2,y=2,s=r"$\Delta y$",fontsize=15)


    """------------------- first 2rd depgree poly plot------------------------------------"""
    
    x_span = 2.3
    x_vec = np.linspace(start=-x_span,stop=x_span,num=100)
    a = 0-1
    b = -3
    c = 0#-6
    d = 0#2
    x_shift = -1
    Y2_neg = [a*(x - x_shift)**3 + b*(x - x_shift)**2 + c*(x - x_shift) + d for x in x_vec]

    plt.plot(x_vec,Y2_neg)

    plt.ylim(-2,10)
    plt.xlim(-x_span, x_span +4)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.legend()
    plt.grid()
    plt.show()




def Investigating_chosen_configuration_New_data():

    data_path ="2D sim results\\Data for thesis\\multi processor result\\" + "cross sims\\T,dt,sigma,tau=(1e-07, 1.25e-13,5.8e+03,2.0e+03)\\"
    #folder_names = "cross sims\\T,dt,sigma,tau=(1e-07, 1.25e-13,5.8e+03,2.0e+03)\\"
    save_name_list = ["-perturbed\\","+perturbed\\","unperturbed\\"]
    data_name = "2D Surface.pkl"
    #2D sim results\Data for thesis\multi processor result\cross sims\T,dt,sigma,tau=(1e-07, 1.25e-13,5.8e+03,2.0e+03)\+perturbed
    df_sim_minus_perturb = pd.read_pickle(data_path  +save_name_list[0]+  data_name)
    df_sim_plus_perturb = pd.read_pickle(data_path + save_name_list[1]  + data_name)
    df_sim_unperturb = pd.read_pickle(data_path + save_name_list[2]  +  data_name)


    #print(df_sim_minus_perturb.info())
    

    """----------------------------------------- Minus perturbed data-----------------------------------------------------"""
    X2_minus_perturb = df_sim_minus_perturb["Chi squared test"][0]
    sim_steps_minus_perturb = df_sim_minus_perturb["sim_steps"][0]
    dt_minus_perturb = df_sim_minus_perturb["dt"][0]
    r_minus_perturb = df_sim_minus_perturb["r"][0]
    z_minus_perturb = df_sim_minus_perturb["z"][0]
    time_minus_perturb = np.linspace(0,sim_steps_minus_perturb*dt_minus_perturb,sim_steps_minus_perturb-1)


    """----------------------------------------- Plus perturbed data-----------------------------------------------------"""
    X2_plus_perturb = df_sim_plus_perturb["Chi squared test"][0]
    sim_steps_plus_perturb = df_sim_plus_perturb["sim_steps"][0]
    dt_plus_perturb = df_sim_plus_perturb["dt"][0]
    r_plus_perturb = df_sim_plus_perturb["r"][0]
    z_plus_perturb = df_sim_plus_perturb["z"][0]
    time_plus_perturb = np.linspace(0,sim_steps_plus_perturb*dt_plus_perturb,sim_steps_plus_perturb-1)


    """----------------------------------------- Unperturbed data-----------------------------------------------------"""
    X2_unperturb = df_sim_unperturb["Chi squared test"][0]
    sim_steps_unperturb = df_sim_unperturb["sim_steps"][0]
    dt_unperturb = df_sim_unperturb["dt"][0]
    r_unperturb = df_sim_unperturb["r"][0]
    z_unperturb = df_sim_unperturb["z"][0]
    r_unperturb_init = df_sim_unperturb["r unperturbed"][0]
    z_unperturb_init = df_sim_unperturb["z unperturbed"][0]
    ds = df_sim_unperturb["ds"][0]
    N = df_sim_unperturb["N"][0]
    time_unperturb = np.linspace(0,sim_steps_unperturb*dt_unperturb,sim_steps_unperturb-1)
    sigma = df_sim_unperturb["sigma"][0]
    tau = df_sim_unperturb["tau"][0]
    

    """----------------------------------------- Compare all X^2 results -----------------------------------------------------"""
    fig, ax = plt.subplots()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    
    line_width = 2
    plt.plot(
        time_minus_perturb,X2_minus_perturb
        ,label="-perturb"
        ,linestyle="-."
        ,linewidth=line_width
        )
    plt.plot(
        time_plus_perturb,X2_plus_perturb
        ,label="+perturb"
        ,linestyle="--"
        ,linewidth=line_width
        )
    plt.plot(
        time_unperturb,X2_unperturb
        ,label="unperturb"
        ,linestyle="-"
        ,linewidth=line_width
        )

    ymax = max([
         max(X2_minus_perturb)
        ,max(X2_plus_perturb)
        ,max(X2_unperturb)
        ])*1.05
    
    #xmax = sim_steps_unperturb*dt_unperturb*1.01 
    #xmin = min([-dt_minus_perturb,-dt_plus_perturb,-dt_unperturb])/1e3
    ymin = min([
         min(X2_minus_perturb)
        ,min(X2_plus_perturb)
        ,min(X2_unperturb)
        ])*(1-0.05)

    #plt.xlim(xmin, xmax)
    #plt.ylim(ymin,ymax)
    plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    plt.legend(fontsize=12)
    plt.xlabel("t [s]",fontsize=15)
    plt.ylabel(r"$\chi^2$ [$\mu m^2$]",fontsize=15)
    plt.title(
        r"$\chi^2$ test for the sum of the difference of each of the points along the links"
        +f"\n"
        +r"$\sigma \approx $" +f"{sigma/1000:0.1f}" + r" $nN/\mu m$ and $\tau \approx$" +f"{tau/1000:0.1f} nN"
        ,fontsize=15)
    plt.grid()
    plt.draw()
    plt.pause(2)

    save_name_1 = "chisqrt for all 3"
    plt.savefig(data_path + save_name_1 +".png")
    plt.savefig(data_path + save_name_1 +".svg")
    
    
    """----------------------------------------- chi^2 test for the end of unperturbed with the others -----------------------------------------------------"""    
    fig, ax = plt.subplots()
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    X2_minus_perturbed_2 = np.zeros(sim_steps_unperturb-1)
    X2_plus_perturbed_2 = np.zeros(sim_steps_plus_perturb-1)

    for t in range(sim_steps_minus_perturb-1):
        X2_minus_perturbed_2[t] = Xsqaured_test(
            N=N
            ,r_init=r_unperturb[sim_steps_unperturb-1]
            ,z_init=z_unperturb[sim_steps_unperturb-1]
            ,r=r_minus_perturb[t]
            ,z=z_minus_perturb[t]
            )
        
    for t in range(sim_steps_plus_perturb-1):
        X2_plus_perturbed_2[t] = Xsqaured_test(
            N=N
            ,r_init=r_unperturb[sim_steps_unperturb-1]
            ,z_init=z_unperturb[sim_steps_unperturb-1]
            ,r=r_plus_perturb[t]
            ,z=z_plus_perturb[t]
            )

    plt.plot(
        time_minus_perturb#[dt_minus_perturb*t for t in range(sim_steps_minus_perturb)]
        ,X2_minus_perturbed_2
        ,label="-perturbed"
        )
    plt.plot(
        time_plus_perturb#[dt_plus_perturb*t for t in range(sim_steps_plus_perturb)]
        ,X2_plus_perturbed_2
        ,label="+perturbed"
        )

    plt.hlines(
        y=0,xmin=-1,xmax=1#max(time_minus_perturb)
        ,label="target"
        ,linestyles="dashed"
        ,color = "black"
        )

    xmax = max([max(time_minus_perturb), max(time_plus_perturb),max(time_unperturb)])
    plt.xlim(xmax*(-5e-3) ,xmax*1.05)
    plt.ylim()
    plt.grid()
    plt.legend()
    plt.show()
    exit()
    """----------------------------------------- Compare all final positions results -----------------------------------------------------"""    
    fig, ax = plt.subplots()
    ax.set_aspect("equal",adjustable="box")
    #wm = plt.get_current_fig_manager()
    #wm.window.state('zoomed')
    line_width = 2.5
    marker_size = 10
    plt.plot(
        r_unperturb_init,z_unperturb_init
        ,marker="o"
        ,linestyle="-"
        ,label="unperturbed initial pos"
        ,linewidth=line_width
        ,markersize=marker_size
        )
    
    plt.plot(
        r_minus_perturb[sim_steps_minus_perturb-1],z_minus_perturb[sim_steps_minus_perturb-1]
        ,marker="s"
        ,linestyle="--"
        ,label="-perturb"
        ,linewidth=line_width
        ,markersize=marker_size
        )
    
    plt.plot(
        r_plus_perturb[sim_steps_unperturb-1] ,z_plus_perturb[sim_steps_unperturb-1]
        ,marker="^"
        ,linestyle="-."
        ,label="+perturb"
        ,linewidth=line_width
        ,markersize=marker_size
        )
    
    plt.plot(
        r_unperturb[sim_steps_unperturb-1],z_unperturb[sim_steps_unperturb-1]
        ,marker="*"
        ,linestyle="dotted"
        ,label="unperturb"
        ,linewidth=line_width
        ,markersize=marker_size
        )
    
    
    xmax = max(r_unperturb_init) #+ ds
    xmin =  min([
        min(r_unperturb_init)
        ,min(r_minus_perturb[sim_steps_minus_perturb-1])
        ,min(r_plus_perturb[sim_steps_unperturb-1])
        ,min(r_unperturb[sim_steps_unperturb-1])
    ])
    ymin =  min([
        min(z_unperturb_init)
        ,min(z_minus_perturb[sim_steps_minus_perturb-1])
        ,min(z_plus_perturb[sim_steps_unperturb-1])
        ,min(z_unperturb[sim_steps_unperturb-1])
    ])
    ymax = xmax - xmin 
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.xlabel(r"r [$\mu m$]",fontsize=15)
    plt.ylabel(r"z [$\mu m$]",fontsize=15)
    plt.legend(fontsize=15)
    plt.grid()
    plt.show()


def plot_test_3d(data_path:str,df_name:str,output_path:str
        ):
    if data_path == "" or df_name== "" or output_path=="":
        print(f" No paths were given, in the plot_Epot_Ekin function")
        exit()
    save_name = df_name.replace(".avi","")

    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    #exit()
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    r_unperturbed = df_sim['r unperturbed'][0]
    z_unperturbed = df_sim['z unperturbed'][0]
    psi = df_sim['psi'][0]
    Area = df_sim['area list'][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    k = df_sim['k'][0]
    kG = df_sim['kG'][0]
    ds = df_sim['ds'][0]
    sigma = df_sim['sigma'][0]
    tau = df_sim['tau'][0]
    dt = df_sim["dt"][0]
    sim_steps = df_sim["sim_steps"][0]
    S = df_sim['Epot'][0]
    T = df_sim['Ekin'][0]
    """---------------------------------Show different positons a 5 different t---------------------------------------------------------------"""
    #fig, ax = plt.subplots()
    #font_size= 15
    #wm = plt.get_current_fig_manager()
    #wm.window.state('zoomed')
    #ax = plt.figure().add_subplot(projection='3d')
    t_ref_list = [0, int(sim_steps/4), int(2*sim_steps/4), int(3*sim_steps/4), int(sim_steps)-1 ]
    style_list =["solid","dotted","dashed","dashdot","solid"]
    color_list = ["k","g","r","b","m","c"]
    marker_list= ["o","v","s","^","*",]


    Plot_3D_mesh(
        N=N ,r=r,z=z
        ,time_pos_list=t_ref_list
        ,linestyle_list=style_list
        ,marker_list=marker_list
        ,color_list=color_list
    )


def plot_comparison_of_plus_minus_un_perturbed_results(path):
    from pathlib import Path
    path = path + "..\\"  # To go one folder back, as the system works from finding the data files and not the folders.
    directory_list = list()
    root_list = list()
    for root, dirs, files in os.walk(path, topdown=False):
        for df_name in files:
            if ".pkl" in df_name:
                data_path = root + "\\"
                directory_list.append(data_path+ df_name)#files[0])
                root_list.append(root)

    ref_folder_list = ["+perturbed","-perturbed","unperturbed"]
    place_ref = ["cross sims","triangle sims","plus sims"]
    #print(Path(root_list[0]).parts)
    
    do_plots = False
    for i in range(len(place_ref)):
        for j in range(len(root_list)):
            A = Path(root_list[j]).parts
            if place_ref[i] in A:
                B = Path(path).parts
                if ref_folder_list[0] in B:
                    #print(f"B={B}")
                    #print("do plots = true")
                    do_plots = True
    
    #print(f"do_plots ={do_plots}")
    if do_plots == True:
        for i in range(len(directory_list)):
            folders = Path(directory_list[i]).parent.parts
            if folders[len(folders)-1] == ref_folder_list[0]:
                df_sim_plus_perturb = pd.read_pickle(directory_list[i])
                #print(directory_list[i])

            if folders[len(folders)-1] == ref_folder_list[1]:
                df_sim_minus_perturb =pd.read_pickle(directory_list[i])
                #print(directory_list[i])

            if folders[len(folders)-1] == ref_folder_list[2]:
                df_sim_unperturb = pd.read_pickle(directory_list[i])
                #print(directory_list[i])


        """----------------------------------------- Minus perturbed data-----------------------------------------------------"""
        X2_minus_perturb = df_sim_minus_perturb["Chi squared test"][0]
        sim_steps_minus_perturb = df_sim_minus_perturb["sim_steps"][0]
        dt_minus_perturb = df_sim_minus_perturb["dt"][0]
        r_minus_perturb = df_sim_minus_perturb["r"][0]
        z_minus_perturb = df_sim_minus_perturb["z"][0]
        S_minus_perturb = df_sim_minus_perturb["Epot"][0]
        time_minus_perturb = np.linspace(0,sim_steps_minus_perturb*dt_minus_perturb,sim_steps_minus_perturb-1)

        color_minus_perturb = "orange"
        linestyle_minus_perturb = "--"
        marker_minus_perturb = "s"
        label_minus_perturb = "-perturbed"

        """----------------------------------------- Plus perturbed data-----------------------------------------------------"""
        X2_plus_perturb = df_sim_plus_perturb["Chi squared test"][0]
        sim_steps_plus_perturb = df_sim_plus_perturb["sim_steps"][0]
        dt_plus_perturb = df_sim_plus_perturb["dt"][0]
        r_plus_perturb = df_sim_plus_perturb["r"][0]
        z_plus_perturb = df_sim_plus_perturb["z"][0]
        S_plus_perturb = df_sim_plus_perturb["Epot"][0]
        time_plus_perturb = np.linspace(0,sim_steps_plus_perturb*dt_plus_perturb,sim_steps_plus_perturb-1)

        color_plus_perturb = "green"
        linestyle_plus_perturb = "dashdot"
        marker_plus_perturb = "^"
        label_plus_perturb = "+perturbed"

        """----------------------------------------- Unperturbed data-----------------------------------------------------"""
        X2_unperturb = df_sim_unperturb["Chi squared test"][0]
        sim_steps_unperturb = df_sim_unperturb["sim_steps"][0]
        dt_unperturb = df_sim_unperturb["dt"][0]
        r_unperturb = df_sim_unperturb["r"][0]
        z_unperturb = df_sim_unperturb["z"][0]
        S_unperturb = df_sim_unperturb["Epot"][0]

        r_unperturb_init = df_sim_unperturb["r unperturbed"][0]
        z_unperturb_init = df_sim_unperturb["z unperturbed"][0]

        time_unperturb = np.linspace(0,sim_steps_unperturb*dt_unperturb,sim_steps_unperturb-1)
        ds = df_sim_unperturb["ds"][0]
        sigma = df_sim_unperturb["sigma"][0]
        tau = df_sim_unperturb["tau"][0]
        N = df_sim_unperturb["N"][0]

        color_unperturb = "red"
        linestyle_unperturb = "dotted"
        marker_unperturb = "*"
        label_unperturb = "unperturbed"

        color_unperturb_init = "blue"
        linestyle_unperturb_init = "-"
        marker_unperturb_init = "o"
        label_unperturb_init = "unperturbed initial position"

        """----------------------------------------- Compare all X^2 results -----------------------------------------------------"""
        fig, ax = plt.subplots()
        wm = plt.get_current_fig_manager()
        wm.window.state('zoomed')
        
        line_width = 2
        plt.plot(
            time_minus_perturb,X2_minus_perturb
            ,label= label_minus_perturb
            ,linewidth=line_width
            ,linestyle=linestyle_minus_perturb
            ,color=color_minus_perturb
            )
        plt.plot(
            time_plus_perturb,X2_plus_perturb
            ,label= label_plus_perturb
            ,linewidth=line_width
            ,linestyle=linestyle_plus_perturb
            ,color=color_plus_perturb
            )
        plt.plot(
            time_unperturb,X2_unperturb
            ,label= label_unperturb
            ,linewidth=line_width
            ,linestyle=linestyle_unperturb
            ,color=color_unperturb
            )

        ymax = max([
            max(X2_minus_perturb)
            ,max(X2_plus_perturb)
            ,max(X2_unperturb)
            ])*1.05
        
        #xmax = sim_steps_unperturb*dt_unperturb*1.01 
        #xmin = min([-dt_minus_perturb,-dt_plus_perturb,-dt_unperturb])/1e3
        ymin = min([
            min(X2_minus_perturb)
            ,min(X2_plus_perturb)
            ,min(X2_unperturb)
            ])*(1-0.05)

        #plt.xlim(xmin, xmax)
        #plt.ylim(ymin,ymax)
        plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
        plt.legend(fontsize=20)
        plt.xlabel("t [s]",fontsize=15)
        plt.ylabel(r"$\chi^2$ [$\mu m^2$]",fontsize=15)
        plt.title(
            r"$\chi^2$ test with respect to the inital position of the unperturbed configuration"
            +f"\n"
            #+r"$\sigma \approx $" +f"{sigma/1000:0.1f}" + r" $nN/\mu m$ and $\tau \approx$" +f"{tau/1000:0.1f} nN"
            +r"$\sigma \approx $" +f"{round(sigma/1000,3)}" + r" $nN/\mu m$ and $\tau \approx$" +f"{round(tau/1000,3)} nN"
            ,fontsize=15)
        plt.grid()
        plt.draw()
        plt.pause(2)

        save_name_1 = "X^2 to unperturb init state"
        plt.savefig(path + save_name_1 +".png")
        plt.savefig(path + save_name_1 +".svg")
        
        
        """----------------------------------------- chi^2 test for the end of unperturbed with the others -----------------------------------------------------"""    
        fig, ax = plt.subplots()
        wm = plt.get_current_fig_manager()
        wm.window.state('zoomed')
        X2_minus_perturbed_2 = np.zeros(sim_steps_unperturb-1)
        X2_plus_perturbed_2 = np.zeros(sim_steps_plus_perturb-1)

        for t in range(sim_steps_minus_perturb-1):
            X2_minus_perturbed_2[t] = Xsqaured_test(
                N=N
                ,r_init=r_unperturb[sim_steps_unperturb-1]
                ,z_init=z_unperturb[sim_steps_unperturb-1]
                ,r=r_minus_perturb[t]
                ,z=z_minus_perturb[t]
                )
            
        for t in range(sim_steps_plus_perturb-1):
            X2_plus_perturbed_2[t] = Xsqaured_test(
                N=N
                ,r_init=r_unperturb[sim_steps_unperturb-1]
                ,z_init=z_unperturb[sim_steps_unperturb-1]
                ,r=r_plus_perturb[t]
                ,z=z_plus_perturb[t]
                )

        plt.plot(
            time_minus_perturb ,X2_minus_perturbed_2
            ,label=label_minus_perturb
            ,linestyle=linestyle_minus_perturb
            ,color=color_minus_perturb
            )
        plt.plot(
            time_plus_perturb ,X2_plus_perturbed_2
            ,label=label_plus_perturb
            ,linestyle=linestyle_plus_perturb
            ,color=color_plus_perturb
            )

        plt.hlines(
            y=0,xmin=-1,xmax=1
            ,label="target"
            ,linestyles="dashed"
            ,color = "black"
            )

        xmax = max([max(time_minus_perturb), max(time_plus_perturb),max(time_unperturb)])
        plt.xlim(xmax*(-5e-3) ,xmax*1.05)
        plt.ylim()
        plt.grid()
        plt.legend(fontsize=20)
        plt.xlabel("t [s]",fontsize=15)
        plt.ylabel(r"$\chi^2$ [$\mu m^2$]",fontsize=15)
        plt.title(
            r"$\chi^2$ with respect to the final position of the dynamical simulation of the unperturbed state"
            +f"\n"
            #+r"$\sigma \approx $" +f"{sigma/1000:0.1f}" + r" $nN/\mu m$ and $\tau \approx$" +f"{tau/1000:0.1f} nN"
            +r"$\sigma \approx $" +f"{round(sigma/1000,3)}" + r" $nN/\mu m$ and $\tau \approx$" +f"{round(tau/1000,3)} nN"
            ,fontsize=15
            )
        plt.draw()
        plt.pause(2)

        save_name_2 = "X^2 to dynamic unperturbed final state"
        plt.savefig(path + save_name_2 +".png")
        plt.savefig(path + save_name_2 +".svg")
        
        
        """----------------------------------------- Compare all final positions results -----------------------------------------------------"""    
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        fig, ax = plt.subplots()
        ax.set_aspect("equal",adjustable="box")
        wm = plt.get_current_fig_manager()
        wm.window.state('zoomed')
        sub_ax = inset_axes(
        parent_axes=ax,
        width=4,
        height=4,
        loc='upper right',
        axes_kwargs={
            'facecolor':"#C9C9C9"
        })

        
        line_width = 2.5
        marker_size = 10
        ax.plot(
            r_unperturb_init,z_unperturb_init
            ,marker=marker_unperturb_init
            ,linestyle=linestyle_unperturb_init
            ,color=color_unperturb_init
            ,label=label_unperturb_init
            ,linewidth=line_width
            ,markersize=marker_size
            )
        sub_ax.plot(
            r_unperturb_init,z_unperturb_init
            ,marker=marker_unperturb_init
            ,linestyle=linestyle_unperturb_init
            ,color=color_unperturb_init
            ,label=label_unperturb_init
            ,linewidth=line_width
            ,markersize=marker_size
        )
        ax.plot(
            r_minus_perturb[sim_steps_minus_perturb-1],z_minus_perturb[sim_steps_minus_perturb-1]
            ,marker=marker_minus_perturb
            ,linestyle=linestyle_minus_perturb
            ,color=color_minus_perturb
            ,label=label_minus_perturb
            ,linewidth=line_width
            ,markersize=marker_size
            )
        
        sub_ax.plot(
            r_minus_perturb[sim_steps_minus_perturb-1],z_minus_perturb[sim_steps_minus_perturb-1]
            ,marker=marker_minus_perturb
            ,linestyle=linestyle_minus_perturb
            ,color=color_minus_perturb
            ,label=label_minus_perturb
            ,linewidth=line_width
            ,markersize=marker_size
            )
        
        ax.plot(
            r_plus_perturb[sim_steps_unperturb-1] ,z_plus_perturb[sim_steps_unperturb-1]
            ,marker=marker_plus_perturb
            ,linestyle=linestyle_plus_perturb
            ,color=color_plus_perturb
            ,label=label_plus_perturb
            ,linewidth=line_width
            ,markersize=marker_size
            )
        
        sub_ax.plot(
            r_plus_perturb[sim_steps_unperturb-1] ,z_plus_perturb[sim_steps_unperturb-1]
            ,marker=marker_plus_perturb
            ,linestyle=linestyle_plus_perturb
            ,color=color_plus_perturb
            ,label=label_plus_perturb
            ,linewidth=line_width
            ,markersize=marker_size
            )
        
        ax.plot(
            r_unperturb[sim_steps_unperturb-1],z_unperturb[sim_steps_unperturb-1]
            ,marker=marker_unperturb
            ,linestyle=linestyle_unperturb
            ,color=color_unperturb
            ,label=label_unperturb
            ,linewidth=line_width
            ,markersize=marker_size
        )
        sub_ax.plot(
            r_unperturb[sim_steps_unperturb-1],z_unperturb[sim_steps_unperturb-1]
            ,marker=marker_unperturb
            ,linestyle=linestyle_unperturb
            ,color=color_unperturb
            ,label=label_unperturb
            ,linewidth=line_width
            ,markersize=marker_size
        )
        


        ymax =  max([
            max(z_unperturb_init)
            ,max(z_minus_perturb[sim_steps_minus_perturb-1])
            ,max(z_plus_perturb[sim_steps_unperturb-1])
            ,max(z_unperturb[sim_steps_unperturb-1])
        ]) + ds
        ymin = min([z_minus_perturb[sim_steps_minus_perturb-1][0] , z_plus_perturb[sim_steps_plus_perturb-1][0] , z_unperturb[sim_steps_unperturb-1][0]]) - ds
        r_edge_init = (r_unperturb[0][0]+ r_unperturb[sim_steps_unperturb-1][0] + r_minus_perturb[sim_steps_minus_perturb-1][0] + r_plus_perturb[sim_steps_plus_perturb-1][0])/4 - ds/2
        r_edge_init = 0
        edge_points = 3
        for i in range(edge_points):
            r_edge_init += (r_unperturb[0][i]+ r_unperturb[sim_steps_unperturb-1][i] + r_minus_perturb[sim_steps_minus_perturb-1][i] + r_plus_perturb[sim_steps_plus_perturb-1][i])/(4*edge_points)
        r_edge_init += - ds/2
        xmin = r_edge_init - (ymax-ymin)/2
        xmax =  r_edge_init + (ymax-ymin)/2
        if xmin == xmin and xmax == xmax and ymin==ymin and ymax == ymax:
            sub_ax.set_xlim(xmin, xmax)
            sub_ax.set_ylim(ymin,ymax)
        #sub_ax.set_xlabel(r"r [$\mu m$]",fontsize=15)
        #sub_ax.set_ylabel(r"z [$\mu m$]",fontsize=15)
        #sub_ax.set_xticks([])
        #sub_ax.set_yticks([])
        #sub_ax.legend(fontsize=20)
        sub_ax.grid()
        plt.draw()
        plt.pause(2)

        ax.vlines(ymax=ymax,ymin=ymin,x=xmin,colors="k")
        ax.vlines(ymax=ymax,ymin=ymin,x=xmax,colors="k")
        ax.hlines(y=ymax,xmax=xmax,xmin=xmin,colors="k")
        ax.hlines(y=ymin,xmax=xmax,xmin=xmin,colors="k",label="Zoom in area")

        #ax.plot([xmin,0.1614],[ymax,0.2569],color="k")
        #ax.plot([xmax,0.3185],[ymin,0.0997],color="k")

        xmax = max(r_unperturb_init) #+ ds
        xmin =  min([
            min(r_unperturb_init)
            ,min(r_minus_perturb[sim_steps_minus_perturb-1])
            ,min(r_plus_perturb[sim_steps_unperturb-1])
            ,min(r_unperturb[sim_steps_unperturb-1])
        ])
        ymin =  min([
            min(z_unperturb_init)
            ,min(z_minus_perturb[sim_steps_minus_perturb-1])
            ,min(z_plus_perturb[sim_steps_unperturb-1])
            ,min(z_unperturb[sim_steps_unperturb-1])
        ])
        ymax = xmax - xmin 
        ax.set_xlim(xmin-2*ds,xmax+2*ds)
        ax.set_ylim(ymin-2*ds,ymax+2*ds)
        ax.set_xlabel(r"r [$\mu m$]",fontsize=15)
        ax.set_ylabel(r"z [$\mu m$]",fontsize=15)
        #ax.legend(fontsize=15,loc="lower right")
        ax.legend(
            bbox_to_anchor=(1.02, 1)
            , loc='upper left'
            , borderaxespad=0
            ,fontsize=15
                    )
        ax.grid()

        plt.title(
            f"The unperturbed intial configuration and the final positions \n"
            +f"of the three different membrane initial conditions \n"
            #+r"$\sigma \approx$" +f"{round(sigma,3)},   " + r"$\tau \approx$"+f"{round(tau,3)}"
            #+r"$\sigma \approx $" +f"{sigma/1000:0.1f}" + r" $nN/\mu m$ and $\tau \approx$" +f"{tau/1000:0.1f} nN"
            +r"$\sigma \approx $" +f"{round(sigma/1000,3)}" + r" $nN/\mu m$ and $\tau \approx$" +f"{round(tau/1000,3)} nN"
            ,fontsize=15
            ,x=0
            ,y=1.01
        )

        plt.pause(1)
        plt.draw()
        save_name_3 = "Compare final pos of membrane positions"
        plt.savefig(path + save_name_3 +".png")
        plt.savefig(path + save_name_3 +".svg")

        
        """----------------------------------------- Compare Potential Energy  -----------------------------------------------------""" 
        fig, ax = plt.subplots()
        wm = plt.get_current_fig_manager()
        wm.window.state('zoomed')
        line_width = 3
        ax.plot(
            time_minus_perturb,S_minus_perturb
            ,linestyle=linestyle_minus_perturb
            ,color=color_minus_perturb
            ,label=label_minus_perturb
            ,linewidth=line_width
        )

        ax.plot(
            time_plus_perturb,S_plus_perturb
            ,linestyle=linestyle_plus_perturb
            ,color=color_plus_perturb
            ,label=label_plus_perturb
            ,linewidth=line_width
        )

        ax.plot(
            time_unperturb,S_unperturb
            ,linestyle=linestyle_unperturb
            ,color=color_unperturb
            ,label=label_unperturb
            ,linewidth=line_width
        )

        plt.legend(fontsize=15)
        plt.title(
            f"The potential energy of the three different initial positions configurations"
            +f"\n"
            +r"$\sigma \approx $" +f"{round(sigma/1000,3)}" + r" $nN/\mu m$ and $\tau \approx$" +f"{round(tau/1000,3)} nN"
            ,fontsize=15
            #,x=0
            ,y=1.01
        )
        plt.xlabel(f"time [s]",fontsize=15)
        plt.ylabel(r"Potential Energy [zJ]",fontsize=15)
        plt.grid()

        plt.pause(1)
        plt.draw()

        save_name_4 = "Compare potential energy"
        plt.savefig(path + save_name_4 +".png")
        plt.savefig(path + save_name_4 +".svg")


def plot_multiprocessing_results():
    path = "2D sim results\\Data for thesis\\multi processor result\\" + "triangle sims\\T,dt,sigma,tau=(1.0e-08,1.1e-13,1.3e+03,2.6e+03)\\"#T,dt,sigma,tau=(1.0e-08,1.0e-13,1.3e+03,2.6e+03)\\"
    path = "2D sim results\\Data for thesis\\" + "fewpoints but low dt\\"+"triangle sims\\T,dt,sigma,tau=(2.0e-08,1.2e-13,1.3e+03,2.6e+03)\\"
    path = "C:\\Users\\adams\\Desktop\\labtop data\\"
    path = "C:\\Users\\adams\\Desktop\\præsentations data\\N=20\\"
    path = "2D sim results\\Data for thesis\\fewpoints but low dt\\triangle sims\\N,T,dt,sigma,tau=(20,2.0e-08,1.0e-13,1.3e+03,2.6e+03)\\"
    path = "2D sim results\\Data for thesis\\really long"
    path = "2D sim results\\Data for thesis\\Verification\\"
    path = "2D sim results\\Data for thesis\\"#new test for N\\"
    #print(path)
    directory_list = list()
    data_files = list()
    make_movie= False#True
    make_figures = False#True
    make_comparison_figs = True
    for root, dirs, files in os.walk(path, topdown=False):
        for df_name in files:
            if ".pkl" in df_name:
                data_path = root + "\\"
                directory_list.append(data_path+ files[0])
                #print(data_path + df_name)
                if make_movie == True:
                    Make_frames(
                        data_path=data_path
                        ,figs_save_path=data_path + "figues for video\\"
                        ,df_name= df_name
                        ,tot_frames= 100
                    )
                    Make_video(
                        output_path=data_path
                        ,input_path=data_path + "figues for video\\"
                        ,video_name= "surface video"
                        ,fps=24
                    )
                if make_figures == True:
                    plot_Epot_Ekin(
                        data_path=data_path
                        ,df_name=df_name
                        ,output_path=data_path
                    )
                    plot_tot_area(
                        data_path=data_path
                        ,df_name=df_name
                        ,output_path=data_path
                    )

                if make_comparison_figs == True:
                    plot_comparison_of_plus_minus_un_perturbed_results(
                        path=data_path
                    )
                #plt.show()
                plt.draw()
                plt.pause(2)
                plt.close("all")


def Front_page_plot():
    from two_d_data_processing import Plot_3D_mesh
    data_path = f"2D sim results\\Data for thesis\\Verification\\c0=c0 tau=0\\"
    df_name = "2D surface N,ds,dt,T,tau,c0=(40, 0.015, 1e-11, 1e-06, 0, 25).pkl"
    output_path = f"2D sim results\\Data for thesis\\"

    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    #exit()
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    N = df_sim["N"][0]
    ds = df_sim['ds'][0]
    dt = df_sim["dt"][0]
    r0 = df_sim["r0"][0]
    sim_steps = df_sim["sim_steps"][0]
    
    r = [r - r0 for r in r]
    r_opposite = [-r for r in r]
    """---------------------------------Show different positons a 5 different t---------------------------------------------------------------"""
    fig, ax = plt.subplots()
    font_size= 15
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    ax.set_aspect("equal",adjustable="box")
    t_ref_list = [int(sim_steps/100), int(10*sim_steps/100), int(40*sim_steps/100), int(70*sim_steps/100), int(sim_steps)-1 ]
    style_list =["solid","dotted","dashed","dashdot","solid"]
    color_list = ["k","g","r","b","m","c"]
    marker_list= ["o","v","s","^","*",]
    for i in range(len(t_ref_list)):
        plt.plot(r[t_ref_list[i]],z[t_ref_list[i]]
        #,marker=marker_list[i]
        ,linestyle="-"#style_list[i]
        ,color="k"#color_list[i]
        #,label=r"t$\approx$"+f"{t_ref_list[i]*dt:0.1e}"
        )
        plt.plot(r_opposite[t_ref_list[i]],z[t_ref_list[i]]
        #,marker=marker_list[i]
        ,linestyle="-"#style_list[i]
        ,color="k"#color_list[i]
        #,label=r"t$\approx$"+f"{t_ref_list[i]*dt:0.1e}"
        )

    #plt.xlabel(r"r [$\mu m$]",fontsize=font_size)
    #lt.ylabel(r"z [$\mu m$]",fontsize=font_size)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    xmin = min([min(r[t]) for t in range(sim_steps)]) - ds
    xmax  = r[0][N] + ds
    ymin = min([min(z[t]) for t in range(sim_steps)]) - ds
    ymax = ymin + (xmax - xmin)


    #plt.xlim(xmin,xmax)
    plt.ylim(ymin,0.2)#ymax)
    ax.set_aspect("equal",adjustable="box")
    #plt.title("show difference from start and end positions",fontsize=15)
    #plt.legend(fontsize=font_size)
    #plt.grid()
    save_name_3 = df_name + "init&end scaled"
    save_name_3 =  "4 different postions in time"
    plt.draw()
    plt.pause(2)
    plt.savefig(output_path + save_name_3 +".png")
    plt.savefig(output_path + save_name_3 +".svg")
    

    """---------------------------------Show different positons a 5 different t 3D plot---------------------------------------------------------------"""
    t_ref_list = [int(10*sim_steps/100), int(40*sim_steps/100) ,sim_steps-1]
    Plot_3D_mesh(
        N=N, r=r ,z=z
        ,linestyle_list = ["-"]
        ,marker_list = ["" for i in range(len(t_ref_list))]
        ,color_list = ["k" for i in range(len(t_ref_list))]
        ,time_pos_list = t_ref_list
    )




if __name__ == "__main__":
    data_path = "C:\\Users\\adams\\Desktop\\praesentations data\\"
    data_path = "2D sim results\\Data for thesis\\Verification\\c0=c0 tau=0\\"
    output_path = data_path
    file_name = "2D Surface.pkl"
    file_name = "2D surface N,ds,dt,T,tau,c0=(40, 0.015, 1e-11, 1e-06, 0, 25)"

    #plot_tot_area()
    #plot_Epot_Ekin()
    #plot_reference_fig_for_finding_what_to_simulate()
    #Find_the_circle_radius_of_rolling_test()
    #Investigating_chosen_configuration_1()
    #figure_3_potential_energy_landscape_cases()
    #Investigating_chosen_configuration_New_data()
    #plot_test_3d(data_path=data_path,df_name=file_name,output_path=output_path)
    #plot_comparison_of_plus_minus_un_perturbed_results(
    #    path="2D sim results\\Data for thesis\\multi processor result\\triangle sims\\N=20\\N,T,dt,sigma,tau=(20,2.0e-08,1.0e-13,1.3e+03,2.6e+03)\\-perturbed\\"
        #path  ="2D sim results\\Data for thesis\\multi processor result\\cross sims\\T,dt,sigma,tau=(1e-07, 1.25e-13,5.8e+03,2.0e+03)\\+perturbed\\"
        #path = "2D sim results\\Data for thesis\\multi processor result\\plus sims\\T,dt,sigma,tau=(1e-07, 1.25e-13,1.494e+04,8.9e+03)\\+perturbed\\"
        #path = "2D sim results\\Data for thesis\Verification\\c0=0 tau=0\\"
    #)
    plot_multiprocessing_results()
    #Front_page_plot()
    #plt.show()