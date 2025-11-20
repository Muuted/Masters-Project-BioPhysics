from Two_D_constants import Two_D_Constants, Two_D_paths
from Two_D_simulation_function import Two_D_simulation
from Make_movie import Make_frames, Make_video
from two_d_data_processing import tot_area, E_pot, E_kin
from Two_D_functions import Langrange_multi, Epsilon_values
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


    
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
    fig, ax=plt.subplots()

    #wm = plt.get_current_fig_manager()
    #wm.window.state('zoomed')
    plt.plot(time,Area_change,'-')
    plt.xlabel("time [s]")
    plt.ylabel("total area")
    plt.title(
        f"Ratio of Amax/Amin={Aratio} \n "
        +f"Amax - AMin={Amax-Amin}"
        )
    ax.ticklabel_format(useOffset=False)
    
    plt.draw()
    plt.pause(0.2)
    save_name_1 = save_name + "Atot"# Total area over time"
    save_name_1 = "Atot"
    plt.savefig(output_path + save_name_1 + ".png")
    plt.pause(0.2)



    fig,ax = plt.subplots()
    plt.plot([i*dt for i in range(len(corr_count))],corr_count,".-")
    plt.title("correction counts pr time")
    plt.xlabel("t[s]")
    plt.ylabel("number of variables corrections")

    plt.draw()
    plt.pause(0.2)
    save_name_1 = save_name + " var corr"# Total area over time"
    save_name_1 =  " var corr"
    plt.savefig(output_path + save_name_1 + ".png")
    plt.pause(0.2)



    dA = np.zeros(sim_steps-1)
    for t in range(sim_steps-1):
        dA[t] = Area_change[t+1] - Area_change[t]


    fig, ax = plt.subplots()
    #wm = plt.get_current_fig_manager()
    #wm.window.state('zoomed')
    plt.plot([i*dt for i in range(len(dA))],dA[0:sim_steps-1],'-')
    #ax.ticklabel_format(useOffset=False)
    plt.title("Change in Area")
    plt.xlabel("t [s]")
    plt.ylabel("Area")
    plt.draw()
    plt.pause(0.3)
    save_name_2 = save_name + " dA"
    save_name_2 = " dA"
    plt.savefig(output_path + save_name_2 + ".png")

    plt.pause(0.3)


    fig,ax = plt.subplots()
    plt.plot(time[0:sim_steps-1],Xsqrt,label=r"$\chi^2$ test")
    plt.title(r"$\chi^2$ test for deviation from the unperturbed state, so $\sigma_i$=1")
    plt.xlabel(" t[s]")
    plt.ylabel(r"$\chi^2$")
    plt.draw()
    plt.pause(0.3)
    save_name_3 = save_name + " chisqrt"
    save_name_3 = " chisqrt"
    plt.savefig(output_path + save_name_3 + ".png")

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
    """T = np.zeros(sim_steps-1)
    S = np.zeros(sim_steps-1)

    for t in range(sim_steps-1):
        T[t ]= (
            E_kin(N=N,t=t,dt=dt,r=r,z=z,Area=Area)
            )
        S[t] = (
            E_pot(N=N,k=k,kG=kG,tau=tau,c0=c0
                  ,r=r[t],psi=psi[t],Area=Area)
        )"""

    t_vec = [dt*i for i in range(sim_steps-1)]
    fontsize = 10
    
    fig, ax = plt.subplots(2,1)
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')

    ax[0].plot(t_vec,T,"-",label="Kinetic energy")
    ax[0].set_xlabel("time [s]",fontsize=fontsize)
    ax[0].set_ylabel(r"$E_{kin}$",fontsize=fontsize)
    ax[0].set_title("Kinetic energy, not scale properly \n" +f"min(E_kin)={min(T)}",fontsize=fontsize)
    ax[0].legend()
    ax[0].ticklabel_format(useOffset=False)
    #plt.figure()
    ax[1].plot(t_vec,S,"-",label="Potential energy")
    ax[1].set_xlabel("time [s]",fontsize=fontsize)
    ax[1].set_ylabel(r"$E_{pot}$",fontsize=fontsize)
    ax[1].set_title("Potential energy  " +r"$min(E_{pot}) \approx$"+f"{round(min(S),3)}  and " +r"$\Delta E_{pot} \approx$"+f"{max(S)-min(S):0.1e}",fontsize=fontsize)
    ax[1].legend()
    ax[1].ticklabel_format(useOffset=False)

    plt.draw()
    plt.pause(0.1)
    save_name_1 = "potential and kinetic energy" + df_name 
    save_name_1 = df_name +" S&T"
    save_name_1 = "S&T"
    plt.savefig(output_path + save_name_1 + ".png")
    
    
    fig, ax = plt.subplots()
    font_size= 15
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.plot(
        r_unperturbed,z_unperturbed
        ,marker="o",linestyle="-"
        ,color="b"
        ,label="Unperturbed initial pos"
        )
    plt.plot(
        r[0],z[0]
        ,marker="o",linestyle="-"
        ,color="k"
        ,label="intial positions"
        )
    if r[sim_steps-1].all() == 0:
        plt.plot(
        r[sim_steps-1],z[sim_steps-1]
        ,marker="o",linestyle="--"
        ,color="g"
        ,label="Error occured, and end isnt there"
        )
    else:
        plt.plot(
        r[sim_steps-1],z[sim_steps-1]
        ,marker="o",linestyle="--"
        ,color="g"
        ,label="end positions"
        )
    plt.xlabel("r",fontsize=font_size)
    plt.ylabel("z",fontsize=font_size)
    plt.title("show difference from start and end positions")
    plt.legend()

    save_name_2 = df_name + "init&end"
    save_name_2 =  "init&end"
    plt.savefig(output_path + save_name_2 +".png")



    fig, ax = plt.subplots()
    font_size= 15
    wm = plt.get_current_fig_manager()
    wm.window.state('zoomed')
    plt.plot(
        r_unperturbed,z_unperturbed
        ,marker="o",linestyle="-"
        ,color="b"
        ,label="Unperturbed initial pos"
        )
    plt.plot(
        r[0],z[0]
        ,marker="o",linestyle="-"
        ,color="k"
        ,label="intial positions"
        )
    if r[sim_steps-1].all() == 0:
        plt.plot(
        r[sim_steps-1],z[sim_steps-1]
        ,marker="o",linestyle="--"
        ,color="g"
        ,label="Error occured, and end isnt there"
        )
    else:
        plt.plot(
        r[sim_steps-1],z[sim_steps-1]
        ,marker="o",linestyle="--"
        ,color="g"
        ,label="end positions"
        )
    plt.xlabel("r",fontsize=font_size)
    plt.ylabel("z",fontsize=font_size)
    xmin,xmax  = max(r[0]) - (N+1)*ds , np.max(r[0]) - (N+1)*ds/2
    ymin,ymax = -ds, xmax - xmin - ds
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    plt.title("show difference from start and end positions")
    plt.legend()

    save_name_3 = df_name + "init&end scaled"
    save_name_3 =  "init&end scaled"
    plt.savefig(output_path + save_name_3 +".png")




def plot_reference_fig_for_finding_what_to_simulate():
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib import cycler
    import os

    np.set_printoptions(legacy='1.25')

    df_name = "Matlab data"
    data_path2 = "2D sim results\\"

    if os.path.exists(data_path2 + df_name):
        df = pd.read_pickle(data_path2 + df_name)
        print("hello")

    print(df.info())
    #print(len(df["tau_list_Flat"]))
    #print(len(df["tau_list_Curved"]))
    #print(df["tau_list_Curved"])


    max_tau = 0
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

    """for n in range(len(df["tau_list_Flat"])):
        for i in range(len(df["tau_list_Flat"][n])):
            if df["tau_list_Flat"][n][i] < min_tau:
                min_tau = df["tau_list_Flat"][n][i]
            
            if df["tau_list_Flat"][n][i] > max_tau:
                max_tau = df["tau_list_Flat"][n][i]
            
            if df["tau_list_Flat"][n][i] not in diff_taus:
                diff_taus.append(df["tau_list_Flat"][n][i])"""

    #print(diff_taus)
    #print(len(diff_taus))
    #print(f"mintau={min_tau} and max tau = {max_tau}")

    fig, ax = plt.subplots()
    cmap = plt.cm.tab20b#cool#coolwarm
    plot_tau_ref = []
    i_ref,n_ref = 0,0
    for n in range(len(df["ExcessAreaCurved"])):
        if len(df["tau_list_Curved"][n]) > 0 :
            i = (df["tau_list_Curved"][n][0]-1)/(max_tau - 1)
            if df["tau_list_Curved"][n][0] not in plot_tau_ref:
                plt.plot(df["ExcessAreaCurved"][n],df["NeckRadiusCurved"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Curved"][n][0]:0.1f}")
                plot_tau_ref.append(df["tau_list_Curved"][n][0])
            else:
                plt.plot(df["ExcessAreaCurved"][n],df["NeckRadiusCurved"][n],".-")#,color=cmap(i))

            for i in range(len(df["ExcessAreaCurved"][n])):
                if 39.0 < df["ExcessAreaCurved"][n][i] < 40 and  1.9 < df["NeckRadiusCurved"][n][i] < 2.1:
                    print(
                        f"sigma = {df["sigma_list_Curved"][n][i]}   "
                        +f"tau = {df["tau_list_Curved"][n][i]}   "
                        +f"sigma = {df["psi_L_list_Curved"][n][i]}   "
                    )
                    plt.plot(df["ExcessAreaCurved"][n][i],df["NeckRadiusCurved"][n][i],"*")#,color=cmap(i))
    """i = 0 
    for n in range(len(df["ExcessAreaFlat"])):
        if len(df["tau_list_Flat"][n]) > 0:
            i = (df["tau_list_Flat"][n][0]-1)/(max_tau -1)
            if df["tau_list_Flat"][n][0] not in plot_tau_ref:
                plt.plot(df["ExcessAreaFlat"][n],df["NeckRadiusFlat"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Flat"][n][0]:0.1f}")
                plot_tau_ref.append(df["tau_list_Flat"][n][0])
            else:
                plt.plot(df["ExcessAreaFlat"][n],df["NeckRadiusFlat"][n],".-")
    """

    font_size = 15
    plt.title("")
    plt.ylabel("Neck Radius",fontsize=font_size)
    plt.xlabel(r"Excess area $\Delta \tilde{A} = \tilde{A}_{neck} - \tilde{A}_{disc}$",fontsize=font_size)
    
    plt.legend(fontsize=font_size-1)
    plt.vlines(x=0,ymin=-1,ymax=8,colors="k",linestyles="--")
    #plt.xlim(-60, 50)
    plt.xlim(-0.1, 50)
    plt.ylim(0 ,7)

    
    plt.show()

if __name__ == "__main__":
    #plot_tot_area()
    #plot_Epot_Ekin()
    plot_reference_fig_for_finding_what_to_simulate()