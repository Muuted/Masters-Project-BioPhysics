import numpy as np
import matplotlib.pyplot as plt
from One_D_Constants import One_D_Constants,gamma
from One_D_Functions import Lagran_multi, dPsidt_RungeKutta_4
from plotting_functions import plot_from_psi_V2 
import os
import pandas as pd
import progressbar
from Make_movie import Make_frames,Make_video
from data_processing import make_circle, total_energy




def show_radius():
    args = One_D_Constants()
    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    df_name = args[15]
    
    time_vec = np.linspace(
        start=0
        ,stop=T
        ,num=sim_steps-1
    )
    c0_list = [ c0/2 , c0 , c0*2]

    df_name_0 = df_name + f" c0={c0_list[0]}  sim time={T}s"
    df_name_1 = df_name + f" c0={c0_list[1]}  sim time={T}s"
    df_name_2 = df_name + f" c0={c0_list[2]}  sim time={T}s"


    df_sim_0 = pd.read_pickle(data_path + df_name_0)
    df_sim_1 = pd.read_pickle(data_path + df_name_1)
    df_sim_2 = pd.read_pickle(data_path + df_name_2)

    radius_0 = df_sim_0["circle radius"][0]
    radius_1 = df_sim_1["circle radius"][0]
    radius_2 = df_sim_2["circle radius"][0]

    length = len(radius_0)
    fig = plt.figure()
    fig.canvas.manager.window.showMaximized()

    plt.plot(time_vec,radius_0
             ,label=f"rad0  " + f"1/c0={1/c0_list[0]}"
             ,linestyle="-"
             ,marker="*"
             ,color="b"
             )
    plt.hlines(
        y=1/c0_list[0]
        ,xmin=0, xmax=T
        ,linestyles="--",colors="r"
        ,label=f"1/c0={1/c0_list[0]}"
        )
    
    plt.plot(time_vec,radius_1
             ,label= f"rad1  " + f"1/c0={1/c0_list[1]}"
             ,linestyle="-"
             ,marker="*"
             ,color="g"
             )
    plt.hlines(
        y=1/c0_list[1]
        ,xmin=0, xmax=T
        ,linestyles="--",colors="k"
        ,label=f"1/c0={1/c0_list[1]}"
        )
    
    plt.plot(time_vec,radius_2
             ,label="rad2  " + f"1/c0={1/c0_list[2]}"
             ,linestyle="-"
             ,marker="*"
             ,color="c"
             )
    plt.hlines(
        y=1/c0_list[2]
        ,xmin=0, xmax=T
        ,linestyles="--",colors="m"
        ,label=f"1/c0={1/c0_list[2]}"
        )

    plt.xlim(xmin=0, xmax=T+0.1)
    plt.ylim(ymin=0,ymax=radius_0[len(radius_0)-1]*1.5)

    plt.title("time evolution of the radius of the fitted cirles")
    plt.xlabel("Time (s)")
    plt.ylabel("Radius (nm)")
    plt.legend(
        fontsize=15
    )
    plt.show()





def plot_end_result_curve():
    args = One_D_Constants()
    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    df_name = args[15]
    
    c0_list = [ c0/2 , c0 , c0*2]

    df_name_0 = df_name + f" c0={c0_list[0]}  sim time={T}s"
    df_name_1 = df_name + f" c0={c0_list[1]}  sim time={T}s"
    df_name_2 = df_name + f" c0={c0_list[2]}  sim time={T}s"


    df_sim_0 = pd.read_pickle(data_path + df_name_0)
    df_sim_1 = pd.read_pickle(data_path + df_name_1)
    df_sim_2 = pd.read_pickle(data_path + df_name_2)

    print(df_sim_0.info())
    x_pos_0 = df_sim_0["x pos"][0][sim_steps-1]
    z_pos_0 = df_sim_0["z pos"][0][sim_steps-1]
    R_0 = df_sim_0["circle radius"][0][sim_steps-2]
    x_circle_0, z_circle_0 = make_circle(
        xc = df_sim_0["x circle center"][0][sim_steps-2]
        ,zc = df_sim_0["z circle center"][0][sim_steps-2]
        ,R=R_0
        ,xlim=[max(x_pos_0),min(x_pos_0)]
        ,zlim=[max(z_pos_0),min(z_pos_0)]
        , ds=ds
    )

    x_pos_1 = df_sim_1["x pos"][0][sim_steps-1]
    z_pos_1 = df_sim_1["z pos"][0][sim_steps-1]
    R_1 = df_sim_1["circle radius"][0][sim_steps-2]
    x_circle_1, z_circle_1 = make_circle(
        xc = df_sim_1["x circle center"][0][sim_steps-2]
        ,zc = df_sim_1["z circle center"][0][sim_steps-2]
        ,R=R_1
        ,xlim=[max(x_pos_1),min(x_pos_1)]
        ,zlim=[max(z_pos_1),min(z_pos_1)]
        , ds=ds
    )

    x_pos_2 = df_sim_2["x pos"][0][sim_steps-1]
    z_pos_2 = df_sim_2["z pos"][0][sim_steps-1]
    R_2 = df_sim_2["circle radius"][0][sim_steps-2]
    x_circle_2, z_circle_2 = make_circle(
        xc = df_sim_2["x circle center"][0][sim_steps-2]
        ,zc = df_sim_2["z circle center"][0][sim_steps-2]
        ,R=R_2
        ,xlim=[max(x_pos_2),min(x_pos_2)]
        ,zlim=[max(z_pos_2),min(z_pos_2)]
        , ds=ds
    )

    fig = plt.figure()
    fig.canvas.manager.window.showMaximized()

    plt.plot(
        x_pos_0,z_pos_0
        ,label=f"c0={c0_list[0]} " 
            + r", Radius$ \approx $"
            +f"{round(R_0,2)}"
            +f"  , 1/c0= {1/c0_list[0]}"
        ,linestyle="-"
        ,marker="o"
        )
    plt.plot(
        x_pos_1,z_pos_1
        ,label=f"c0={c0_list[1]} "
            + r", Radius$ \approx $"
            +f"{round(R_1,2)}"
            +f"  , 1/c0= {1/c0_list[1]}"
        ,linestyle="-"
        ,marker="o"
        )
    plt.plot(
        x_pos_2,z_pos_2
        ,label=f"c0={c0_list[2]} "
            + r", Radius$ \approx $"
            +f"{round(R_2,2)}"
            +f"  , 1/c0= {1/c0_list[2]}"
        ,linestyle="-"
        ,marker="o"
        )
    plt.plot(
        x_circle_0,z_circle_0,label=f"circle fit, c0={c0_list[0]} "
    )
    plt.plot(
        x_circle_1,z_circle_1,label=f"circle fit, c0={c0_list[1]} "
    )
    plt.plot(
        x_circle_2,z_circle_2,label=f"circle fit, c0={c0_list[2]} "
    )
    plt.xlabel("x",fontsize=25)
    
    plt.ylabel("z",fontsize=25)
    plt.title(
        "Final state of membrane in simulation \n "
        +"for different spontaneous curvatures \n"
        +f"simulation time ={T}s with dt={dt:0.1e}"
        ,fontsize=20
        )

    plt.legend(
        fontsize=20
    )

    """
    fig,ax = plt.subplots(2,2)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    fig.canvas.manager.window.showMaximized()

    ax[0,0].plot(
        x_pos_0,z_pos_0
        ,label=f"c0={c0_list[0]} " 
            + r", Radius$ \approx $"
            +f"{round(R_0,2)}"
            +f"  , 1/c0= {1/c0_list[0]}"
        ,linestyle="-"
        ,marker="o"
        )
    ax[0,0].plot(
        x_circle_0,z_circle_0
        ,label="circle fit"
        )
    ax[0,0].legend(
        fontsize=15
    )


    ax[0,1].plot(
        x_pos_1,z_pos_1
        ,label=f"c0={c0_list[1]} " 
            + r", Radius$ \approx $"
            +f"{round(R_1,2)}"
            +f"  , 1/c0= {1/c0_list[1]}"
        ,linestyle="-"
        ,marker="o"
        )
    ax[0,1].plot(
        x_circle_1,z_circle_1
        ,label="circle fit"
        )
    ax[0,1].legend(
        fontsize=15
    )



    ax[1,0].plot(
        x_pos_2,z_pos_2
        ,label=f"c0={c0_list[2]} " 
            + r", Radius$ \approx $"
            +f"{round(R_2,2)}"
            +f"  , 1/c0= {1/c0_list[2]}"
        ,linestyle="-"
        ,marker="o"
        )
    ax[1,0].plot(
        x_circle_2,z_circle_2
        ,label="circle fit"
        )
    ax[1,0].legend(
        fontsize=15
    )"""

    plt.show()



def plot_energies_multiple():
    args = One_D_Constants()
    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    df_name_orginal = args[15]
    
    time_vec = np.linspace(
        start=0
        ,stop=T
        ,num=sim_steps-1
    )
    c0_list = [ c0/2 , c0 , c0*2]
    m = 1 #grams

    for i in range(len(c0_list)):
        # Get the results from sim 0 ----------------------------------------------------
        df_name = df_name_orginal + f" c0={c0_list[i]}  sim time={T}s"
        df_sim = pd.read_pickle(data_path + df_name)
        psi = df_sim["psi"][0]
        x = df_sim["x pos"][0]
        z = df_sim["z pos"][0]

        E_pot , E_kin ,E_tot = total_energy(
            k=k,c0=c0_list[i]
            ,sigma=k*c0_list[i]**2
            ,ds=ds,m=m
            ,psi=psi,x=x,z=z
        )

        time_vec = np.linspace(
        start=0,stop=T
        ,num=len(E_pot)
        )

        plt.figure()
        plt.plot(
            time_vec, E_tot, label=f"c0={c0_list[i]}"
        )

        plt.title(
            r"The Total energy $E_{tot}$ = T + V" + f"\n"
            +r"$\frac{\Delta E_{tot}}{E_{max}}$="#+"max(E_tot)-min(E_tot)="
            +f"{(max(E_tot)-min(E_tot))/max(E_tot)}"
            )
        plt.xlabel("time [s]")
        plt.ylabel("Energy [units?]")
        plt.legend()

    plt.figure()
    for i in range(len(c0_list)):
        # Get the results from sim 0 ----------------------------------------------------
        df_name = df_name_orginal + f" c0={c0_list[i]}  sim time={T}s"
        df_sim = pd.read_pickle(data_path + df_name)
        psi = df_sim["psi"][0]
        x = df_sim["x pos"][0]
        z = df_sim["z pos"][0]
        
        E_pot , E_kin ,E_tot = total_energy(
            k=k,c0=c0_list[i]
            ,sigma=k*c0_list[i]**2
            ,ds=ds,m=m
            ,psi=psi,x=x,z=z
        )
        time_vec = np.linspace(
        start=0,stop=T
        ,num=len(E_pot)
        )
        plt.plot(
            time_vec, E_tot, label=f"c0={c0_list[i]}"
        )

    plt.title(
        r"The Total energy $E_{tot}$ = T + V"
        )
    plt.xlabel("time [s]")
    plt.ylabel("Energy [units?]")
    plt.legend()
    
    plt.show()


def plot_energies_One():
    args = One_D_Constants()
    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0,sim_steps  =args[6:10]
    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    df_name = args[15]
    
    time_vec = np.linspace(
        start=0
        ,stop=T
        ,num=sim_steps-1
    )
    m = 1 #grams

    fig,ax = plt.subplots()

    # Get the results from sim 0 ----------------------------------------------------
    df_sim = pd.read_pickle(data_path + df_name)
    psi = df_sim["psi"][0]
    x = df_sim["x pos"][0]
    z = df_sim["z pos"][0]
    
    E_pot , E_kin ,E_tot = total_energy(
        k=k,c0=c0
        ,sigma=k*c0**2
        ,ds=ds,m=m
        ,psi=psi,x=x,z=z
    )
    time_vec = np.linspace(
    start=0,stop=T
    ,num=len(E_pot)
    )
    plt.plot(
        time_vec, E_tot, label=f"c0={c0}"
    )

    textstr = "\n".join((
        f"dt= {dt:0.1e}",
        f"ds={ds:0.1e}",
        f"N={N}",
        r" $ T_{tot} $ =" + f"{T}s",
    ))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(0.88, 0.95, textstr
                 ,transform=ax.transAxes
                 , fontsize=14
                 ,verticalalignment='top'
                 , bbox=props
                 )
    
    plt.title(
        r"The Total energy $E_{tot}$ = T + V"
        )
    plt.xlabel("time [s]")
    plt.ylabel("Energy [units?]")
    plt.legend()

    plt.show()



if __name__ == "__main__":
    #show_radius()
    #plot_end_result_curve()
    plot_energies_multiple()
    plot_energies_One()

    