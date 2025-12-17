import matplotlib.pyplot as plt
import numpy as np
import cv2
from pathlib import Path
import shutil
import os
import pandas as pd
from Two_D_constants import Two_D_Constants, Two_D_paths
import progressbar
import glob 
#from plotting_functions import plot_from_psi_V2

def Make_video(
        output_path
        ,input_path
        ,video_name
        ,fps
        ):
       
    #input_path = input_path +"\\"
    # Create a list of all the input image files
    FILES = []
    num = 0
    for frame in os.listdir(input_path):
        num += 1#FILES.append(frame)
    
    FILES = [f"{i}.png" for i in range(num-1)]
    
    # Get the filename from the output path
    filename = video_name + ".avi"
    #print(f'Creating video "{filename}" from images "{FILES}"')

    # Load the first image to get the frame size
    frame = cv2.imread(input_path + FILES[0],cv2.IMREAD_UNCHANGED)
    

    #print(f" first file name = {FILES[0]}")
    
    
    height, width, layers = np.shape(frame)
    
    # Create a VideoWriter object to write the video file
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')

    video = cv2.VideoWriter(filename=filename
                            , fourcc=fourcc
                            , fps= fps
                            , frameSize=(width, height)
                            )

    # Loop through the input images and add them to the video
    for image_path in FILES:
        #print(f'Adding image "{image_path}" to video "{output_path}"... ')
        video.write(cv2.imread( input_path + image_path))
        plt.pause(0.05)

    

    # Release the VideoWriter and move the output file to the specified location
    cv2.destroyAllWindows()
    video.release() 

    if not os.path.exists(output_path):
            os.makedirs(output_path)
    
    # remove exisiting file
    if os.path.isfile(output_path + filename):
        os.remove(output_path + filename )
    
    
    shutil.move(filename, output_path)

    files = glob.glob(input_path + "/*")
    for f in files:
        os.remove(f)
    print(
        "\n ---------------------- \n"
        +" Movie has been released   ,   name=" + filename
        +"\n ---------------------- \n"
    )




def Make_frames(
        data_path: str
        ,figs_save_path: str
        ,df_name: str
        ,tot_frames: int = 200
        ):

    if not os.path.exists(figs_save_path):
            os.makedirs(figs_save_path)
    
    files = glob.glob(figs_save_path + "/*")
    for f in files:
        os.remove(f)
   
    df_sim = pd.read_pickle(data_path + df_name)
    #print(df_sim.info())
    
    r = df_sim['r'][0]
    z = df_sim['z'][0]
    r_unperturbed = df_sim['r unperturbed'][0]
    z_unperturbed = df_sim['z unperturbed'][0]
    dt = df_sim['dt'][0]
    ds = df_sim['ds'][0]
    L = df_sim["L"][0]
    tau = df_sim["tau"][0]
    sigma = df_sim["sigma"][0]
    N = df_sim["N"][0]
    c0 = df_sim["c0"][0]
    gam2 = df_sim["gam(i>0)"][0]
    sim_steps = df_sim["sim_steps"][0]
    r0 = df_sim["r0"][0]
    #T_tot = df_sim["Total time [sec]"][0]
    T_tot = dt*sim_steps

    textstr = "\n".join((
        f"dt= {dt:0.1e}",
        f"c0={c0:0.1e}",
        f"N={N}",
        f"gam(i>1)={gam2:0.1e}",
        r" $ T_{tot} $ =" + f"{T_tot:0.2e}s",
    ))

    if int(np.shape(r)[0]/tot_frames) < 2:
        frame_vec= [i for i in range(np.shape(r)[0])]
    else:
        frame_vec = [i  for i in range(np.shape(r)[0]) if i%int(np.shape(r)[0]/tot_frames)==0]

    
    fig,ax = plt.subplots()
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #fig.canvas.manager.window.full_screen_toggle()
    #mng = plt.get_current_fig_manager()
    #mng.full_screen_toggle()
    #wm = plt.get_current_fig_manager()
    #wm.window.state('zoomed')

    rmin,rmax = np.min([np.min(i) for i in r]), np.max([np.max(i) for i in r])
    zmin, zmax = np.min([np.min(i) for i in z]) , np.max([np.max(i) for i in z])
    k = 0
    print("\n Making frames progressbar")
    b = progressbar.ProgressBar(maxval=len(frame_vec)-1)
    b.start()
    for t in frame_vec:
        b.update(k)
        plt.plot(r_unperturbed,z_unperturbed,"-s",label="unperturbed state")
        plt.plot(r[t],z[t],'-o',label="dynamics")
        xmax = max(r[0])
        xmin = 0
        #plt.xlim([max(r[0]) - ds*(20), max(r[0])])
        #plt.xlim([xmin,xmax])
        #plt.ylim([-ds, xmax -ds])
        #plt.ylim([zmin*0.99,zmax*1.01])
        plt.xlim(min(r[0])-ds, max(r[0])+ds)
        ceil = max(r[0])-min(r[0]) + 2*ds
        plt.ylim(-ceil/10, 9*ceil/10)
        plt.xlabel(r"r [$\mu m$]")
        plt.ylabel(r"z [$\mu m$]")
        sig_dig = 1/dt
        plt.title(f"Dynamics for time={t*dt:0.2e}s  \n and frame ={k} of {len(frame_vec)}")
        plt.legend()
        plt.text(0.7, 0.95, textstr
                 ,transform=ax.transAxes
                 ,fontsize=12
                 ,verticalalignment='top'
                 ,bbox=props
                 )
        plt.grid()

        plt.draw()
        plt.pause(0.1)

        fig.savefig(figs_save_path +f"{k}")
        fig.clear()
        k += 1
    b.finish()
    plt.close()
    print("\n")



def Get_directory_files(path):
    directory_list = list()
    data_files = list()
    make_movie= True
    make_figures = True
    file_names = []
    folder_path =[]

    for root, dirs, files in os.walk(path, topdown=False):
        for df_name in files:
            if ".pkl" in df_name:
                 file_names.append(df_name)
                 folder_path.append(root + "\\")
        
    return file_names,folder_path

if __name__=="__main__":
    const_args = Two_D_Constants(
        print_val=True
    )

    L,r0,N,ds,T,dt = const_args[0:6]
    k,c0,sim_steps = const_args[6:9]
    sigma, tau, kG = const_args[9:12]
    Area_list, psi_list = const_args[12:14]
    radi_list,z_list = const_args[14:16]

    #sim_steps = 1000

    path_args = Two_D_paths()
    data_path, fig_save_path = path_args[0:2]
    video_save_path,figs_for_video_path = path_args[2:4]
    df_name, fps_movie ,num_frames = path_args[4:7]

    Make_frames(
        data_path=data_path
        ,figs_save_path=figs_for_video_path
        ,df_name=df_name
    )
    Make_video(
        output_path=video_save_path
        ,input_path=figs_for_video_path
        ,video_name= df_name
        ,fps=fps_movie
    )
    
    