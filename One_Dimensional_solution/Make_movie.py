import matplotlib.pyplot as plt
import numpy as np
import cv2
from pathlib import Path
import shutil
import os
import pandas as pd
from One_D_Constants import One_D_Constants
import progressbar

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
    filename = video_name
    #print(f'Creating video "{filename}" from images "{FILES}"')

    # Load the first image to get the frame size
    frame = cv2.imread(input_path + FILES[0],cv2.IMREAD_UNCHANGED)
    

    print(f" first file name = {FILES[0]}")
    
    
    height, width, layers = np.shape(frame)
    
    # Create a VideoWriter object to write the video file
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')

    video = cv2.VideoWriter(filename=filename, fourcc=fourcc
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

    # remove exisiting file
    if os.path.isfile(output_path + filename):
        os.remove( output_path + filename )
    
    shutil.move(filename, output_path)

    print(
        "\n ---------------------- \n"
        +" Movie has been released"
        +"\n ---------------------- \n"
    )



def Make_frames(
        data_path
        ,figs_save_path
        ,df_name
        ,tot_frames = 216
        ):

    if not os.path.exists(figs_save_path):
            os.makedirs(figs_save_path)

    df_sim = pd.read_pickle(data_path + df_name)
    print(df_sim.info())
    x = df_sim['x pos'][0]
    z = df_sim['z pos'][0]
    dt = df_sim['dt'][0]
    print("dt=",dt)
    tot_time = df_sim['Total time [sec]']
    frame_vec = [i  for i in range(np.shape(x)[0]) if i%int(np.shape(x)[0]/tot_frames)==0]
    
    fig = plt.figure()
    fig.canvas.manager.window.showMaximized()
    xmin,xmax = min([min(i) for i in x]), max([max(i) for i in x])
    zmin, zmax = min([min(i) for i in z]) , max([max(i) for i in z])
    k = 0
    print("\n Plotting progressbar")
    b = progressbar.ProgressBar(maxval=len(frame_vec)-1)
    for t in frame_vec:
        b.update(k)
        plt.plot(x[t],z[t])
        plt.xlim([xmin*0.9, xmax*1.1])
        plt.ylim([zmin*0.9, zmax*1.1])
        plt.xlabel(f"x")
        plt.ylabel(f"z")
        plt.title(f"Dynamics for time={t*dt}s and k={k}")
        plt.pause(0.1)
        plt.draw()
        plt.pause(0.1)

        fig.savefig(figs_save_path +f"{k}")
        fig.clear()
        k += 1
    
    plt.close()

def Make_movie():
    args = One_D_Constants()

    data_path= args[11]
    video_save_path,video_fig_path = args[13:15]
    
    
    
    Make_frames(
        data_path=data_path
        ,figs_save_path=video_fig_path
        ,df_name="1D surface membrane dynamics"
    )
    

    Make_video(
        output_path = video_save_path
        ,input_path = video_fig_path
        ,video_name = "dynamics movie.avi"
        ,fps=12
    )


if __name__=="__main__":
    args = One_D_Constants()

    save_path, data_path, fig_save_path = args[10:13]
    video_save_path,video_fig_path = args[13:15]
    
    making_frame = True
    making_video = True

    if making_frame==True:
        Make_frames(
            data_path=data_path
            ,figs_save_path=video_fig_path
            ,df_name="1D surface membrane dynamics"
        )
    
    if making_video == True:
        Make_video(
            output_path = video_save_path
            ,input_path = video_fig_path
            ,video_name = "dynamics movie.avi"
            ,fps=24
    )
    