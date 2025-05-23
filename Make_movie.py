import matplotlib.pyplot as plt
import numpy as np
import cv2
from pathlib import Path
import shutil
import os
import pandas as pd
from One_D_Constants import One_D_Constants


def Make_video2(
        output_path
        ,input_path
        ,video_name
        ,fps
        ):
       
    
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




def Make_frames(
        data_path,figs_save_path,df_name
        ):

    df_sim = pd.read_pickle(data_path + df_name)
    
    x = df_sim['x pos'][0]
    z = df_sim['z pos'][0]
    tot_time = df_sim['Total time [sec]']

    plt.figure()
    plt.plot(x[0],z[0])
    plt.show()


if __name__=="__main__":
    args = One_D_Constants()

    L,r0,N,ds,T,dt = args[0:6]
    psi_list,k,c0  =args[6:9]
    save_path, data_path, fig_save_path = args[9:12]
    video_save_path,video_fig_path = args[12:14]
    pass