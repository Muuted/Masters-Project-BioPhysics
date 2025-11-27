import numpy as np
from Two_D_functions import Kronecker , c_diff_f,c_diff_g,constraint_f,constraint_g
#import cProfile
#import re
import os
from Two_D_constants import Two_D_paths
from Make_movie import Make_frames,Make_video
from two_d_plot_data import plot_Epot_Ekin, plot_tot_area
import matplotlib.pyplot as plt

def plot_everytyhing():
    path = "2D sim results\\Data for thesis\\multi test\\"

    directory_list = list()
    data_files = list()
    make_movie= True
    make_figures = False
    for root, dirs, files in os.walk(path, topdown=False):
        for df_name in files:
            if ".pkl" in df_name:
                
                data_path = root + "\\"
                directory_list.append(data_path+ files[0])
                print(data_path + df_name)
                if make_movie == True:
                    Make_frames(
                        data_path=data_path
                        ,figs_save_path=data_path + "figues for video\\"
                        ,df_name= df_name
                        ,tot_frames= 250
                    )
                    Make_video(
                        output_path=data_path
                        ,input_path=data_path + "figues for video\\"
                        ,video_name= "surface video"
                        ,fps=12
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

                #plt.show()
                plt.draw()
                plt.pause(2)
                plt.close("all")


if __name__ == "__main__":
    plot_everytyhing()