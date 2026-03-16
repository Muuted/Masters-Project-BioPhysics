import os 

def speed_dt_test():
    from Surface_class import Surface_membrane
    import numpy as np
    import os
    import pandas as pd
    import matplotlib.pyplot as plt

    membrane_choice = ["triangle","plus","cross"]
    save_path =  "2D sim results\\object results speed test\\for dt\\"
    #data_name = "2D Surface sim.pkl"
    N_choice= [30,30,30]
    sim_time = 1e-9
    dt_list = np.linspace(1e-12,2.5e-11,10)

    for index in range(len(membrane_choice)):
        for dt in dt_list:
            membrane = Surface_membrane(
                N = N_choice[index]
                ,T = sim_time
                ,dt = dt
                ,const_index = index
                ,save_path = save_path + membrane_choice[index] + "\\"
            )
            #membrane.run_sim()

    data_paths_list = [[],[],[]]
    directory_list = list()
    for root, dirs, files in os.walk(save_path, topdown=False):
        for df_name in files:
            if ".pkl" in df_name:
                data_path = root + "\\"
                directory_list.append(data_path+ files[0])

                if membrane_choice[0] in data_path:
                    data_paths_list[0].append(data_path + df_name)

                if membrane_choice[1] in data_path:
                    data_paths_list[1].append(data_path + df_name)

                if membrane_choice[2] in data_path:
                    data_paths_list[2].append(data_path + df_name)
                
    for i in range(3):
        for _ in range(len(data_paths_list[i])):
            for l in range(len(data_paths_list[i])-1):
                df1 = pd.read_pickle(data_paths_list[i][l])
                dt_df1 = df1["dt"][0]
                df2 = pd.read_pickle(data_paths_list[i][l+1])
                dt_df2 = df2["dt"][0]

                if dt_df2 < dt_df1:
                    data_paths_list[i][l], data_paths_list[i][l+1] = data_paths_list[i][l+1], data_paths_list[i][l]

    dt_val_lists = [[],[],[]]
    sim_time_list = [[],[],[]]
    for i in range(3):
        for j in range(len(data_paths_list[i])):
            df  = pd.read_pickle(data_paths_list[i][j])
            dt_val_lists[i].append(df["dt"][0])
            sim_time_list[i].append(df["simulation time [s]"][0])#/60)

    markers = ["o","s","^"]
    fig, ax = plt.subplots()
    for i in range(3):
        ax.plot(
            dt_val_lists[i], sim_time_list[i]#[:]
            ,label=membrane_choice[i] + f", N={N_choice[i]}"
            ,marker= markers[i]
        )    

    i = 0
    x0 = dt_val_lists[i][0]
    x1 = dt_val_lists[i][len(dt_val_lists[i])-1]
    y0 = sim_time_list[i][0]
    y1 = sim_time_list[i][len(sim_time_list[i])-1]
    tau = -(x1 - x0)/np.log(y1/y0)
    A = y0

    exp_fit = [A*np.exp(-(l- x0)/tau) for l in dt_val_lists[i]]

    ax.plot(
        dt_val_lists[i], exp_fit
        ,marker="."
        ,linestyle="-"
        ,label= r"exp fit y, $\tau$=" + f"{tau:0.2e}"
    )

    plt.legend()
    plt.grid()
    plt.xlabel("dt [s]")
    plt.ylabel("time spend simulating [min]")
    plt.title(f"Real time simulating for T={sim_time:0.1e}"
              f"\n" + r"$ y = A \cdot \exp ({-(x-x_0 )/\tau} )$" 
              ,fontsize=9)
    plt.draw()
    plt.pause(0.5)
    plt.savefig(save_path + "time comparison plot.png")
    plt.show()


def speed_N_test():
    from Surface_class import Surface_membrane
    import numpy as np
    import os
    import pandas as pd
    import matplotlib.pyplot as plt

    membrane_choice = ["triangle","plus","cross"]
    save_path =  "2D sim results\\object results speed test\\for N\\"
    #data_name = "2D Surface sim.pkl"<
    N_choice= [30,20,20]
    sim_time = 1e-9
    dt_list = [1.25e-11,1.25e-11,1.25e-11]#np.linspace(1e-12,2.5e-11,10)

    for index in range(len(membrane_choice)):
        for i in range(21):
            if i < 11 and index == 0:
                membrane = Surface_membrane(
                    N = N_choice[index] + i
                    ,T = sim_time
                    ,dt = dt_list[index]
                    ,const_index = index
                    ,save_path = save_path + membrane_choice[index] + "\\"
                )
                #membrane.run_sim()
            elif index != 0:
                membrane = Surface_membrane(
                    N = N_choice[index] + i
                    ,T = sim_time
                    ,dt = dt_list[index]
                    ,const_index = index
                    ,save_path = save_path + membrane_choice[index] + "\\"
                )
                #membrane.run_sim()
            
    data_paths_list = [[],[],[]]
    directory_list = list()
    for root, dirs, files in os.walk(save_path, topdown=False):
        for df_name in files:
            if ".pkl" in df_name:
                data_path = root + "\\"
                directory_list.append(data_path+ files[0])

                if membrane_choice[0] in data_path:
                    data_paths_list[0].append(data_path + df_name)

                if membrane_choice[1] in data_path:
                    data_paths_list[1].append(data_path + df_name)

                if membrane_choice[2] in data_path:
                    data_paths_list[2].append(data_path + df_name)
                
    for i in range(3):
        for _ in range(len(data_paths_list[i])):
            for l in range(len(data_paths_list[i])-1):
                df1 = pd.read_pickle(data_paths_list[i][l])
                N_df1 = df1["N"][0]
                df2 = pd.read_pickle(data_paths_list[i][l+1])
                N_df2 = df2["N"][0]

                if N_df2 < N_df1:
                    data_paths_list[i][l], data_paths_list[i][l+1] = data_paths_list[i][l+1], data_paths_list[i][l]

    N_val_lists = [[],[],[]]
    sim_time_list = [[],[],[]]
    for i in range(3):
        for j in range(len(data_paths_list[i])):
            df  = pd.read_pickle(data_paths_list[i][j])
            N_val_lists[i].append(df["N"][0])
            sim_time_list[i].append(df["simulation time [s]"][0]/60)

    markers = ["o","s","^"]
    fig, ax = plt.subplots()
    for i in range(3):
        ax.plot(
            N_val_lists[i], sim_time_list[i][:]
            ,label=membrane_choice[i] 
            ,marker= markers[i]
        )
        a,b = np.polyfit(N_val_lists[i], sim_time_list[i][:],deg=1)
        linreq = [a*n + b for n in N_val_lists[2]]
        ax.plot(N_val_lists[2],linreq
                ,label=f"lin fit: " + membrane_choice[i]
                )


    plt.legend()
    plt.grid()
    plt.xlabel("dt [s]")
    plt.ylabel("time spend simulating [min]")
    plt.title(f"Real time simulating for T={sim_time:0.1e} and dt={dt_list[i]:0.1e}")
    plt.draw()
    plt.pause(0.5)
    plt.savefig(save_path + "N comparison plot.png")
    plt.show()



def speed_var_corr_tol_test():
    import numpy as np
    from Surface_class import Surface_membrane
    import pandas as pd
    import matplotlib.pyplot as plt
    from two_d_data_processing import get_files

    tot_sim_time = 1e-8
    tol_list = np.linspace(1e-5,1e-3,20)
    #tol_list = np.linspace(1e-5,1e-1,10)
    init_config = ["triangle","plus","cross"]
    save_path = f"2D sim results\\object results speed test\\tol test\\T={tot_sim_time}\\"

    for tol in tol_list:
        for i in range(3):
            membrane = Surface_membrane(
                N = 30
                ,T = tot_sim_time
                ,dt = 1e-11
                ,const_index = i
                ,save_path =save_path + init_config[i] + "\\" + f"tol={tol:0.3e}\\"
            )
            membrane.var_corr_tol = tol
            membrane.show_stationary_state = False
            membrane.print_constants = False

            membrane.run_sim()
    
    dict_list = get_files(save_path)
    data_lists = [[],[],[]]

    for data in dict_list:
        if "triangle" in data:
            data_lists[0].append(data)
        if "plus" in data:
            data_lists[1].append(data)
        if "cross" in data:
            data_lists[2].append(data)
 
    for i in range(3):
        for _ in range(len(data_lists[i])):
            for j in range(len(data_lists[i])-1):
                df1 = pd.read_pickle(data_lists[i][j])
                df2 = pd.read_pickle(data_lists[i][j+1])

                tol1 , tol2 = df1["Tolerance"][0] , df2["Tolerance"][0]

                if tol1 > tol2 :
                    data_lists[i][j], data_lists[i][j+1] = data_lists[i][j+1], data_lists[i][j]


    sim_time_list = [[],[],[]]
    for i in range(3):
        for j in range(len(data_lists[i])):
            df = pd.read_pickle(data_lists[i][j])
            sim_time = df["simulation time [s]"][0]/60 # make it into minutes
            sim_time_list[i].append(sim_time)

    fig, ax = plt.subplots()
    for i in range(3):
        ax.plot(
            tol_list, sim_time_list[i]
            ,marker="."
            ,linestyle="-"
            ,label=init_config[i]
        )

    plt.xlabel("Tolerance for when to correct variables")
    plt.ylabel("sim time [min]")
    plt.title(f"compare total sim time given T={tot_sim_time:0.1e}")
    plt.ticklabel_format(style="sci",scilimits=[0,1])
    plt.grid()
    plt.legend()
    plt.draw()
    plt.pause(0.5)
    plt.savefig(save_path + "compare sim time with Tolerances.png")
    plt.show()
    

def Overflow_for_dt_test(
    tot_sim_time = 1e-8
    ,N_val = [20,20,20]
    ,plot_fig = True
    ):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from Surface_class import Surface_membrane
    from two_d_data_processing import get_files

    #tot_sim_time = 1e-8
    #N_val = [20,20,20]
    #N_val = [30,30,30]
    dt_val_list = np.linspace(1e-11,1e-9,10)
    init_config = ["triangle","plus","cross"]
    save_path = f"2D sim results\\object results speed test\\dt overflow test\\T={tot_sim_time} N={N_val[1]}\\"

    for dt in dt_val_list:
        for i in range(3):
            membrane = Surface_membrane(
                N = N_val[i]
                ,T = tot_sim_time
                ,dt = dt
                ,const_index = i
                ,save_path =save_path + init_config[i] + "\\" + f"dt={dt:0.2e}\\"
            )
            membrane.show_stationary_state = False
            membrane.print_constants = False

            membrane.run_sim()

    dict_list = get_files(save_path)
    dict = [[],[],[]]
    for d in range(len(dict_list)):
        file = dict_list[d]
        if init_config[0] in file:
            dict[0].append(file)
        if init_config[1] in file:
            dict[1].append(file)
        if init_config[2] in file:
            dict[2].append(file)

    for i in range(len(init_config)):
        length = len(dict[i])
        for _ in range(length):
            for j in range(length-1):
                df1 = pd.read_pickle(dict[i][j])
                dt1 = df1["dt"][0]
                df2 = pd.read_pickle(dict[i][j+1])
                dt2 = df2["dt"][0]

                if dt1 > dt2:
                    dict[i][j], dict[i][j+1] = dict[i][j+1] , dict[i][j]    

    data_sim_time, data_dt_vals = [[],[],[]], [[],[],[]]
    for i in range(3):
        for j in range(len(dt_val_list)):
            df = pd.read_pickle(dict[i][j])
            if df["Overflow err"][0] == False:
                data_sim_time[i].append( df["simulation time [s]"][0] )
                data_dt_vals[i].append( df["dt"][0] )


    fig, ax = plt.subplots()
    for i in range(3):
        ax.plot(
            data_dt_vals[i], data_sim_time[i]
            ,marker=".",linestyle="-"
            ,label= init_config[i] + f", N={N_val[i]}"
        )
    ref_frame = [0 for i in range(len(dt_val_list))]
    ax.plot(dt_val_list,ref_frame
            ,marker="|",linestyle="-"
            ,label="possible val of dt")    
    
    plt.xlabel("dt [s]")
    plt.ylabel("sim time [s]")
    plt.title(
        "Testing Overflow error from dt."
        +f"\n T={tot_sim_time}"
        )
    plt.grid()
    plt.legend()
    
    plt.draw()
    plt.pause(0.5)
    plt.savefig(save_path + "Overflow test graph" +".png")

    if plot_fig == True:    
        plt.show()
    




def compare_thesis_data_and_new_data():
    import numpy as np
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    pass

if __name__ == "__main__":
    #speed_dt_test()
    #speed_N_test()
    #compare_thesis_data_and_new_data()
    #speed_var_corr_tol_test()
    Overflow_for_dt_test()