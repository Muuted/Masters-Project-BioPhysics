
def speed_dt_test():
    from Surface_class import Surface_membrane
    import numpy as np
    import os
    import pandas as pd
    import matplotlib.pyplot as plt

    membrane_choice = ["triangle","plus","cross"]
    save_path =  "2D sim results\\object results speed test\\for dt\\"
    #data_name = "2D Surface sim.pkl"
    N_choice= [30,20,20]
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
            sim_time_list[i].append(df["simulation time [s]"][0]/60)

    markers = ["o","s","^"]
    fig, ax = plt.subplots()
    for i in range(3):
        ax.plot(
            dt_val_lists[i], sim_time_list[i][:]
            ,label=membrane_choice[i] + f", N={N_choice[i]}"
            ,marker= markers[i]
        )
    plt.legend()
    plt.grid()
    plt.xlabel("dt [s]")
    plt.ylabel("time spend simulating [min]")
    plt.title(f"Real time simulating for T={sim_time:0.1e}")
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
    #data_name = "2D Surface sim.pkl"
    N_choice= [30,20,20]
    sim_time = 1e-9
    dt_list = [1.25e-11,1.25e-11,1.25e-11]#np.linspace(1e-12,2.5e-11,10)

    for index in range(len(membrane_choice)):
        for i in range(10):
            membrane = Surface_membrane(
                N = N_choice[index] + i
                ,T = sim_time
                ,dt = dt_list[index]
                ,const_index = index
                ,save_path = save_path + membrane_choice[index] + "\\"
            )
            membrane.run_sim()
            
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
    plt.legend()
    plt.grid()
    plt.xlabel("dt [s]")
    plt.ylabel("time spend simulating [min]")
    plt.title(f"Real time simulating for T={sim_time:0.1e} and dt={dt_list[i]:0.1e}")
    plt.draw()
    plt.pause(0.5)
    plt.savefig(save_path + "time comparison plot.png")
    plt.show()

def compare_thesis_data_and_new_data():
    import numpy as np
    import os
    import pandas as pd
    import matplotlib.pyplot as plt
    pass

if __name__ == "__main__":
    #speed_dt_test()
    speed_N_test()
    #compare_thesis_data_and_new_data()