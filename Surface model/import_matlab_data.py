import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cycler
import os

np.set_printoptions(legacy='1.25')

df_name = "Matlab data"
data_path2 = "2D sim results\\"

if not os.path.exists(data_path2 + df_name):
    data_path = "2D sim results\\matlab data transfer files\\" 
    name_list= [
        "sigma_list_Curved"
        ,"tau_list_Curved"
        ,"psi_L_list_Curved"
        ,"ExcessAreaCurved"
        ,"DeltaSACurved"
        ,"NeckRadiusCurved"
        ,"r1Curved"
        ,"sigma_list_Flat"
        ,"tau_list_Flat"
        ,"psi_L_list_Flat"
        ,"ExcessAreaFlat"
        ,"DeltaSAFlat"
        ,"NeckRadiusFlat"
        ,"r1Flat"
        ]

    save_dict = {}

    for name in name_list:

        df_data = pd.read_csv(data_path + name + ".dat")

        names = [i  for i in df_data["index"]]

        placeholderlist = []
        for n in range(len(names)):
            placeholderlist.append([i for i in df_data.iloc[n][1:] if i == i])

        
        save_dict[name] = placeholderlist


    df = pd.DataFrame(save_dict)


    print(data_path2 + df_name)
    if not os.path.exists(data_path2):
        os.makedirs(data_path2)
    df.to_pickle(data_path2 + df_name)

elif os.path.exists(data_path2 + df_name):
    df = pd.read_pickle(data_path2 + df_name)
    print("hello")

print(df.info())
print(len(df["tau_list_Flat"]))
print(len(df["tau_list_Curved"]))
print(df["tau_list_Curved"])


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

for n in range(len(df["tau_list_Flat"])):
    for i in range(len(df["tau_list_Flat"][n])):
        if df["tau_list_Flat"][n][i] < min_tau:
            min_tau = df["tau_list_Flat"][n][i]
        
        if df["tau_list_Flat"][n][i] > max_tau:
            max_tau = df["tau_list_Flat"][n][i]
        
        if df["tau_list_Flat"][n][i] not in diff_taus:
            diff_taus.append(df["tau_list_Flat"][n][i])

print(diff_taus)
print(len(diff_taus))
print(f"mintau={min_tau} and max tau = {max_tau}")

fig, ax = plt.subplots()
cmap = plt.cm.coolwarm
 
index = [i for i in range(len(diff_taus))]
i = 0
plot_tau_ref = []
for n in range(len(df["ExcessAreaCurved"])):
    if len(df["tau_list_Curved"][n]) > 0 :
        i = (df["tau_list_Curved"][n][0]-1)/(max_tau - 1)
        if df["tau_list_Curved"][n][0] not in plot_tau_ref:
            plt.plot(df["ExcessAreaCurved"][n],df["NeckRadiusCurved"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Curved"][n][0]:0.1f}")
            plot_tau_ref.append(df["tau_list_Curved"][n][0])
        else:
            plt.plot(df["ExcessAreaCurved"][n],df["NeckRadiusCurved"][n],".-")#,color=cmap(i))

for n in range(len(df["ExcessAreaFlat"])):
    if len(df["tau_list_Flat"][n]) > 0:
        i = (df["tau_list_Flat"][n][0]-1)/(max_tau -1)
        if df["tau_list_Flat"][n][0] not in plot_tau_ref:
            plt.plot(df["ExcessAreaFlat"][n],df["NeckRadiusFlat"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Flat"][n][0]:0.1f}")
            plot_tau_ref.append(df["tau_list_Flat"][n][0])
        else:
            plt.plot(df["ExcessAreaFlat"][n],df["NeckRadiusFlat"][n],".-")
    

#mpl.rcParams['axes.prop_cycle'] = cycler(color=cmap(np.linspace(min(diff_taus), max(diff_taus), len(diff_taus))))

plt.legend()
plt.vlines(x=0,ymin=-1,ymax=8,colors="k",linestyles="--")
plt.xlim(-60, 50)
plt.ylim(0 ,7)
plt.show()