import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

np.set_printoptions(legacy='1.25')

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
df_name = "Matlab data"
data_path2 = "2D sim results\\"

print(data_path2 + df_name)
if not os.path.exists(data_path2):
    os.makedirs(data_path2)
df.to_pickle(data_path2 + df_name)

plt.figure()
for n in range(len(df["ExcessAreaCurved"])):
    plt.plot(df["ExcessAreaCurved"][n],df["NeckRadiusCurved"][n],".-")

for n in range(len(df["ExcessAreaFlat"])):
    plt.plot(df["ExcessAreaFlat"][n],df["NeckRadiusFlat"][n],".-")

plt.vlines(x=0,ymin=-1,ymax=8,colors="k",linestyles="--")
plt.xlim(-60, 50)
plt.ylim(0 ,7)
plt.show()
