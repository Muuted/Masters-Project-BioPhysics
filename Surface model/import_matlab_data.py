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
    #print("hello")

print(df.info())
#print(len(df["tau_list_Flat"]))
#print(len(df["tau_list_Curved"]))
#print(df["tau_list_Curved"])


max_tau = -1e20
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


pos_A_min_x1 ,pos_A_min_x2 ,neg_A_min_x3 = 36 ,18.6 ,-12
pos_A_max_x1 ,pos_A_max_x2 ,neg_A_max_x3 = 37 ,19 ,-11
pos_A_min_y1 ,pos_A_min_y2 ,neg_A_min_y3 = 3 ,6.15 ,2.725
pos_A_max_y1 ,pos_A_max_y2 ,neg_A_max_y3 = 3.05 ,6.25 ,2.8

n_neg_A ,n_pos_A = [] ,[]
i_neg ,i_pos = [],[]

for n in range(len(df["ExcessAreaCurved"])):
    for i in range(len(df["ExcessAreaCurved"][n])):
        if len(df["ExcessAreaCurved"]) > 0:
            A_curve = df["ExcessAreaCurved"][n][i]
            r1_curve = df["r1Curved"][n][i]
            if pos_A_min_x1 < A_curve < pos_A_max_x1 and pos_A_min_y1 < r1_curve < pos_A_max_y1:
                n_pos_A.append(n)
                i_pos.append(i)
                print("first")
            if pos_A_min_x2 < A_curve < pos_A_max_x2 and pos_A_min_y2 < r1_curve < pos_A_max_y2:
                n_pos_A.append(n)
                i_pos.append(i)
                print("2nd")

for n in range(len(df["ExcessAreaFlat"])):
    for i in range(len(df["ExcessAreaFlat"][n])):
        if len(df["ExcessAreaFlat"]) > 0:
            A_Flat = df["ExcessAreaFlat"][n][i]
            r1_Flat = df["r1Flat"][n][i]
            if neg_A_min_x3 < A_Flat < neg_A_max_x3 and neg_A_min_y3 < r1_Flat < neg_A_max_y3:
                n_neg_A.append(n)
                i_neg.append(i)

print(f"pos n={n_pos_A} i pos={i_pos}")
print(f"pos n={n_neg_A} i pos={i_neg}")

#print(diff_taus)
#print(len(diff_taus))
print(f"mintau={min_tau} and max tau = {max_tau}")

fig, ax = plt.subplots()
cmap = plt.cm.coolwarm
 
index = [i for i in range(len(diff_taus))]
i = 0
plot_tau_ref = []

plt.plot(
    df["ExcessAreaCurved"][n_pos_A[0]][i_pos[0]]
    ,df["r1Curved"][n_pos_A[0]][i_pos[0]]
    #,label=""
    ,marker="o"
    ,color="k"
    )
plt.plot(
    df["ExcessAreaCurved"][n_pos_A[1]][i_pos[1]]
    ,df["r1Curved"][n_pos_A[1]][i_pos[1]]
    ,marker="o"
    ,color="k"
    )
plt.plot(
    df["ExcessAreaFlat"][n_neg_A[0]][i_neg[0]]
    ,df["r1Flat"][n_neg_A[0]][i_neg[0]]
    ,marker="o"
    ,color="k"
    )

for n in range(len(df["ExcessAreaCurved"])):
    if len(df["tau_list_Curved"][n]) > 0 :
        i = (df["tau_list_Curved"][n][0]-1)/(max_tau - 1)
        if df["tau_list_Curved"][n][0] not in plot_tau_ref:
            #plt.plot(df["ExcessAreaCurved"][n],df["NeckRadiusCurved"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Curved"][n][0]:0.1f}")
            plt.plot(df["ExcessAreaCurved"][n],df["r1Curved"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Curved"][n][0]:0.1f}")
            plot_tau_ref.append(df["tau_list_Curved"][n][0])
        else:
            #plt.plot(df["ExcessAreaCurved"][n],df["NeckRadiusCurved"][n],".-")#,color=cmap(i))
            plt.plot(df["ExcessAreaCurved"][n],df["r1Curved"][n],".-",color=cmap(i))#,color=cmap(i))


for n in range(len(df["ExcessAreaFlat"])):
    if len(df["tau_list_Flat"][n]) > 0:
        i = (df["tau_list_Flat"][n][0]-1)/(max_tau -1)
        if df["tau_list_Flat"][n][0] not in plot_tau_ref:
            #plt.plot(df["ExcessAreaFlat"][n],df["NeckRadiusFlat"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Flat"][n][0]:0.1f}")
            plt.plot(df["ExcessAreaFlat"][n],df["r1Flat"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Flat"][n][0]:0.1f}")
            plot_tau_ref.append(df["tau_list_Flat"][n][0])
        else:
            #plt.plot(df["ExcessAreaFlat"][n],df["NeckRadiusFlat"][n],".-")
            plt.plot(df["ExcessAreaFlat"][n],df["r1Flat"][n],".-",color=cmap(i))
    

print(
    f"sigma ={df["sigma_list_Curved"][n_pos_A[0]][i_pos[0]]}    "
    ,f"tau ={df["tau_list_Curved"][n_pos_A[0]][i_pos[0]]}   "
    ,f"psi2 ={df["psi_L_list_Curved"][n_pos_A[0]][i_pos[0]]}    "
)
print(
    f"sigma ={df["sigma_list_Curved"][n_pos_A[1]][i_pos[1]]}    "
    ,f"tau ={df["tau_list_Curved"][n_pos_A[1]][i_pos[1]]}   "
    ,f"psi2 ={df["psi_L_list_Curved"][n_pos_A[1]][i_pos[1]]}    "
)
print(
    f"sigma ={df["sigma_list_Curved"][n_neg_A[0]][i_neg[0]]}    "
    ,f"tau ={df["tau_list_Curved"][n_neg_A[0]][i_neg[0]]}   "
    ,f"psi2 ={df["psi_L_list_Curved"][n_neg_A[0]][i_neg[0]]:e}    "
)


plt.legend()
plt.vlines(x=0,ymin=-1,ymax=8,colors="k",linestyles="--")
plt.xlim(-60, 50)
plt.ylim(0 ,7)
plt.xlabel(r"$\Delta \tilde{A} = \tilde{A}_{membrane} - \tilde{A}_{disc} $")
plt.ylabel(r"edge radius")
plt.show()