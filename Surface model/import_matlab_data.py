import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cycler
import os

from two_d_continues_integration import find_init_stationary_state
from Two_D_constants import Two_D_Constants_stationary_state


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


pos_A_min_x1 ,pos_A_min_x2 ,neg_A_min_x3 = 36 ,3.4 ,-1.08#-12
pos_A_max_x1 ,pos_A_max_x2 ,neg_A_max_x3 = 37 ,3.6 ,-1.06#-11
pos_r1_min_y1 ,pos_r1_min_y2 ,neg_r1_min_y3 = 3 ,7.35 , 0.64#2.725
pos_r1_max_y1 ,pos_r1_max_y2 ,neg_r1_max_y3 = 3.05 ,7.45 ,0.66#2.8

n_neg_A ,n_pos_A = [] ,[]
i_neg ,i_pos = [],[]

for n in range(len(df["ExcessAreaCurved"])):
    for i in range(len(df["ExcessAreaCurved"][n])):
        if len(df["ExcessAreaCurved"]) > 0:
            A_curve = df["ExcessAreaCurved"][n][i]
            r1_curve = df["r1Curved"][n][i]
            if pos_A_min_x1 < A_curve < pos_A_max_x1 and pos_r1_min_y1 < r1_curve < pos_r1_max_y1:
                n_pos_A.append(n)
                i_pos.append(i)
                print("first")
            if pos_A_min_x2 < A_curve < pos_A_max_x2 and pos_r1_min_y2 < r1_curve < pos_r1_max_y2:
                n_pos_A.append(n)
                i_pos.append(i)
                print("2nd")

for n in range(len(df["ExcessAreaFlat"])):
    for i in range(len(df["ExcessAreaFlat"][n])):
        if len(df["ExcessAreaFlat"]) > 0:
            A_Flat = df["ExcessAreaFlat"][n][i]
            r1_Flat = df["r1Flat"][n][i]
            if neg_A_min_x3 < A_Flat < neg_A_max_x3 and neg_r1_min_y3 < r1_Flat < neg_r1_max_y3:
                n_neg_A.append(n)
                i_neg.append(i)

print(f"pos n={n_pos_A} i pos={i_pos}")
print(f"pos n={n_neg_A} i pos={i_neg}")

#print(diff_taus)
#print(len(diff_taus))
print(f"mintau={min_tau} and max tau = {max_tau}")

"""----------------------------------- Matlab data ----------------------------------------------"""
const_args = Two_D_Constants_stationary_state(
        print_val = False
        ,show_stationary_state = False
        ,start_flat = False
        ,perturb = False
    )

L,r0,N,ds,T,dt = const_args[0:6]
k,c0,sim_steps = const_args[6:9]
sigma, tau, kG = const_args[9:12]

lc = 1/c0
sigma_c = k*c0**2
tau_c = k*c0

rs2 = 20*lc 
zs2 = 0
s0, sN = 0, 50*lc

fig, ax = plt.subplots(1,2)
cmap = plt.cm.coolwarm
wm = plt.get_current_fig_manager()
wm.window.state('zoomed')

index = [i for i in range(len(diff_taus))]
i = 0
plot_tau_ref = []

for n in range(len(df["ExcessAreaCurved"])):
    if len(df["tau_list_Curved"][n]) > 0 :
        i = (df["tau_list_Curved"][n][0]-1)/(max_tau - 1)
        if df["tau_list_Curved"][n][0] not in plot_tau_ref:
            #plt.plot(df["ExcessAreaCurved"][n],df["NeckRadiusCurved"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Curved"][n][0]:0.1f}")
            ax[0].plot(
                [i/lc**2 for i in df["ExcessAreaCurved"][n]]
                ,[i/lc for i in df["r1Curved"][n]],".-",color=cmap(i)
                ,label=r"$\tau \approx$"+f"{df["tau_list_Curved"][n][0]/lc**2:0.1f}"
                )
            plot_tau_ref.append(df["tau_list_Curved"][n][0])
        else:
            #plt.plot(df["ExcessAreaCurved"][n],df["NeckRadiusCurved"][n],".-")#,color=cmap(i))
            ax[0].plot(
                [i/lc**2 for i in df["ExcessAreaCurved"][n]]
                , [i/lc for i in df["r1Curved"][n]]
                ,".-",color=cmap(i)
                )#,color=cmap(i))


for n in range(len(df["ExcessAreaFlat"])):
    if len(df["tau_list_Flat"][n]) > 0:
        i = (df["tau_list_Flat"][n][0]-1)/(max_tau -1)
        if df["tau_list_Flat"][n][0] not in plot_tau_ref:
            #plt.plot(df["ExcessAreaFlat"][n],df["NeckRadiusFlat"][n],".-",color=cmap(i),label=r"$\tau$"+f"={df["tau_list_Flat"][n][0]:0.1f}")
            ax[0].plot(
                [i/lc**2 for i in df["ExcessAreaFlat"][n]]
                ,[i/lc for i in df["r1Flat"][n]]
                ,".-",color=cmap(i),label=r"$\tau \approx$"+f"{df["tau_list_Flat"][n][0]/lc**2:0.1f}"
                )
            plot_tau_ref.append(df["tau_list_Flat"][n][0])
        else:
            #plt.plot(df["ExcessAreaFlat"][n],df["NeckRadiusFlat"][n],".-")
            ax[0].plot(
                [i/lc**2 for i in df["ExcessAreaFlat"][n]]
                ,[i/lc for i in df["r1Flat"][n]]
                ,".-",color=cmap(i))
    

print(
    f"   symbol:  ^"
    ,f"sigma ={df["sigma_list_Curved"][n_pos_A[0]][i_pos[0]]}    "
    ,f"tau ={df["tau_list_Curved"][n_pos_A[0]][i_pos[0]]}   "
    ,f"psi2 ={df["psi_L_list_Curved"][n_pos_A[0]][i_pos[0]]}    \n"
    ,f"sigma ={df["sigma_list_Curved"][n_pos_A[0]][i_pos[0]]*sigma_c}    "
    ,f"tau ={df["tau_list_Curved"][n_pos_A[0]][i_pos[0]]*tau_c}   "
)
print("---------------------------------")
print(
    f"    symbol x"
    ,f"sigma ={df["sigma_list_Curved"][n_pos_A[1]][i_pos[1]]}    "
    ,f"tau ={df["tau_list_Curved"][n_pos_A[1]][i_pos[1]]}   "
    ,f"psi2 ={df["psi_L_list_Curved"][n_pos_A[1]][i_pos[1]]}    \n"
    ,f"sigma ={df["sigma_list_Curved"][n_pos_A[1]][i_pos[1]]*sigma_c}    "
    ,f"tau ={df["tau_list_Curved"][n_pos_A[1]][i_pos[1]]*tau_c}   "
)
print("---------------------------------")
print(
    f"    symbol +"
    ,f"sigma ={df["sigma_list_Curved"][n_neg_A[0]][i_neg[0]]}    "
    ,f"tau ={df["tau_list_Curved"][n_neg_A[0]][i_neg[0]]}   "
    ,f"psi2 ={df["psi_L_list_Curved"][n_neg_A[0]][i_neg[0]]:e}    \n"
    ,f"sigma ={df["sigma_list_Curved"][n_neg_A[0]][i_neg[0]]*sigma_c}    "
    ,f"tau ={df["tau_list_Curved"][n_neg_A[0]][i_neg[0]]*tau_c}   "
)

markers = ["^","+","x"]
markers_latex = [r"$\triangle$",r"$\plus$" ,r"$\times$"]
ax[0].plot(
    df["ExcessAreaCurved"][n_pos_A[0]][i_pos[0]]/lc**2
    ,df["r1Curved"][n_pos_A[0]][i_pos[0]]/lc
    ,marker=markers[0] 
    ,color="k"
    ,markersize = 15
    ,mfc = "none"
    )

ax[0].plot(
    df["ExcessAreaCurved"][n_pos_A[1]][i_pos[1]]/lc**2
    ,df["r1Curved"][n_pos_A[1]][i_pos[1]]/lc
    ,marker=markers[1] 
    ,color="k"
    ,markersize = 12
    )
ax[0].plot(
    df["ExcessAreaFlat"][n_neg_A[0]][i_neg[0]]/lc**2
    ,df["r1Flat"][n_neg_A[0]][i_neg[0]]/lc
    ,marker=markers[2] 
    ,color="k"
    ,markersize = 12
    ,mfc = "none"
    )

plt.subplots_adjust(
    #left=0.05
    #,right=0.5
    wspace=0.5
)
#plt.legend(fontsize=13)
ax[0].legend(
    bbox_to_anchor=(1.02, 1)
    , loc='upper left'
    , borderaxespad=0
    ,fontsize=15
             )
ax[0].vlines(x=0,ymin=-100,ymax=300,colors="k",linestyles="--")
#plt.xlim(-60, 50)
ax[0].set_xlim(-1.2e4, 2.9e4)
ax[0].set_ylim(0 ,230)
ax[0].set_title(
    "edge radius vs Excess Area   ,"
    +r"[$\tau]=\frac{\mu g\cdot \mu m }{s^2}$"
    ,fontsize=15)
ax[0].set_xlabel(r"$\Delta A = A_{membrane} - A_{disc} $  [$\mu m^2$]",fontsize=15)
ax[0].set_ylabel(r"edge radius [$\mu m$]",fontsize=15)



"""---------------------------------------------- Initial configurations ---------------------------------------------------"""
#fig,ax = plt.subplots()

print(f"n pos A :{n_pos_A} and i pos = {i_pos}")
print(f"n neg A :{n_neg_A} and i neg = {i_neg}")

sigma_list = [
    df["sigma_list_Curved"][n_pos_A[0]][i_pos[0]]
    ,df["sigma_list_Curved"][n_pos_A[1]][i_pos[1]]
    ,df["sigma_list_Flat"][n_neg_A[0]][i_neg[0]]
    ]
tau_list = [
    df["tau_list_Curved"][n_pos_A[0]][i_pos[0]]
    ,df["tau_list_Curved"][n_pos_A[1]][i_pos[1]]
    ,df["tau_list_Flat"][n_neg_A[0]][i_neg[0]]
    ]
psi_L_list = [
    df["psi_L_list_Curved"][n_pos_A[0]][i_pos[0]]
    ,df["psi_L_list_Curved"][n_pos_A[1]][i_pos[1]]
    ,df["psi_L_list_Flat"][n_neg_A[0]][i_neg[0]]
    ]

i = 0
psi,r,z, r_contin_0, z_contin_0, alpha = find_init_stationary_state(
            k=k ,c0=c0 ,ds=ds, kG=kG
            ,r_L=rs2 ,z_L=zs2 ,s0=s0 ,sN=sN
            ,total_points = N
            ,tau=tau_list[i]*tau_c
            ,sigma=sigma_list[i]*sigma_c
            ,psi_L=psi_L_list[i]
        )
i = 1
psi,r,z, r_contin_1, z_contin_1, alpha = find_init_stationary_state(
            k=k ,c0=c0 ,ds=ds, kG=kG
            ,r_L=rs2 ,z_L=zs2 ,s0=s0 ,sN=sN
            ,total_points = N
            ,tau=tau_list[i]*tau_c
            ,sigma=sigma_list[i]*sigma_c
            ,psi_L=psi_L_list[i]
        )
i = 2
psi,r,z, r_contin_2, z_contin_2, alpha = find_init_stationary_state(
            k=k ,c0=c0 ,ds=ds, kG=kG
            ,r_L=rs2 ,z_L=zs2 ,s0=s0 ,sN=sN
            ,total_points = N
            ,tau=tau_list[i]*tau_c
            ,sigma=sigma_list[i]*sigma_c
            ,psi_L=psi_L_list[i]
        )

spacing = max(z_contin_2) + max(z_contin_1) + max(z_contin_0)

ax[1].plot(
    r_contin_0
    ,z_contin_0 + spacing
    ,label = markers_latex[0] + f",r1={df["r1Curved"][n_pos_A[0]][i_pos[0]]/lc:.1f}"
    ,linestyle="-"
    ,linewidth=2
         )

ax[1].plot(
    r_contin_1
    ,z_contin_1 
    ,label = markers_latex[1] + f",r1={df["r1Curved"][n_pos_A[0]][i_pos[0]]/lc:.1f}"
    ,linestyle="dashed"
    ,linewidth=2
    )

ax[1].plot(
    r_contin_2
    ,z_contin_2 - spacing
    ,label = markers_latex[2]  + f",r1={df["r1Flat"][n_neg_A[0]][i_neg[0]]/lc:.1f}"
    ,linestyle="-."
    ,linewidth=2
    )

rmin = min([min(r_contin_0) , min(r_contin_1) , min(r_contin_2)]) - ds
rmax = max([max(r_contin_0) , max(r_contin_1) , max(r_contin_2)]) + ds

ax[1].set_xlim(rmin,rmax)

ax[1].set_ylim(-(rmax-rmin)/2,(rmax-rmin)/2)
ax[1].legend(fontsize=25)
ax[1].set_xlabel(r"r [$\mu m$]",fontsize=15)
ax[1].set_ylabel(r"z [$\mu m$]",fontsize=15)
ax[1].set_title(
    "Initial configurations of the three chosen \n points in the phase space diagram on the left"
    ,fontsize=15
    )
ax[1].set_yticks([])
#ax.set_aspect("equal",adjustable="box")

plt.show()