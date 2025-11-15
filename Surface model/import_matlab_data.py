import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.set_printoptions(legacy='1.25')

data_path = "2D sim results\\" 
name_list= [
    "sigma curved","tau curved","psi2 curved" ,"A excess curved" ,"dSA curved" ,"Neck r0 curved" 
    ,"r1 curved" ,"sigma flat" ,"tau flat" ,"psi2 flat" ,"A excess flat" ,"dSA flat" ,"Neck r0 flat" ,"r1 flat"
    ]



save_dict = {}

for name in name_list:

    df_data = pd.read_csv(data_path + name + ".dat")

    names = [i  for i in df_data["index"]]

    print(names)
    placeholderlist = []
    for n in range(len(names)):
        placeholderlist.append([i for i in df_data.iloc[n][1:] if i == i])

    
    save_dict[name] = placeholderlist


df = pd.DataFrame(save_dict)
print(df.info())
exit()
plt.figure()
for n in range(len(names)):
    plt.plot(ExcessAreaCurved[n],NeckRadiusCurved[n],".-",label=f"n={n+1}")

plt.legend()
plt.show()
