import numpy as np
import pandas as pd



data_path = "2D sim results\\" 
data_names = ["tau_list_curved.dat","tau list curved.txt"]

data = data_path + data_names[0]

df_data= pd.read_csv(data)


print(
    np.shape(df_data)
)

print(df_data["Unnamed: 1"])#["sigma curved"][:])
