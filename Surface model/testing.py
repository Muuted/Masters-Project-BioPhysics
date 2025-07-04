import numpy as np
import matplotlib.pyplot as plt


def find_point(data,dx,dy,start_index):
    data_len = np.shape(data) 
    dist_ref = 1e10
    for i in range(data_len):
        if i != start_index:
            if data[start_index][0] - data[i][0] <= dx and  data[start_index][1]-data[i][1] <= dy:
                dist = (data[start_index][0] - data[i][0] )**2  + (data[start_index][1]-data[i][1])**2
            




if __name__ == "__main__":

    x_data = np.linspace(0,10,100)
    y_data = [np.sin(i) for i in x_data]

    plt.figure()
    plt.plot(x_data,y_data)
    plt.show()
