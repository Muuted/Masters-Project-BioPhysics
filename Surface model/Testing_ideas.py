import numpy as np
from Two_D_functions import Kronecker , c_diff_f,c_diff_g,constraint_f,constraint_g
#import cProfile
#import re
import os



path = "2D sim results\\Data for thesis\\multi test\\"


directory_list = list()
data_files = list()
for root, dirs, files in os.walk(path, topdown=False):
    for name in dirs:
        directory_list.append(os.path.join(root, name))
    

for folder in directory_list:
    for root,dirs,files in os.walk(folder,topdown=False):
        
        print(files)