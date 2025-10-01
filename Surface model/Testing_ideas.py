import numpy as np
from Two_D_functions import Kronecker , c_diff_f,c_diff_g,constraint_f,constraint_g
#import cProfile
#import re




theta = [0, round(np.pi/4,2), round(np.pi/2,2) , round(3*np.pi/4,2), round(np.pi,2)]
print(f"Possible  \n angles = [0, np.pi/4, np.pi/2 ,3*np.pi/4, np.pi] \n angles = {theta}")

point_1 = [0,0]
point_2 = [-1,-1]
x1 = [point_2[0]-point_1[0]]
y1 = [point_2[1]-point_1[1]]

angle = np.arctan2(y1,x1)
head_angle =3*np.pi/4
print(f"angle={angle}")
#print(f"head angle={head_angle}")
#print(f"difference ={head_angle- angle[0]}")

def Get_angle(x1,y1,x2,y2):

    x1 = [x2 - x1]
    y1 = [y2 - y1]

    psi = np.pi - np.arctan2(y1,x1)[0]
    return psi
    
