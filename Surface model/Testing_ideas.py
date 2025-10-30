import numpy as np
from Two_D_functions import Kronecker , c_diff_f,c_diff_g,constraint_f,constraint_g
#import cProfile
#import re



a = 1e-50
alist = np.zeros(10,dtype=float)
for i in range(len(alist)):
    alist[i] = a
da = (alist[2]*alist[1])/0.1
alist[0]= da
print(da)
print(alist)