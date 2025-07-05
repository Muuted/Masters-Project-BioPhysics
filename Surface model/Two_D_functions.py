import numpy as np
from Two_D_constants import Two_D_Constants
#import os
#import pandas as pd
#import progressbar


def Delta_s(A:list,r:list,i:int):
    return A[i]/(np.pi(r[i+1] +r[i]))
