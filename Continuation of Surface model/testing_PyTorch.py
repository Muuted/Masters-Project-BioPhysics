import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from scipy.differentiate import derivative


#from Two_D_functions import constraint_f

def constraint_f(r:list,psi:list,i:int,N:int,Area:list[float]) -> float:
    f = ""
    print(i)
    if 0 <= i <= N-1:
        f = np.pi*(r[i+1]**2 - r[i]**2)/Area[i] - np.cos(psi[i])
    else:
        print(f"the value of i is to large, in the constraint equation \n"
              +f"i={i} and N={N}")
        exit()
    if f == "":
        print(f"value of f never took a number")
        exit()
    return f



r = [3,2,1,0]
psi = [0,1,1,0]
Area = [ 1, 1,1,0]
N = 2
i = 1





result = np.gradient(
    f=[
        [constraint_f(r=r
                      ,psi=[psi[0] for _ in range(N)]
                      ,Area=Area,N=N,i=i) for _ in range(N)
                      ],
        [constraint_f(r=[r[0] for _ in range(N)]
                      ,psi=psi
                      ,Area=Area,N=N,i=i) for _ in range(N)
                      ]
    ]
)

print(result)