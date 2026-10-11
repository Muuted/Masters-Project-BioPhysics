import numpy as np
import matplotlib.pyplot as plt
from Two_D_functions import constraint_f, constraint_g, check_constraints_truth
from error_correction_functions import c_diff,c_diff_f,c_diff_g,Epsilon_v2,Make_variable_corrections




def test_constraint_fg(N,r,z,psi,A,expected_f,expected_g):
    """

    N = 2
    r,z,psi = np.zeros((2,N+1),dtype=float),np.zeros((2,N+1),dtype=float),np.zeros((2,N+1),dtype=float)
    A = np.zeros(N,dtype=float)
    
    r[0],r[1],r[2] = 0,1,2
    A[0],A[1] = np.pi,3*np.pi
    
    expected_f,expected_g = [0,0],[0,0] #These were calculated by hand
    """
    results_f, results_g = [],[]

    for i in range(N):
        pass
        results_f.append(constraint_f(i=i,N=N,r=r,psi=psi,Area=A))
        results_g.append(constraint_g(i=i,N=N,r=r,z=z,psi=psi,Area=A))


    f_passes, g_passes = 0,0
    for i in range(N):

        if expected_f[i] == results_f[i]:
            f_passes += 1

        if expected_g[i] == results_g[i]:
            g_passes += 1

    if f_passes == N:
        print(f"The f constraint equations passed")


    if f_passes == N and g_passes == N:
        return "pass"
    else:
        return "fail"

def run_test_constraint_fg():
    N1 = 2
    r1,z1,psi1 = np.zeros(N1,dtype=float),np.zeros(N1,dtype=float),np.zeros(N1,dtype=float)
    A1 = np.zeros(N1,dtype=float)
    
    r1[0],r1[1],r1[2] = 0,1,2
    A1[0],A1[1] = np.pi,3*np.pi
    
    expected_f,expected_g = np.zeros(N1,dtype=float),np.zeros(N1,dtype=float) #These were calculated by hand
    tests_cases = {
        "test 1":{
            "N": N1,
            "r":r1,
            "z":z1,
            "psi":psi1,
            "A":A1,
            "expected f": expected_f,
            "expected g": expected_g,
        },

        "test 2":{}
    }

if __name__ == "__main__":
    pass