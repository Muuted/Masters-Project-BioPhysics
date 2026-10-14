import numpy as np
import matplotlib.pyplot as plt
from Two_D_functions import constraint_f, constraint_g, check_constraints_truth
from error_correction_functions import c_diff,c_diff_f,c_diff_g,Epsilon_v2,Make_variable_corrections




def test_constraint_fg(testcase:dict)->str:
    print("Testing the constraint equations f,g")
    N = testcase["N"]
    r,z,psi = testcase["r"],testcase["z"], testcase["psi"]
    A = testcase["A"]
    
    
    expected_f,expected_g = testcase["expected f"],testcase["expected g"] #These were calculated by hand
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
        return "passed"
    else:
        return "failed"

def test_c_diff_functions(testcase:dict)->str:
    N = testcase["N"]
    r,z,psi = testcase["r"],testcase["z"], testcase["psi"]
    A = testcase["A"]
    return "not finished"


def test_error_corr_function(testcase:dict)-> str:
    print("Testing the error correction function for epsilon,M,b values")
    N = testcase["N"]
    r,z,psi = testcase["r"],testcase["z"], testcase["psi"]
    A = testcase["A"]

    expected_M = testcase["expected M"]
    expected_x = testcase["expected x"]
    expected_epsilon = testcase["expected epsilon"]
    test_passed = True

    epsilon_f,epsilon_g, M, x = Epsilon_v2(
        N=N,r=r,z=z,psi=psi,Area=A,
        testing=True
        )
    epsilon_found = np.append(epsilon_f,epsilon_g)

    row,col = np.shape(M)  # pyright: ignore[reportAssignmentType]

    " testing the matricies "
    if np.shape(M) != np.shape(expected_M) or np.shape(expected_x) != np.shape(x) or np.shape(expected_epsilon) != np.shape(epsilon_found):
        test_passed = False
    
    for i in range(row):
        for j in range(col):
            if expected_M[i][j] != M[i][j]: # pyright: ignore[reportIndexIssue]
                test_passed = False
                print(
                    f"The test failed in the comparison of the matricies for finding the epsilon values \n"
                    f"on index (i,j)={i,j} with the values \n"
                    f"expected_M_{i,j} = {expected_M[i][j]} \n"
                    f"M_{i,j} = {M[i][j]}" # type: ignore
                )

    "Testing the b vectors"
    for i in range(row):
        if expected_epsilon[i] != epsilon_found[i]:
            test_passed = False

        if  expected_x[i] != x[i]: # type: ignore
            test_passed = False


    if test_passed == True:
        return "Passed"
    else:
        return "Failed"

def run_test_constraint_fg():
    "------------------------------ 1st test case -----------------------------------------------------"
    N1 = 2
    r1,z1,psi1 = np.zeros(N1+1,dtype=float),np.zeros(N1+1,dtype=float),np.zeros(N1+1,dtype=float)
    A1 = np.zeros(N1,dtype=float)
    M1 = np.zeros((2*N1,2*N1),dtype=float)
    x1 = np.zeros(2*N1,dtype=float)
    b1 = np.zeros(2*N1,dtype=float)
    epsilon1 = np.zeros(2*N1,dtype=float)

    r1[0],r1[1],r1[2] = 0,1,2
    z1[0] = 1
    psi1[0] = 0.5 
    A1[0],A1[1] = np.pi,3*np.pi
    M1[0][0],M1[0][1] = 4 , 4/3
    M1[1][0],M1[1][1] = 4/4 , 4/6
    M1[2][2],M1[2][3] = 3 , -1
    M1[3][2],M1[3][3] = -1 , 2

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
            "expected M": M1,
            "expected x": x1,
            "expected b": b1,
            "expected epsilon": epsilon1,
            "description": "Testing the constraint equations (f,g) \n "
            "and their differentiated versions (dfdr,dfdz,dfdpsi, dgdr,dgdz,dgdpsi)"
        }
    }


    "------------------------------ 2nd test case -----------------------------------------------------"


    "------------------------------ Running the test cases -----------------------------------------------------"

    print("\n \n Running tests for the error correction function.")
    for i in tests_cases.keys():
        print(f"test description : {tests_cases[i]["description"]}")
        res = test_constraint_fg(testcase=tests_cases[i])
        print(f"The test was : {res} \n ")

        res = test_error_corr_function(testcase=tests_cases[i])
        print(f"The test was : {res} \n ")

        res = test_c_diff_functions(testcase=tests_cases[i])
        print(f"The test was : {res} \n ")

if __name__ == "__main__":
    run_test_constraint_fg()